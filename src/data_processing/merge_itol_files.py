#!/usr/bin/env python3
"""
itol_merge_datasets.py

Read several iTOL template/dataset files (e.g. DATASET_COLORSTRIP / simple DATA blocks)
and produce a single iTOL-compliant file that contains one DATASET_COLORSTRIP block per
input file. This follows the iTOL DATASET_COLORSTRIP template and is intended to be
robust to common variations (TAB/COMMA/SPACE separator, optional labels, comments).

Based on iTOL documentation: DATASET_COLORSTRIP template (itol.embl.de/help/dataset_color_strip_template.txt)
and general iTOL template rules.

Usage:
    python itol_merge_datasets.py  out_file dataset1.txt dataset2.txt dataset3.txt

The script will try to preserve dataset labels (DATASET_LABEL) if present in the input
files; otherwise it will use the input filename (basename) as the DATASET_LABEL.

"""
import csv
import sys
from pathlib import Path
from typing import List, Dict, Any


def detect_separator_and_metadata(lines: List[str]) -> Dict[str, Any]:
    """Scan header lines (before DATA) to detect SEPARATOR and optional metadata.

    Returns dict with keys: sep (literal char), sep_token (TAB/COMMA/SPACE), dataset_label (or None), other meta
    """
    meta = {
        "sep_token": None,
        "sep": "\t",  # default to TAB
        "dataset_label": None,
        "color": None,
        "color_branches": None,
    }
    for ln in lines:
        s = ln.strip()
        if not s or s.startswith("#"):
            continue
        parts = s.split()
        # detect SEPARATOR line
        if s.upper().startswith("SEPARATOR"):
            # could be: SEPARATOR TAB  or SEPARATOR \t  or SEPARATOR COMMA
            # try to extract token after SEPARATOR
            tokens = s.split(None, 1)
            if len(tokens) > 1:
                token = tokens[1].strip()
                if token.upper() == "TAB":
                    meta["sep_token"] = "TAB"
                    meta["sep"] = "\t"
                elif token.upper() == "COMMA" or token == ",":
                    meta["sep_token"] = "COMMA"
                    meta["sep"] = ","
                elif token.upper() == "SPACE":
                    meta["sep_token"] = "SPACE"
                    meta["sep"] = " "
                else:
                    # if literal char provided (e.g. \t), try to interpret
                    meta["sep_token"] = token
                    meta["sep"] = token
        # dataset label
        if s.upper().startswith("DATASET_LABEL"):
            # could be: DATASET_LABEL\tlabel or DATASET_LABEL label
            label = s.split(None, 1)[1] if len(s.split(None, 1)) > 1 else None
            if label:
                meta["dataset_label"] = label.strip()
        if s.upper().startswith("COLOR ") or s.upper().startswith("COLOR\t"):
            # capture default dataset color if present
            parts = s.split(None, 1)
            if len(parts) > 1:
                meta["color"] = parts[1].strip()
        if s.upper().startswith("COLOR_BRANCHES"):
            parts = s.split(None, 1)
            if len(parts) > 1:
                meta["color_branches"] = parts[1].strip()
        # stop scanning once we hit DATA (caller will handle)
        if s.upper().startswith("DATA"):
            break
    return meta


def read_itol_data_block(path: Path) -> Dict[str, Any]:
    """Read an iTOL template/dataset file and return a dict with metadata and rows.

    Returned dict keys:
      - dataset_label (str or None)
      - sep (one of '\t', ',', or ' ')
      - rows: list of tuples (id, color, label_or_empty)
    """
    text = path.read_text(encoding="utf-8")
    lines = text.splitlines()

    # split header (before DATA) and data lines (after DATA)
    header_lines = []
    data_lines = []
    in_data = False
    for ln in lines:
        if not in_data and ln.strip().upper().startswith("DATA"):
            in_data = True
            continue
        if not in_data:
            header_lines.append(ln)
        else:
            data_lines.append(ln)

    meta = detect_separator_and_metadata(header_lines)
    sep = meta["sep"]

    # parse data_lines according to sep
    rows = []
    for raw in data_lines:
        if not raw.strip() or raw.lstrip().startswith("#"):
            continue
        # if sep is space we want to split on any whitespace but keep label text intact
        if sep == " ":
            parts = raw.split()
            # id is first, color second if present, the rest is label
            if len(parts) == 1:
                rows.append((parts[0], "", ""))
            elif len(parts) == 2:
                rows.append((parts[0], parts[1], ""))
            else:
                rows.append((parts[0], parts[1], " ".join(parts[2:])))
        else:
            # use csv reader to properly handle commas and tabs
            reader = list(csv.reader([raw], delimiter=sep))
            if not reader or not reader[0]:
                continue
            parts = reader[0]
            # strip surrounding whitespace
            parts = [p.strip() for p in parts]
            if len(parts) == 1:
                rows.append((parts[0], "", ""))
            elif len(parts) == 2:
                rows.append((parts[0], parts[1], ""))
            else:
                # more than 2 columns: treat as id,color,label (label may itself contain separators if originally quoted)
                rows.append((parts[0], parts[1], "\t".join(parts[2:])))

    dataset_label = meta.get("dataset_label") or path.stem
    return {"dataset_label": dataset_label, "sep": sep, "color": meta.get("color"), "color_branches": meta.get("color_branches"), "rows": rows}


def write_combined_colorstrip(out_path: Path, blocks: List[Dict[str, Any]]):
    """Write a multi-dataset iTOL file containing one DATASET_COLORSTRIP block per input block."""
    with out_path.open("w", encoding="utf-8") as out:
        for b in blocks:
            out.write("DATASET_COLORSTRIP\n")
            # always use TAB separator in the combined file for simplicity (safe unless somebody used RGB with commas)
            out.write("SEPARATOR TAB\n")
            out.write(f"DATASET_LABEL\t{b['dataset_label']}\n")
            if b.get("color"):
                out.write(f"COLOR\t{b['color']}\n")
            # preserve color branches setting if present
            if b.get("color_branches") is not None:
                out.write(f"COLOR_BRANCHES\t{b['color_branches']}\n")
            out.write("DATA\n")
            for (id_, color, label) in b["rows"]:
                # id and color are mandatory for color strips; if color missing, write blank
                if label:
                    out.write(f"{id_}\t{color}\t{label}\n")
                else:
                    out.write(f"{id_}\t{color}\n")
            out.write("\n")


def main(argv: List[str]):
    if len(argv) < 3:
        print("Usage: python itol_merge_datasets.py out_file dataset1.txt [dataset2.txt ...]")
        return 1
    out_file = Path(argv[1])
    input_paths = [Path(p) for p in argv[2:]]

    blocks = []
    for p in input_paths:
        if not p.exists():
            raise FileNotFoundError(f"Input file not found: {p}")
        block = read_itol_data_block(p)
        blocks.append(block)

    write_combined_colorstrip(out_file, blocks)
    print(f"Wrote combined file: {out_file}")
    return 0


if __name__ == '__main__':
    raise SystemExit(main(sys.argv))
