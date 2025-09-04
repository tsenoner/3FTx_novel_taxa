#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on:  Mon 14 Aug 2023 14:53:54
Description: A script to generate iTOL (Interactive Tree Of Life) color files
             from a CSV input. It allows automatic color generation or custom
             color mapping from a JSON file.
Usage:       python color_itol.py <input_file> --id_column <id_column> --group_column <group_column> [--color_file <color_file_path>]

@author: tsenoner
"""

import json
import logging
from argparse import ArgumentParser
from colorsys import rgb_to_hsv
from pathlib import Path
from typing import Dict, List

import numpy as np
import pandas as pd
from distinctipy import distinctipy


def generate_colors(num_unique_colors: int, seed: int = 42) -> List[str]:
    """Generate a list of distinct colors sorted by hue."""
    colorblind_type = "Normal"
    colors = distinctipy.get_colors(
        num_unique_colors, colorblind_type=colorblind_type, rng=seed
    )

    # Convert to HEX and sort by hue
    hex_colors = [
        "#" + "".join([f"{int(c * 255):02X}" for c in color]) for color in colors
    ]
    sorted_colors = sorted(
        hex_colors,
        key=lambda color: rgb_to_hsv(
            int(color[1:3], 16) / 255.0,
            int(color[3:5], 16) / 255.0,
            int(color[5:7], 16) / 255.0,
        )[0],
    )
    return sorted_colors


def load_custom_colors(color_file_path: Path) -> Dict[str, str]:
    """Load custom colors from a JSON file."""
    with open(color_file_path, "r") as file:
        return json.load(file)


def write_itol_file(
    df: pd.DataFrame,
    itol_out_path: Path,
    id_column: str,
    group_column: str,
    color_file: Path = None,
) -> None:
    """Write the iTOL color file based on the given DataFrame."""

    # Check for valid columns
    if id_column not in df.columns or group_column not in df.columns:
        raise ValueError(
            f"The specified columns '{id_column}' or '{group_column}' do not"
            " exist in the input file."
        )

    # Generate or load custom colors
    unique_groups = df[group_column].dropna().unique()
    if color_file:
        color_mapping = load_custom_colors(color_file)
    else:
        colors = generate_colors(len(unique_groups))
        color_mapping = dict(zip(sorted(unique_groups), colors))
        if df[group_column].isna().any():
            color_mapping[np.nan] = "#E6E1DB"

    df["color"] = df[group_column].map(color_mapping)

    # Write the iTOL file
    with open(itol_out_path, "w") as handle:
        handle.write("TREE_COLORS\n")
        handle.write("SEPARATOR TAB\n")
        handle.write("DATA\n")
        for _, row in df.iterrows():
            uid = row[id_column]
            color = row["color"]
            group = row[group_column]
            if group == "none":
                continue
            handle.write(f"{uid}\trange\t{color}\t{group}\n")


def parse_args() -> ArgumentParser:
    """Parse command-line arguments."""
    parser = ArgumentParser(description="Create iTOL color file.")
    parser.add_argument("input_file", type=Path, help="Path to the input CSV file.")
    parser.add_argument(
        "-id",
        "--id_column",
        type=str,
        required=True,
        help="Column name for the identifier.",
    )
    parser.add_argument(
        "-g",
        "--group_column",
        type=str,
        required=True,
        help="Column name for the group.",
    )
    parser.add_argument(
        "-c",
        "--color_file",
        type=Path,
        help=(
            "Path to the custom color JSON file (optional). The file should"
            " contain a JSON object with group names as keys and hex colors as"
            " values. Example: {'group1': '#FF0000', 'group2': '#00FF00'}."
        ),
    )
    return parser.parse_args()


def main() -> None:
    """Main function to read input file and create iTOL color file."""

    # Parse arguments and validate input file
    args = parse_args()
    input_file = args.input_file
    if not input_file.exists():
        logging.error(f"The specified input file '{input_file}' does not exist.")
        return

    # Output file path
    itol_out_path = (
        input_file.parent / f"iTOL_{input_file.stem}_{args.group_column}.txt"
    )

    # Read CSV and write iTOL file
    try:
        df = pd.read_csv(input_file)
        write_itol_file(
            df,
            itol_out_path,
            args.id_column,
            args.group_column,
            args.color_file,
        )
        logging.info(
            f"iTOL color file successfully created and saved to: {itol_out_path}"
        )
    except Exception as e:
        logging.error(f"An error occurred: {e}")


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO)
    main()