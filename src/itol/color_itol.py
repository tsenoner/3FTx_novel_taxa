#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on:  Mon 14 Aug 2023 14:53:54
Description: A script to generate iTOL (Interactive Tree Of Life) color files
             from a CSV input. It allows automatic color generation or custom
             color mapping from a JSON file. Now includes legend functionality
             with configurable shapes and positioning.
Usage:       python color_itol.py -i <input_file> --id_column <id_column> --group_column <group_column>
             [--color_file <color_file_path>] [--legend_shape <shape>] [--legend_position <x,y>]

@author: tsenoner
"""

import json
import logging
from argparse import ArgumentParser
from colorsys import rgb_to_hsv
from pathlib import Path
from typing import Dict, List, Optional, Tuple

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
    legend_shape: Optional[int] = None,
    legend_position: Optional[Tuple[int, int]] = None,
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
        # handle.write("TREE_COLORS\n")
        handle.write("DATASET_COLORSTRIP\n")
        handle.write("SEPARATOR TAB\n")
        handle.write(f"DATASET_LABEL\t{group_column}\n")

        # Add legend configuration (always generate legend)
        handle.write(f"LEGEND_TITLE\t{group_column}\n")

        if legend_position:
            handle.write(f"LEGEND_POSITION_X\t{legend_position[0]}\n")
            handle.write(f"LEGEND_POSITION_Y\t{legend_position[1]}\n")

        # Get unique groups for legend
        unique_groups = df[group_column].dropna().unique()
        unique_groups = sorted([g for g in unique_groups if g != "none"])

        # Legend shape (default to square if not provided)
        if legend_shape is None:
            legend_shape = 1  # Default to square
        legend_shapes = [legend_shape] * len(unique_groups)

        # Legend colors
        legend_colors = [color_mapping[group] for group in unique_groups]

        # Legend labels
        legend_labels = [str(group) for group in unique_groups]

        # Write legend configuration
        handle.write(f"LEGEND_SHAPES\t{'\t'.join(map(str, legend_shapes))}\n")
        handle.write(f"LEGEND_COLORS\t{'\t'.join(legend_colors)}\n")
        handle.write(f"LEGEND_LABELS\t{'\t'.join(legend_labels)}\n")

        handle.write("DATA\n")
        # Sort by group name alphabetically before writing
        df_sorted = df.sort_values(by=group_column, na_position="last")
        for _, row in df_sorted.iterrows():
            uid = row[id_column]
            color = row["color"]
            group = row[group_column]
            if group == "none":
                continue
            # handle.write(f"{uid}\trange\t{color}\t{group}\n")
            handle.write(f"{uid}\t{color}\t{group}\n")


def parse_args() -> ArgumentParser:
    """Parse command-line arguments."""
    parser = ArgumentParser(description="Create iTOL color file.")
    parser.add_argument(
        "-i",
        "--input",
        type=Path,
        required=True,
        help="Path to the input CSV file.",
    )
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
    parser.add_argument(
        "--legend_shape",
        type=int,
        choices=[1, 2, 3, 4, 5, 6],
        help=(
            "Shape for all legend entries (optional). "
            "Shapes: 1=square, 2=circle, 3=star, 4=right triangle, 5=left triangle, 6=checkmark. "
            "Example: '2' for circles."
        ),
    )
    parser.add_argument(
        "--legend_position",
        type=str,
        help=(
            "Legend position as 'X,Y' coordinates (optional). "
            "Example: '100,200' for position (100, 200). "
            "If not provided, automatic positioning is used."
        ),
    )
    return parser.parse_args()


def main() -> None:
    """Main function to read input file and create iTOL color file."""

    # Parse arguments and validate input file
    args = parse_args()
    input_file = args.input
    if not input_file.exists():
        logging.error(f"The specified input file '{input_file}' does not exist.")
        return

    # Output file path
    itol_out_path = (
        input_file.parent / f"iTOL_{input_file.stem}_{args.group_column}.txt"
    )

    # Parse legend arguments
    legend_position = None
    if args.legend_position:
        try:
            x, y = args.legend_position.split(",")
            legend_position = (int(x.strip()), int(y.strip()))
        except ValueError:
            logging.error("Invalid legend_position format. Use 'X,Y' coordinates.")
            return

    # Read CSV and write iTOL file
    try:
        df = pd.read_csv(input_file)
        write_itol_file(
            df,
            itol_out_path,
            args.id_column,
            args.group_column,
            args.color_file,
            args.legend_shape,
            legend_position,
        )
        logging.info(
            f"iTOL color file successfully created and saved to: {itol_out_path}"
        )
    except Exception as e:
        logging.error(f"An error occurred: {e}")


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO)
    main()
