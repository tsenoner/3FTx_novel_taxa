#!/usr/bin/env python3
"""
Convert aligned FASTA files to PHYLIP-relaxed format.
Extracts clean sequence IDs (part before first space) from FASTA headers.
"""

import argparse
import sys
from pathlib import Path
from Bio import SeqIO, AlignIO
from Bio.Align import MultipleSeqAlignment


def convert_fasta_to_phylip(input_fasta, output_phylip):
    """
    Convert aligned FASTA to PHYLIP-relaxed format.
    
    Args:
        input_fasta (str): Path to input FASTA file
        output_phylip (str): Path to output PHYLIP file
    """
    # Check if input file exists
    if not Path(input_fasta).exists():
        raise FileNotFoundError(f"Input file not found: {input_fasta}")
    
    records = []
    
    try:
        for record in SeqIO.parse(input_fasta, "fasta"):
            # Extrahiere den Teil vor dem ersten Leerzeichen aus der description
            # (description enthält die komplette Zeile ohne '>')
            new_id = record.description.split()[0]  # A0A3M7SVS3_IPR045860_1_133-194
            
            # Setze die neuen Werte
            record.id = new_id
            record.name = new_id  # Setze name auch, für Konsistenz
            record.description = ""  # Leere description für sauberes PHYLIP
            records.append(record)
        
        if not records:
            raise ValueError("No sequences found in input file")
        
        print(f"Processed {len(records)} sequences")
        
        # Überprüfe, ob alle Sequenzen gleich lang sind
        seq_lengths = set(len(rec.seq) for rec in records)
        if len(seq_lengths) > 1:
            print(f"Error: Found sequences with different lengths: {seq_lengths}")
            raise ValueError("Nicht alle Sequenzen sind gleich lang! Alignment ist inkonsistent.")
        
        print(f"All sequences have length: {list(seq_lengths)[0]}")
        
        # Write custom PHYLIP format to ensure single-line entries
        num_sequences = len(records)
        sequence_length = len(records[0].seq)
        
        with open(output_phylip, "w") as handle:
            # Write header line
            handle.write(f"{num_sequences} {sequence_length}\n")
            
            # Write each sequence on a single line
            for record in records:
                handle.write(f"{record.id} {str(record.seq)}\n")
        
        print(f"✅ Fertig! Gespeichert als: {output_phylip}")
        
    except Exception as e:
        print(f"❌ Error: {e}", file=sys.stderr)
        sys.exit(1)


def main():
    parser = argparse.ArgumentParser(
        description="Convert aligned FASTA files to PHYLIP-relaxed format",
        epilog="""
Examples:
  %(prog)s input.fasta                    # Output to input.phy
  %(prog)s input.fasta output.phy         # Specify output file
  %(prog)s aligned.fasta -o result.phylip # Using -o flag
        """,
        formatter_class=argparse.RawDescriptionHelpFormatter
    )
    
    parser.add_argument(
        "input_fasta",
        help="Input FASTA file (aligned sequences)"
    )
    
    parser.add_argument(
        "output_phylip",
        nargs="?",
        help="Output PHYLIP file (optional - defaults to input name with .phy extension)"
    )
    
    parser.add_argument(
        "-o", "--output",
        dest="output_phylip_alt",
        help="Alternative way to specify output file"
    )
    
    parser.add_argument(
        "-v", "--verbose",
        action="store_true",
        help="Verbose output"
    )
    
    args = parser.parse_args()
    
    # Determine output file
    if args.output_phylip_alt:
        output_file = args.output_phylip_alt
    elif args.output_phylip:
        output_file = args.output_phylip
    else:
        # Default: change extension to .phy
        input_path = Path(args.input_fasta)
        output_file = input_path.with_suffix(".phy")
    
    if args.verbose:
        print(f"Input file: {args.input_fasta}")
        print(f"Output file: {output_file}")
        print()
    
    # Convert the file
    convert_fasta_to_phylip(args.input_fasta, str(output_file))


if __name__ == "__main__":
    main()