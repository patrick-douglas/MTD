#!/usr/bin/env python3

import re
import argparse

def extract_clade_names_from_tree(tree_file):
    with open(tree_file, 'r') as file:
        tree_data = file.read()
        
    # Extract clade names (assuming Newick format with parentheses)
    clade_names = set(re.findall(r'\b(\w+)\b', tree_data))
    
    return clade_names

def filter_annotations(annot_file, valid_clades, corrected_file):
    with open(annot_file, 'r') as file:
        lines = file.readlines()
    
    with open(corrected_file, 'w') as file:
        for line in lines:
            parts = line.strip().split()
            if len(parts) >= 3:
                clade = parts[0]
                if clade in valid_clades:
                    file.write(line)
            else:
                # Write lines that don't match the expected format (e.g., configuration lines)
                file.write(line)

def main():
    parser = argparse.ArgumentParser(description="Verify and correct annotation file based on tree.")
    parser.add_argument('tree_file', type=str, help="Path to the tree file.")
    parser.add_argument('annot_file', type=str, help="Path to the annotation file.")
    parser.add_argument('corrected_file', type=str, help="Path to the corrected annotation file.")
    
    args = parser.parse_args()
    
    # Extract clade names from the tree
    valid_clades = extract_clade_names_from_tree(args.tree_file)
    
    # Filter and correct the annotation file
    filter_annotations(args.annot_file, valid_clades, args.corrected_file)
    
    print(f"Annotation file has been corrected and saved to {args.corrected_file}")

if __name__ == "__main__":
    main()

