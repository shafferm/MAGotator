import argparse
import csv
from collections import defaultdict
import os
import glob
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, FeatureLocation
from Bio.Seq import Seq

def parse_arguments():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(description="Generate GFF and/or GBK files from raw annotations, with specified databases formatting.")
    parser.add_argument("--gff", action='store_true', help="Generate GFF file")
    parser.add_argument("--gbk", action='store_true', help="Generate GBK file")
    parser.add_argument("--samples_paths", nargs='+', help="Alternating list of sample names and paths to their .fna files.")
    parser.add_argument("--database_list", type=str, help="List of databases to include in the annotations. Use 'empty' for all.", default="empty")
    parser.add_argument("--annotations", required=True, help="Path to the raw annotations file")
    args = parser.parse_args()
    args.database_list = None if args.database_list == "empty" else args.database_list.split()
    return args

def parse_samples_and_paths(samples_paths):
    """Parses the alternating list of sample names and .fna file paths into a structured dictionary."""
    # This function already expects samples_paths to be a list
    iterator = iter(samples_paths)
    return dict(zip(iterator, iterator))

def sanitize_description(description):
    """Replace semicolons in descriptions to avoid parsing issues."""
    return description.replace(';', ',')

def format_attributes(annotation, database_list):
    """Format and order database-specific annotations for the GFF attributes column, with customized formatting."""
    attributes = []
    for key, value in sorted(annotation.items()):
        if key.endswith('_id') or key.endswith('_description'):
            db_name = key.split('_')[0]
            if database_list is None or db_name in database_list:
                upper_db_name = db_name.upper()
                if key.endswith('_id'):
                    description_key = f"{db_name}_description"
                    desc = sanitize_description(annotation.get(description_key, "NA"))
                    attributes.append(f"({upper_db_name}) {key} {value}; {desc}")
    return "; ".join(attributes)

def generate_gff(samples_annotations, database_list):
    """Generate GFF files for each sample, filtered by specified databases."""
    for sample, annotations in samples_annotations.items():
        with open(f"GFF/{sample}.gff", "w") as gff_file:
            gff_file.write(f"##gff-version 3\n")
            metadata = annotations[0]  # Assuming shared metadata across each sample's annotations
            gff_file.write(f"# Completeness: {metadata['Completeness']}\n")
            gff_file.write(f"# Contamination: {metadata['Contamination']}\n")
            gff_file.write(f"# Taxonomy: {metadata['taxonomy']}\n")
            
            for annotation in annotations:
                attributes_str = format_attributes(annotation, database_list)
                strand = '+' if annotation['strandedness'] == '+1' else '-'
                gff_line = f"{annotation['query_id']}\t.\tgene\t{annotation['start_position']}\t{annotation['stop_position']}\t.\t{strand}\t.\t{attributes_str}\n"
                gff_file.write(gff_line)

def parse_fna_sequence(fna_file_path):
    """Parse the .fna file to get sequences indexed by their header name."""
    sequences = {}
    for seq_record in SeqIO.parse(fna_file_path, "fasta"):
        # The header name is directly used as the key
        header_name = seq_record.id
        sequences[header_name] = seq_record.seq
    return sequences

def format_qualifiers(annotation, database_list):
    """
    Format database-specific annotations into qualifiers.
    """
    qualifiers = {}
    for key, value in annotation.items():
        if key.endswith('_id') or key.endswith('_description'):
            db_name = key.split('_')[0]
            upper_db_name = db_name.upper()
            if database_list is None or db_name in database_list:
                description_key = f"{db_name}_description"
                desc = annotation.get(description_key, "NA")
                qualifiers[f"({upper_db_name}) {db_name}_id"] = value
                qualifiers[f"({upper_db_name}) description"] = desc
    return qualifiers

def find_sample_fna_files(sample, fna_directory):
    """
    Find all .fna files starting with the sample name in the specified directory.
    """
    pattern = os.path.join(fna_directory, f"{sample}*.fna")
    return glob.glob(pattern)

def aggregate_sample_sequences(sample_files):
    """
    Aggregate sequences from multiple .fna files for a sample.
    """
    sequences = {}
    for file_path in sample_files:
        for seq_record in SeqIO.parse(file_path, "fasta"):
            sequences[seq_record.id] = seq_record.seq
    return sequences

def generate_gbk(samples_annotations, database_list, samples_and_paths):
    """
    Generate GBK files for each sample, containing all annotations for that sample,
    using a direct mapping of samples to .fna file paths.
    """
    print("Starting GBK generation...")
    print(f"All samples and paths: {samples_and_paths}")  # Print the entire samples_paths list

    os.makedirs("GBK", exist_ok=True)

    # Convert the list of tuples into a dictionary for easy lookup
    sample_to_fna_path = dict(samples_and_paths)
    
    for sample, annotations in samples_annotations.items():
        print(f"Processing sample: {sample}")
        
        if sample in sample_to_fna_path:
            fna_file_path = sample_to_fna_path[sample]
            print(f"Using .fna file for {sample}: {fna_file_path}")

            # Check if the .fna file exists
            if os.path.exists(fna_file_path):
                print(f"File exists: {fna_file_path}")
                sequences = parse_fna_sequence(fna_file_path)
                print(f"Aggregated sequences for {sample}: {list(sequences.keys())}")
            else:
                print(f"File does not exist: {fna_file_path}")

            # Continue with the rest of your code for generating GBK...
        else:
            print(f"No .fna file path found for sample {sample}")

def main():
    args = parse_arguments()

    # Load annotations and organize by sample
    samples_annotations = defaultdict(list)
    with open(args.annotations, 'r') as file:
        reader = csv.DictReader(file, delimiter='\t')
        for row in reader:
            samples_annotations[row['sample']].append(row)

    # Directly parse the samples and paths passed as arguments
    samples_and_paths = parse_samples_and_paths(args.samples_paths)

    if args.gbk:
        generate_gbk(samples_annotations, args.database_list, samples_and_paths)

    if args.gff:
        generate_gff(samples_annotations, args.database_list)


if __name__ == "__main__":
    main()
