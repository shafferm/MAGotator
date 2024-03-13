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
    parser.add_argument("--samples_paths_file", required=True, help="Path to the file containing sample names and paths to their .fna files.")
    parser.add_argument("--database_list", type=str, help="List of databases to include in the annotations. Use 'empty' for all.", default="empty")
    parser.add_argument("--annotations", required=True, help="Path to the raw annotations file")
    args = parser.parse_args()
    args.database_list = None if args.database_list == "empty" else args.database_list.split()
    return args

def parse_samples_and_paths(annotation_files):
    """
    Parses the list of sample names and .fna file paths into a structured list of tuples.

    Args:
    - annotation_files: A list of strings where each sample name is followed by its corresponding .fna file path.

    Returns:
    - A list of tuples, each containing a sample name and its .fna file path.
    """
    samples_and_paths = [(annotation_files[i], annotation_files[i + 1]) for i in range(0, len(annotation_files), 2)]
    return samples_and_paths

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
    os.makedirs("GBK", exist_ok=True)

    # Convert the list of tuples into a dictionary for easy lookup
    sample_to_fna_path = dict(samples_and_paths)
    
    for sample, annotations in samples_annotations.items():
        print(f"Processing sample: {sample}")
        
        if sample in sample_to_fna_path:
            fna_file_path = sample_to_fna_path[sample]
            print(f"Using .fna file for {sample}: {fna_file_path}")
            sequences = parse_fna_sequence(fna_file_path)
            print(f"Aggregated sequences for {sample}: {list(sequences.keys())}")
            
            seq_record = SeqRecord(Seq(""), id=sample, description=f"Generated GBK file for {sample}", annotations={"molecule_type": "DNA"})
            
            for annotation in annotations:
                query_id = annotation['query_id']
                print(f"Processing annotation {query_id} for {sample}")
                
                if query_id in sequences:
                    print(f"Match found for {query_id} in {sample}")
                    sequence = sequences[query_id]
                    feature_location = FeatureLocation(start=int(annotation['start_position']) - 1, end=int(annotation['stop_position']), strand=1 if annotation['strandedness'] == '+1' else -1)
                    qualifiers = format_qualifiers(annotation, database_list)
                    feature = SeqFeature(feature_location, type="gene", qualifiers=qualifiers)
                    seq_record.features.append(feature)
                else:
                    print(f"No match found for {query_id} in {sample}")
                    
            if sequences:
                representative_seq = next(iter(sequences.values()))
                seq_record.seq = representative_seq
                print(f"Set representative sequence for {sample}: {repr(representative_seq[:50])}...")
            else:
                print(f"No sequences found for {sample}, GBK file will be empty!")
            
            output_filename = f"GBK/{sample}.gbk"
            with open(output_filename, "w") as output_handle:
                SeqIO.write(seq_record, output_handle, "genbank")
            print(f"GBK file generated for {sample}: {output_filename}")
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

    # Here, you need to ensure the fna_directory is correctly defined or passed
    # This could be another command line argument, or it could be a predefined value within your script
    # Parse sample names and .fna file paths from the specified file
    with open(args.samples_paths_file, 'r') as f:
        samples_and_paths = [line.strip().split() for line in f.readlines()]

    if args.gbk:
        generate_gbk(samples_annotations, args.database_list, samples_and_paths)

    if args.gff:
        generate_gff(samples_annotations, args.database_list)

# Make sure this fna_directory path is correctly passed to the generate_gbk function

if __name__ == "__main__":
    main()
