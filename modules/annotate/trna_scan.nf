process TRNA_SCAN {

    errorStrategy 'finish'

    tag { sample }

    input:
    tuple val(sample), path(fasta)

    output:
    tuple val(sample), path("${sample}_processed_trnas.tsv"), emit: trna_scan_out, optional: true

    script:
    """
    #!/usr/bin/env python

    import pandas as pd
    import subprocess

    # Function to process tRNAscan output
    def process_trnascan_output(input_file, output_file, sample_name):
        # Read the input file into a DataFrame
        trna_frame = pd.read_csv(input_file, sep="\\t", skiprows=[0, 2])

        # Strip leading and trailing spaces from column names
        trna_frame.columns = trna_frame.columns.str.strip()

        # Add a new "sample" column and populate it with the sample_name value
        trna_frame.insert(0, "sample", sample_name)

        # Check if "Note" column is present
        if "Note" in trna_frame.columns:
            # Process the "Note" column to update the "Type" column
            trna_frame["Type"] = trna_frame.apply(lambda row: row["Type"] + " (pseudo)" if str(row["Note"]).lower().startswith("pseudo") else row["Type"], axis=1)

            # Drop the processed "Note" column
            trna_frame = trna_frame.drop(columns=["Note"])

        # Keep only the first occurrence of "Begin" and "End" columns
        trna_frame = trna_frame.loc[:, ~trna_frame.columns.duplicated(keep='first')]

        # Remove columns starting with "Begin" or "End"
        trna_frame = trna_frame.loc[:, ~trna_frame.columns.str.match('(Begin|End)\\.')]

        # Reorder columns
        columns_order = ["sample", "Name", "tRNA #", "Begin", "End", "Type", "Codon", "Score"]

        # Rename specified columns
        trna_frame = trna_frame.rename(columns={"Name": "query_id", "Begin": "begin", "End": "end", "Type": "type", "Codon": "codon", "Score": "score"})

        # Create the "gene_id" column by concatenating "type" and "codon"
        trna_frame["gene_id"] = trna_frame["type"] + " (" + trna_frame["codon"] + ")"

        # Check if DataFrame is empty
        if not trna_frame.empty:
            # Write the processed DataFrame to the output file
            trna_frame.to_csv(output_file, sep="\\t", index=False)

    # Run tRNAscan-SE with the necessary input to avoid prompts
    trna_out = "${sample}_trna_out.txt"
    subprocess.run(["tRNAscan-SE", "-G", "-o", trna_out, "--thread", "${params.threads}", "${fasta}"], input=b'O\\n', check=True)

    # Process tRNAscan-SE output
    #process_trnascan_output(trna_out, "${sample}_processed_trnas.tsv", "${sample}")
    """
}
