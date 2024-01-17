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

    def process_trnascan_output(input_file, output_file):
        # Read the input file into a DataFrame
        trna_frame = pd.read_csv(input_file, sep="\t", skiprows=[0, 2])

        # Print column names for debugging
        print("Original Columns:")
        print(trna_frame.columns)

        # Remove the "Note" column if present
        trna_frame = trna_frame.drop(columns=["Note"], errors="ignore")

        # Remove the first and third lines
        trna_frame = trna_frame.iloc[1::2]

        # Remove additional occurrences of "Begin" and "End" columns along with their values
        trna_frame = trna_frame.loc[:, ~trna_frame.columns.duplicated(keep='first')]

        # Reorder columns if "Name" is present
        if "Name" in trna_frame.columns:
            columns_order = ["Name", "tRNA #", "Begin", "End", "Type", "Codon", "Score"]
            trna_frame = trna_frame[columns_order]

        # Write the processed DataFrame to the output file
        trna_frame.to_csv(output_file, sep="\t", index=False)

    # Run tRNAscan-SE
    trna_out = "${sample}_trna_out.txt"
    subprocess.run(["tRNAscan-SE", "-G", "-o", trna_out, "--thread", "${params.threads}", "${fasta}"], check=True)

    # Process tRNAscan-SE output
    process_trnascan_output(trna_out, "${sample}_processed_trnas.tsv")
    """
}
