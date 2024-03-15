process RRNA_COLLECT {

    errorStrategy 'finish'

    input:
    file combined_rrnas

    output:
    path("collected_rrnas.tsv"), emit: rrna_collected_out, optional: true
    path("combined_rrna_scan.tsv"), emit: rrna_combined_out, optional: true

    script:
    """
    #!/usr/bin/env python

    import os
    import pandas as pd
    from collections import Counter

    # List all tsv files in the current directory
    tsv_files = [f for f in os.listdir('.') if f.endswith('.tsv')]

    # Extract sample names from the file names
    samples = [os.path.basename(file).replace("_processed_rrnas.tsv", "") for file in tsv_files]

    # Initialize flag for checking "NULL" content
    all_files_null = True

    # For storing data from each file
    individual_dfs = []

    # For collected_rrnas.tsv
    collected_data = pd.DataFrame(columns=["gene_id", "gene_description", "category", "topic_ecosystem", "subcategory"] + samples)
    sample_counts = {sample: Counter() for sample in samples}

    for file in tsv_files:
        with open(file, 'r') as f:
            contents = f.read().strip()
            if contents == "NULL":
                continue
            all_files_null = False
            df = pd.read_csv(file, sep='\t')
            # Make sure 'sample' column is added to individual DataFrames
            df['sample'] = os.path.basename(file).replace("_processed_rrnas.tsv", "")
            individual_dfs.append(df)
            
            # Update collected_data and sample_counts for collected_rrnas.tsv
            gene_ids = df['type'].unique()
            for gene_id in gene_ids:
                if gene_id not in collected_data['gene_id'].values:
                    collected_data = collected_data.append({'gene_id': gene_id, 'gene_description': gene_id + " gene", 'category': 'rRNA', 'topic_ecosystem': "", 'subcategory': ""}, ignore_index=True)
                sample_counts[os.path.basename(file).replace("_processed_rrnas.tsv", "")][gene_id] += df[df['type'] == gene_id].shape[0]

    # Check and handle all "NULL" files scenario
    if all_files_null:
        with open("collected_rrnas.tsv", "w") as f:
            f.write("NULL")
        with open("combined_rrna_scan.tsv", "w") as f:
            f.write("NULL")
    else:
        # Finalize collected_rrnas.tsv data
        for sample in samples:
            collected_data[sample] = collected_data['gene_id'].apply(lambda x: sample_counts[sample].get(x, 0))
        collected_data.sort_values(by='gene_id', inplace=True)
        collected_data.drop_duplicates(subset=['gene_id'], inplace=True)
        collected_data.to_csv("collected_rrnas.tsv", sep="\t", index=False)

        # Ensure the order of columns in combined_df includes 'query_id' at the beginning
        # This requires individual_dfs to already have 'query_id' or equivalent. If not, further adjustment is needed
        combined_df = pd.concat(individual_dfs, ignore_index=True)
        # Reorder or ensure 'query_id' is the first column if necessary
        desired_order = ['query_id', 'sample', 'begin', 'end', 'strand', 'type', 'e-value', 'note']
        combined_df = combined_df[[col for col in desired_order if col in combined_df.columns]]  # Adjusting to the desired column order
        combined_df.to_csv("combined_rrna_scan.tsv", sep="\t", index=False)

    """
}
