process MERGE_ANNOTATIONS {

    input:
    path( old_annotations, stageAs: "old_annotations.tsv" )
    path( new_annotations, stageAs: "new_annotations.tsv" )


    output:
    path "raw-annotations.tsv", emit: merged_annots_out

    script:
    """
    #!/usr/bin/env python
    import pandas as pd

    # Define the paths to the existing combined_annotations.tsv file and the user-provided file
    existing_file_path = 'old_annotations.tsv'
    user_file_path = 'new_annotations.tsv'

    # Load the existing combined_annotations.tsv file into a DataFrame
    existing_df = pd.read_csv(existing_file_path, sep='\t')

    # Load the user-provided combined_annotations.tsv file into a DataFrame
    user_df = pd.read_csv(user_file_path, sep='\t')

    # Print column names and data types for diagnostic purposes
    print("Existing DataFrame columns:")
    print(existing_df.columns)
    print("User-provided DataFrame columns:")
    print(user_df.columns)

    # Check if the 'query_id' column exists in both DataFrames
    if 'query_id' not in existing_df.columns:
        raise ValueError("The 'query_id' column does not exist in the existing DataFrame.")
    if 'query_id' not in user_df.columns:
        raise ValueError("The 'query_id' column does not exist in the user-provided DataFrame.")

    # Debug statements to check the contents of the 'query_id' column in each DataFrame
    print("Existing DataFrame 'query_id' values:")
    print(existing_df['query_id'])
    print("User-provided DataFrame 'query_id' values:")
    print(user_df['query_id'])

    # Merge the two DataFrames based on the 'query_id' column
    merged_df = pd.merge(existing_df, user_df, on='query_id', how='outer')

    # Save the merged DataFrame to a new file
    merged_file_path = 'merged_combined_annotations.tsv'
    merged_df.to_csv(merged_file_path, sep='\t', index=False)

    print(f"Merged annotations saved to {merged_file_path}")

    """
}
