import pandas as pd
import os
import glob

# Path to the directory containing BLAST .tbl files
blast_results_dir = "output/blast_results"

# Path to the filtered_insertions.csv file
filtered_insertions_csv = "output/filtered_insertions.csv"

# Read filtered_insertions.csv into a DataFrame, no index column
df_insertions = pd.read_csv(filtered_insertions_csv, index_col=False)

# Create a list of .tbl files
tbl_files = glob.glob(os.path.join(blast_results_dir, "*.tbl"))

# Process each .tbl file
for tbl_file in tbl_files:
    # Extract the transcriptome name from the file name
    transcriptome = os.path.basename(tbl_file).replace(".tbl", "")
    
    # Read the .tbl file into a DataFrame
    df_blast = pd.read_csv(tbl_file, sep="\t", header=None, names=["query_id", "subject_id", "percentage_identity", 
                                                                    "alignment_length", "mismatches", "gap_opens", 
                                                                    "query_start", "query_end", "subject_start", 
                                                                    "subject_end", "evalue", "bit_score"])
    
    # coerce query_id to string
    df_blast.query_id = df_blast.query_id.astype(str)
    
    # Filter the DataFrame to keep only the top bit score for each query_id
    df_blast = df_blast.loc[df_blast.groupby("query_id")["bit_score"].idxmax()]
    
    # Merge the BLAST results with the filtered_insertions DataFrame, adding new columns for bit score, percentage identity, alignment length, and subject_id
    df_insertions = pd.merge(df_insertions, df_blast[["query_id", "bit_score", "percentage_identity", "alignment_length", "subject_id"]], 
                             left_on="insertion_id", right_on="query_id", how="left")
    
    # Rename the columns to include the transcriptome name
    df_insertions = df_insertions.rename(columns={
        "bit_score": f"bit_score_{transcriptome}",
        "percentage_identity": f"percentage_identity_{transcriptome}",
        "alignment_length": f"alignment_length_{transcriptome}",
        "subject_id": f"subject_id_{transcriptome}"
    })
    
    # Drop the extra query_id column
    df_insertions = df_insertions.drop(columns=["query_id"])

# Write the merged DataFrame to a new CSV file
df_insertions.to_csv("output/merged_blast_results.csv", index=False)
