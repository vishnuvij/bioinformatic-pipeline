import os
import sys
import glob
import pandas as pd
from Bio import SeqIO

# --- Step 1: Load and Rename the TSV File ---
def load_tsv(tsv_file, output_excel):
    """Read a TSV file, assign column names, and save as an Excel file."""
    column_names = [
        "Protein accession", "Sequence MD5 digest", "Sequence length", "Analysis", 
        "Signature Accession", "Signature description", "Start location", "Stop location",
        "e-value", "Status", "Date", "InterPro annotations- accession", "InterPro annotations- description",
        "GO annotations", "Pathways annotations"
    ]
    
    df = pd.read_csv(tsv_file, sep="\t", header=None, names=column_names)
    df.to_excel(output_excel, index=False)
    print(f"✅ Saved TSV as Excel: {output_excel}")
    return df

# --- Step 2: Filter Data Using the Excel File ---
def filter_protein_accessions(tsv_df, excel_df, output_excel):
    """Match PFAM Code to Signature Accession, extract Protein Accession, and save to Excel."""
    
    # Merge Excel data (PFAM Code) with TSV data (Signature Accession)
    merged_df = excel_df.merge(tsv_df, left_on="PFAM Code", right_on="Signature Accession", how="inner")

    # Select and save relevant columns
    export_df = merged_df[["PFAM Code", "Signature Accession", "Protein accession"]].drop_duplicates()
    export_df.to_excel(output_excel, index=False)
    
    print(f"✅ Saved filtered protein accessions to {output_excel}")

    return export_df["Protein accession"].astype(str).tolist()

# --- Step 3: Process the FASTA File ---
def filter_fasta(input_fasta, output_fasta, accessions):
    """Filter FASTA sequences based on a list of Protein Accession IDs."""
    matched_records = []
    
    with open(output_fasta, "w") as out_fasta:
        for record in SeqIO.parse(input_fasta, "fasta"):
            if record.id in accessions:
                matched_records.append(record)
                SeqIO.write(record, out_fasta, "fasta")

    print(f"✅ Saved matched sequences to {output_fasta}")
    return matched_records

def split_fasta(matched_records, m_fasta, without_m_fasta):
    """Split FASTA sequences into two files: one starting with 'M', one without 'M'."""
    m_records = []
    without_m_records = []

    for record in matched_records:
        if record.seq.startswith("M"):
            m_records.append(record)
        else:
            without_m_records.append(record)

    SeqIO.write(m_records, m_fasta, "fasta")
    SeqIO.write(without_m_records, without_m_fasta, "fasta")

    print(f"✅ Saved {len(m_records)} sequences to {m_fasta}")
    print(f"✅ Saved {len(without_m_records)} sequences to {without_m_fasta}")

def extract_methionine_accessions(m_fasta, protein_accession_file, output_excel):
    """Extract gene IDs from M.fasta, match them to Protein accession column, and save matched rows."""
    
    # Step 1: Extract Gene IDs from M.fasta
    m_gene_ids = set(record.id for record in SeqIO.parse(m_fasta, "fasta"))

    # Step 2: Load the protein_accession.xlsx file
    df = pd.read_excel(protein_accession_file)

    # Step 3: Filter rows where 'Protein accession' matches IDs from M.fasta
    matched_df = df[df["Protein accession"].astype(str).isin(m_gene_ids)]

    # Step 4: Save matched rows to an Excel file
    matched_df.to_excel(output_excel, index=False)

    print(f"✅ Saved methionine-matching protein accessions to {output_excel}")


def generate_filtered_summary(filtered_fasta, output_summary):
    """Generate a summary for filtered.fasta, including actual sequences."""
    
    def get_sequence_stats(fasta_file):
        """Helper function to get sequence stats from a FASTA file."""
        records = list(SeqIO.parse(fasta_file, "fasta"))
        lengths = [len(record.seq) for record in records]
        
        if not lengths:
            return 0, "N/A", "N/A", "N/A", 0, 0, 0.0

        longest_seq_record = max(records, key=lambda rec: len(rec.seq))
        shortest_seq_record = min(records, key=lambda rec: len(rec.seq))

        return (
            len(lengths),  # Number of sequences
            longest_seq_record.id,  # ID of longest sequence
            str(longest_seq_record.seq),  # Longest sequence itself
            shortest_seq_record.id,  # ID of shortest sequence
            str(shortest_seq_record.seq),  # Shortest sequence itself
            len(longest_seq_record.seq),  # Length of longest sequence
            len(shortest_seq_record.seq),  # Length of shortest sequence
            sum(lengths) / len(lengths)  # Average length
        )

    stats = get_sequence_stats(filtered_fasta)

    summary_text = (
        "Filtered FASTA Summary\n"
        "=======================\n\n"
        f"Total sequences: {stats[0]}\n"
        f"Longest sequence ID: {stats[1]}\n"
        f"Longest sequence: {stats[2]}\n"
        f"Shortest sequence ID: {stats[3]}\n"
        f"Shortest sequence: {stats[4]}\n"
        f"Length of longest sequence: {stats[5]} bp\n"
        f"Length of shortest sequence: {stats[6]} bp\n"
        f"Average sequence length: {stats[7]:.2f} bp\n"
    )

    with open(output_summary, "w") as summary_file:
        summary_file.write(summary_text)

    print(f"✅ Saved filtered.fasta summary to {output_summary}")


# --- Step 4: Generate Summary File ---
# --- Step 4: Generate Summary File ---
def generate_summary(m_fasta, without_m_fasta, output_summary):
    """Generate a summary for both M.fasta and withoutM.fasta files, including actual sequences."""
    
    def get_sequence_stats(fasta_file):
        """Helper function to get sequence stats from a FASTA file."""
        records = list(SeqIO.parse(fasta_file, "fasta"))
        lengths = [len(record.seq) for record in records]
        
        if not lengths:
            return 0, "N/A", "N/A", "N/A", 0, 0, 0.0

        longest_seq_record = max(records, key=lambda rec: len(rec.seq))
        shortest_seq_record = min(records, key=lambda rec: len(rec.seq))

        return (
            len(lengths),  # Number of sequences
            longest_seq_record.id,  # ID of longest sequence
            str(longest_seq_record.seq),  # Longest sequence itself
            shortest_seq_record.id,  # ID of shortest sequence
            str(shortest_seq_record.seq),  # Shortest sequence itself
            len(longest_seq_record.seq),  # Length of longest sequence
            len(shortest_seq_record.seq),  # Length of shortest sequence
            sum(lengths) / len(lengths)  # Average length
        )

    # Get stats for sequences starting with 'M'
    m_stats = get_sequence_stats(m_fasta)

    # Get stats for sequences NOT starting with 'M', only if the file isn't empty
    if os.path.getsize(without_m_fasta) > 0:
        without_m_stats = get_sequence_stats(without_m_fasta)
    else:
        # If no sequences exist, return default values for the "without M" section
        without_m_stats = (0, "N/A", "N/A", "N/A", "N/A", 0, 0, 0.0)

    summary_text = (
        "Summary Report\n"
        "==============\n\n"
        "Part 1: Sequences starting with 'M'\n"
        f"Number of sequences: {m_stats[0]}\n"
        f"Longest sequence ID: {m_stats[1]}\n"
        f"Longest sequence: {m_stats[2]}\n"
        f"Shortest sequence ID: {m_stats[3]}\n"
        f"Shortest sequence: {m_stats[4]}\n"
        f"Length of longest sequence: {m_stats[5]} bp\n"
        f"Length of shortest sequence: {m_stats[6]} bp\n"
        f"Average sequence length: {m_stats[7]:.2f} bp\n\n"
        "Part 2: Sequences NOT starting with 'M'\n"
        f"Number of sequences: {without_m_stats[0]}\n"
        f"Longest sequence ID: {without_m_stats[1]}\n"
        f"Longest sequence: {without_m_stats[2]}\n"
        f"Shortest sequence ID: {without_m_stats[3]}\n"
        f"Shortest sequence: {without_m_stats[4]}\n"
        f"Length of longest sequence: {without_m_stats[5]} bp\n"
        f"Length of shortest sequence: {without_m_stats[6]} bp\n"
        f"Average sequence length: {without_m_stats[7]:.2f} bp\n"
    )

    with open(output_summary, "w") as summary_file:
        summary_file.write(summary_text)

    print(f"✅ Saved summary with sequences to {output_summary}")

# --- Find Input Files Automatically ---
def find_files():
    """Automatically find FASTA, Excel, and TSV files in the folder."""
    folder = "."
    
    # Find TSV file (ignore generated files)
    tsv_files = [f for f in glob.glob(os.path.join(folder, "*.tsv")) if "with_column_name" not in f]
    tsv_file = tsv_files[0] if tsv_files else None

    # Find Excel file (ignore generated files)
    excel_files = [f for f in glob.glob(os.path.join(folder, "*.xlsx")) if "with_column_name" not in f]
    excel_file = excel_files[0] if excel_files else None

    # Find FASTA file
    fasta_files = glob.glob(os.path.join(folder, "*.fasta"))+ glob.glob(os.path.join(folder, "*.faa"))
    fasta_file = fasta_files[0] if fasta_files else None

    return tsv_file, excel_file, fasta_file

# --- Main Execution ---
def main():
    if len(sys.argv) != 2:
        print("Usage: python pipeline.py '<Species Name>'")
        sys.exit(1)

    species_name = sys.argv[1]
    tsv_file, excel_file, fasta_file = find_files()

    if not tsv_file or not excel_file or not fasta_file:
        print("❌ Error: Missing required input files.")
        sys.exit(1)

    output_tsv_excel = f"{species_name}_with_column_name.xlsx"
    output_protein_accession = f"{species_name}_protein_accession.xlsx"
    matched_fasta = f"{species_name}_filtered.fasta"
    m_fasta = f"{species_name}_M.fasta"
    without_m_fasta = f"{species_name}_withoutM.fasta"
    summary_file = f"{species_name}_summary.txt"

    tsv_df = load_tsv(tsv_file, output_tsv_excel)
    excel_df = pd.read_excel(excel_file)
    protein_accessions = filter_protein_accessions(tsv_df, excel_df, output_protein_accession)
    matched_records = filter_fasta(fasta_file, matched_fasta, protein_accessions)
    split_fasta(matched_records, m_fasta, without_m_fasta)
    generate_summary(m_fasta, without_m_fasta, summary_file)
    filtered_summary_file = f"{species_name}_filtered_summary.txt"
    generate_filtered_summary(matched_fasta, filtered_summary_file)
    # Generate Protein Accession Methionine Excel File
    methionine_output_excel = f"{species_name}_Protein_Accession_Methionine.xlsx"
    extract_methionine_accessions(m_fasta, output_tsv_excel, methionine_output_excel)



if __name__ == "__main__":
    main()
