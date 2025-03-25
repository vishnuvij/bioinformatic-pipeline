1.Overview

This Python script (pipeline.py) is a bioinformatics pipeline that processes:

A FASTA/FAA file (contains biological sequences).
A TSV file (contains sequence metadata which is intreproscan output tsv file).
An Excel file (contains PFAM codes for filtering).

The script automates file detection, data processing, sequence filtering, and summary generation to output multiple useful files.

2.Input Files

Place the script (pipeline.py) inside a folder containing:

File Type	Example Name	Description
FASTA /FAA file	sequences.fasta	Contains protein sequences with unique IDs.
TSV file	interpro_output.tsv	Metadata table (no headers), includes Protein accession, Signature accession, etc.
Excel file	pfam_data.xlsx	Contains PFAM Codes used for filtering.
2.a File Detection
The script automatically finds the FASTA, TSV, and Excel files, regardless of their names.



3.Running the script

Install all the required dependencies.Place all the input files (fasta , tsv and excel) along with the pipeline.py script in a single folder. Open bash and execute the following command.

python pipeline.py "Species Name"   (Replace 'species name' with the name of the species the fasta file is associated to )

4.Expected Output

After running the script it is expected to generate 8 output files as follows. The run time depends on the size of input files. The following messages are expected to be shown in Terminal if the script is able to process the files correctly. 

✅ Saved TSV as Excel: Species Name_with_column_name.xlsx

✅ Saved filtered protein accessions to Species Name_protein_accession.xlsx

✅ Saved matched sequences to Species Name_filtered.fasta

✅ Saved filtered.fasta summary to Species Name_filtered_summary.txt

✅ Saved X sequences to Species Name_M.fasta

✅ Saved Y sequences to Species Name_withoutM.fasta

✅ Saved summary with sequences to Species Name_summary.txt

✅ Saved methionine-matching protein accessions to Species Name_Protein_Accession_Methionine.xlsx



5.Explanation of Output Files in the Bioinformatics Pipeline
Species Name_with_column_name.xlsx : Saving file as excel file after renaming with column headers for better file processing.
Species Name_protein_accession.xlsx : This file contains a filtered list of protein accessions, based on the PFAM codes.

Species Name_filtered.fasta : This FASTA file is derived from the original input fasta file and contains only sequences that match the extracted protein accession IDs from protein_accession.xlsx.
Species Name_filtered_summary.txt : This summary provides intermediate general statistics for filtered.fasta, allowing quick validation of the filtering process.
Species Name_M.fasta : This file is created from Species_Name_filtered.fasta. This new fasta contains only sequences that start with character 'M' in Species_Name_filtered.fasta.
Species Name_withoutM.fasta : This file is created from Species_Name_filtered.fasta. This new fasta contains only sequences that dosent start with character 'M' in Species_Name_filtered.fasta .
Species Name_summary.txt: This file provides detailed statistics for Species_Name_M.fasta and Species_name_withoutM.fasta, including sequence count, longest/shortest sequences, and their IDs.
Species Name_Protein_Accession_Methionine.xlsx : This is basically a subset of  Species_Name_with_column_name.xlsx. It is derived by matching the id's to Protein accession column of Species_Name_with_column_name.xlsx and all the row's with a match is extracted.

