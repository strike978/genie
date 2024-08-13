import streamlit as st
import pandas as pd
import io
import csv


def read_snp_list(snp_file_path):
    snp_order = []
    snp_set = set()
    with open(snp_file_path, mode='r', encoding='utf-8') as snp_file:
        reader = csv.reader(snp_file)
        next(reader)  # Skip the header row
        for row in reader:
            snp = row[1].strip()
            snp_order.append(snp)
            snp_set.add(snp)
    return snp_order, snp_set


def filter_tsv_data(tsv_file, snps):
    snp_lines = {snp: None for snp in snps}
    tsv_data = io.StringIO(tsv_file.decode('utf-8'))
    for line in tsv_data:
        if line.startswith('#'):  # Skip header lines
            continue
        parts = line.strip().split('\t')
        if len(parts) == 4:
            rsid = parts[0].strip()
            if rsid in snp_lines:
                snp_lines[rsid] = line.strip()
    return snp_lines


def main():
    st.title("SNP Data Filter")

    # Path to the SNP list CSV file
    snp_file_path = 'genes.csv'

    # Upload TSV file
    tsv_file = st.file_uploader("Upload TSV File", type=["txt"])

    if tsv_file is not None:
        # Process files
        snp_order, snp_set = read_snp_list(snp_file_path)
        snp_lines = filter_tsv_data(tsv_file.read(), snp_set)

        # Prepare output
        output_lines = ["# rsid\tchromosome\tposition\tgenotype"]
        for snp in snp_order:
            if snp_lines[snp] is not None:
                output_lines.append(snp_lines[snp])

        # Display output in a text area
        st.text_area("Filtered SNP Data", value="\n".join(
            output_lines), height=300)


if __name__ == "__main__":
    main()
