import streamlit as st
import pandas as pd
import io
import csv


def read_snp_list(snp_file_path):
    gene_snps = {}
    with open(snp_file_path, mode='r', encoding='utf-8') as snp_file:
        reader = csv.reader(snp_file)
        next(reader)  # Skip the header row
        for row in reader:
            gene = row[0].strip()
            snp = row[1].strip()
            ancestral = row[3].strip()
            derived = row[2].strip()
            if gene not in gene_snps:
                gene_snps[gene] = {}
            gene_snps[gene][snp] = {'ancestral': ancestral, 'derived': derived}
    return gene_snps


def filter_tsv_data(tsv_file, gene_snps):
    snp_lines = {gene: {snp: None for snp in snps}
                 for gene, snps in gene_snps.items()}
    tsv_data = io.StringIO(tsv_file.decode('utf-8'))
    for line in tsv_data:
        if line.startswith('#'):  # Skip header lines
            continue
        parts = line.strip().split('\t')
        if len(parts) == 4:
            rsid = parts[0].strip()
            for gene, snps in gene_snps.items():
                if rsid in snps:
                    snp_lines[gene][rsid] = line.strip()
    return snp_lines


def calculate_allele_percentages(snp_lines, gene_snps):
    gene_percentages = {}
    for gene, snps in snp_lines.items():
        total_alleles = 0
        ancestral_count = 0
        derived_count = 0

        # Debugging statement to check SNPs for each gene
        print(f"Gene: {gene}, Total SNPs: {len(snps)}")

        for snp, line in snps.items():
            if line is not None:
                parts = line.split('\t')
                genotype = parts[3].strip()
                ancestral_allele = gene_snps[gene][snp]['ancestral']
                derived_allele = gene_snps[gene][snp]['derived']

                # Debugging statements
                print(f"Gene: {gene}, SNP: {snp}, Genotype: {genotype}, Ancestral: {
                      ancestral_allele}, Derived: {derived_allele}")

                # Count each allele in the genotype
                for allele in genotype:
                    total_alleles += 1
                    if allele == ancestral_allele:
                        ancestral_count += 1
                    elif allele == derived_allele:
                        derived_count += 1

                # Additional debugging to check counts
                print(f"Ancestral Count: {ancestral_count}, Derived Count: {
                      derived_count}, Total Alleles: {total_alleles}")

        ancestral_percentage = (
            ancestral_count / total_alleles) * 100 if total_alleles > 0 else 0
        derived_percentage = (derived_count / total_alleles) * \
            100 if total_alleles > 0 else 0

        gene_percentages[gene] = {
            'ancestral_percentage': ancestral_percentage,
            'derived_percentage': derived_percentage
        }

        # Debugging statement to check final percentages
        print(f"Gene: {gene}, Ancestral Percentage: {
              ancestral_percentage}, Derived Percentage: {derived_percentage}")

    return gene_percentages


def main():
    st.set_page_config(layout="wide")
    st.title("SNP Data Filter")

    # Path to the SNP list CSV file
    snp_file_path = 'genes.csv'

    # Upload TSV file
    tsv_file = st.file_uploader("Upload TSV File", type=["txt"])

    if tsv_file is not None:
        # Process files
        gene_snps = read_snp_list(snp_file_path)
        snp_lines = filter_tsv_data(tsv_file.read(), gene_snps)

        # Prepare output
        output_lines = ["# rsid\tchromosome\tposition\tgenotype"]
        for gene, snps in gene_snps.items():
            for snp in snps:
                if snp_lines[gene][snp] is not None:
                    output_lines.append(snp_lines[gene][snp])

        # Calculate allele percentages
        gene_percentages = calculate_allele_percentages(snp_lines, gene_snps)

        # Display output in a text area
        st.text_area("Filtered SNP Data", value="\n".join(
            output_lines), height=300)

        # Prepare data for datatable
        data = []
        for gene, percentages in gene_percentages.items():
            data.append({
                'Gene': gene,
                'Derived Allele Percentage': f"{percentages['derived_percentage']:.2f}%",
                'Ancestral Allele Percentage': f"{percentages['ancestral_percentage']:.2f}%"
            })

        # Create DataFrame
        df = pd.DataFrame(data)

        # Display DataFrame in Streamlit
        st.dataframe(df, use_container_width=True)


if __name__ == "__main__":
    main()
