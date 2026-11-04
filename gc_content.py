import argparse
from Bio import SeqIO
import pandas as pd
from pathlib import Path

import matplotlib.pyplot as plt
import seaborn as sns

parser = argparse.ArgumentParser(description='Calculate gc content from fasta files')
parser.add_argument('-i','--input_dir', type=str, default='.', help='PATH to the directory cotains fasta files', required=True)
parser.add_argument('-o','--out_dir', type=str, default='.', help='PATH to the output directory', required=True)
args = parser.parse_args()
input_dir = args.input_dir
out_dir = args.out_dir
print("input directory:" + input_dir)
print("output directory:" + out_dir)

####################looking for fasta files############################
def find_fasta(input_dir):
    data_dir = Path(input_dir)
    fasta_files = [
        file
        for file in data_dir.glob("*")
        if file.is_file()
        and file.suffix.lower() in {".fasta", ".fa", ".fna"}
    ]

    return fasta_files
################################################

###################calculation#############################
def calculate_gc(fasta, out_dir):
    out_path = Path(out_dir)
    gc_results = []

    with open(fasta, 'r') as file:
        for record in SeqIO.parse(file, 'fasta'):
            id_part = record.id
            desc_part = record.description
            seq = record.seq

            #change the all sequence letters into capital letters (just in case)
            seq = str(seq).upper()
            contig_len = len(seq)
            gc = seq.count("G") + seq.count("C")
            gc_content = gc/contig_len * 100

            gc_results.append({
                'ID': id_part,
                'Length' : contig_len,
                'GC_count' : gc,
                'GC_content' : gc_content
            })
    df = pd.DataFrame(gc_results)
    file_name = f"{fasta.stem}_gc_content.tsv"
    df.to_csv(out_path/file_name, sep="\t", index=False)

    return df

df = pd.read_csv("gc_content.tsv", sep="\t")
################################################


####################Visualization############################
def plot_gc(df, fasta, out_dir):
    out_path = Path(out_dir)
    plt.figure(figsize=(8, 5))

    sns.histplot(
        data=df,
        x="GC_content",
        binwidth=1
    )

    plt.xlabel("GC content (%)")
    plt.ylabel("Number of sequences")
    plt.title(f"GC content distribution: {fasta.name}")

    # Save PNG
    png_file = f"{fasta.stem}_gc_content.png"
    plt.savefig(out_path/png_file, dpi=300, bbox_inches="tight")

    plt.close()
################################################

################################################

fasta_files = find_fasta(input_dir) #This function automatically detect your fasta files

for fasta in fasta_files:

    print(f"Processing: {fasta.name}")

    df = calculate_gc(fasta,out_dir)

    if df.empty:
        continue

    plot_gc(df, fasta, out_dir)
################################################