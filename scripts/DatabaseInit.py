import os
import argparse
import pandas as pd
from Bio import SeqIO

parser = argparse.ArgumentParser(description='Run Diamond')

parser.add_argument('-i', '--input', dest='DatabaseFasta', type=str, required=False, default="./Database/SILVA_138.1_SSURef_NR99_tax_silva.fasta", 
                    help="the path of the database")
args = parser.parse_args()
databaseFasta = os.path.abspath(args.DatabaseFasta)
DatabaseDir = os.path.dirname(databaseFasta)

idList = []
taxaList = []
for seq in SeqIO.parse(databaseFasta, "fasta"):
    idList.append(seq.id)
    taxaList.append(seq.description.replace(seq.id + " " , ""))

df = pd.DataFrame()
df["ID"] = idList
df["Taxa"] = taxaList
df.to_csv(os.path.join(DatabaseDir, "silva-138-99-index.tsv"), index=None, sep="\t")
