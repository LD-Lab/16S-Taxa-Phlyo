# Conda 
# conda env create -n AssignTaxa -f ./environment.yml

mkdir Database
# Download the database from SILVA
wget -P ./Database https://www.arb-silva.de/fileadmin/silva_databases/release_138_1/Exports/SILVA_138.1_SSURef_NR99_tax_silva.fasta.gz
gunzip ./Database/SILVA_138.1_SSURef_NR99_tax_silva.fasta.gz

# Prepare index file
python ./scripts/DatabaseInit.py -i ./Database/SILVA_138.1_SSURef_NR99_tax_silva.fasta

# make blast db
makeblastdb -in ./Database/SILVA_138.1_SSURef_NR99_tax_silva.fasta -parse_seqids -dbtype nucl

