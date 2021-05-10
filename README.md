# 16S-Taxa-Phlyo Pipeline
Use 16S Sanger sequencing to identify the taxa of culturomics
This pipeline utilize the feature of QIIME 2

![](img/16s-taxa-phylo.png)

## Usage:

### Request

For basic usage of the sever, please refer to [LD Lab BioInfo Wiki](https://github.com/LD-Lab/LD-Lab-BioInfo-Wiki).

You can simply setup the conda env by installing these packages:

```conda env create -n AssignTaxa -f ./environment.yml```

```conda install pandas biopython fastp blast```.
Name                    Version
pandas                    1.2.0
biopython                 1.78
fastp                     0.20.1
blast                     2.10.1

### Data 

For each strain, 4 sequencing result files should be named as follows. ```.seq``` files are not necessary.

```StrainName__27F__*.ab1```

```StrainName__27F__*.seq```

```StrainName__1492R__*.ab1```

```StrainName__1492R__*.seq```

All the files should have the same and uniqe ```StrainName``` as the prefixion. The forward and reverse sequence should have the ```__***F__``` and ```__***R__``` in the filename respectively. You can also specify the ```StrainName``` with a ```metadata``` file.

All sequence results can be put in one ```path```. For example,

```bash
~/16S-Taxa-Phlyo/test_data
```

### Database

The database ```SILVA 16S NR99 Ref``` can be download from [SILVA](https://www.arb-silva.de/) as [SILVA_138.1_SSURef_NR99_tax_silva](https://www.arb-silva.de/fileadmin/silva_databases/release_138_1/Exports/SILVA_138.1_SSURef_NR99_tax_silva.fasta.gz).

Run ```bash init.sh``` to download and initialize the database to ```./Database```.

### Pipeline

The one line command would give the ```taxonomy.tsv``` based on the sequencing result. If you already have the 16S sequence as a ```.fasta``` file, you can skip to QIIME2 taxa classifier section below.

#### One line command

The following command will save the taxa assignment result in ```./test_out``` based on the sequences in ```./test_data```.

```bash
conda activate AssignTaxa
python ./scripts/16S-Taxa-Phylo-Pipeline.py -i ./test_data -o ./test_out
```



