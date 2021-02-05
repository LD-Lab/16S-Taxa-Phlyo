'''
Copyright {2020} Junyu Chen

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

    http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.
'''


import os
import argparse
import collections
import pandas as pd
import subprocess
from Bio import SeqIO

from itertools import repeat
from multiprocessing import Pool, freeze_support

## Generate manifest Table
def manifestGen(InDir, OutDir):
    sampleID = []
    single = {}
    pair = {}
    df = pd.DataFrame()
    SampleList = collections.defaultdict(list)
    for file in os.listdir(InDir):
        filePath = os.path.join(InDir, file)
        if file.endswith(".ab1") and os.path.getsize(filePath) > 0:
            info = file.split("__")
            sampleID.append(info[0])
            SampleList[info[0]].append(filePath)
    for sample in set(sampleID):
        if len(SampleList[sample]) == 1:
            single[sample] = SampleList[sample]
        elif len(SampleList[sample]) == 2:
            pair[sample] = SampleList[sample]
            for end in pair[sample]:
                if "F__" in end:
                    R1 = end
                elif "R__" in end:
                    R2 = end
            df = df.append({"ID": sample, "R1":R1, "R2":R2},ignore_index=True)
    return single, df # return single end dict and pair end dataframe

## Parse pair end abi to fastq file
def parsePairs(df, OutDir):
    R1OutList = []
    R2OutList = []
    df1 = pd.DataFrame()
    fastqDir = os.path.join(OutDir, "fastq")
    if os.path.exists(fastqDir) == 0:
        os.makedirs(fastqDir, 0o777, True)
    for ID in df["ID"]:
        R1OutList.append(os.path.join(fastqDir, ID + "_R1.fastq"))
        R2OutList.append(os.path.join(fastqDir, ID + "_R2.fastq"))
    df1["ID"] = df["ID"]
    df1["R1"] = R1OutList
    df1["R2"] = R2OutList
    return df1

## Parse single end abi to fastq file
def parseSingle(singleDict, OutDir):
    OutList = []
    singleDir = os.path.join(OutDir, "fatsq")
    if os.path.exists(singleDir) == 0:
        os.makedirs(singleDir, 0o777, True)
    for ID in single.keys():
        if "F__" in single[ID][0]:
            OutList.append(os.path.join(singleDir, ID + "_R1.fastq"))
        elif "R__" in single[ID][0]:
            OutList.append(os.path.join(singleDir, ID + "_R2.fastq"))
    return OutList

## abi to fastq
def ABI2Fastq(abi, fastq):
    SeqIO.convert(abi, "abi", fastq, "fastq")
def ABI2FastqParallel(abiList, fastqList):
    pool = Pool(processes = 4)
    pool.starmap(ABI2Fastq, zip(abiList, fastqList))
    pool.close()
    pool.join()
    pool.terminate()

## Fastp pair end tirmming and merging
def RunFastp(R1, R2, prefix, OutDir):
    mergeDir = os.path.join(OutDir, "merge")
    fastpDir = os.path.join(OutDir, "fastp")
    unmergeDir = os.path.join(OutDir, "unmerge", "fastq")
    if os.path.exists(OutDir) == 0:
        os.makedirs(OutDir, 0o777, True)
    if os.path.exists(mergeDir) == 0:
        os.makedirs(mergeDir, 0o777, True)
    if os.path.exists(fastpDir) == 0:
        os.makedirs(fastpDir, 0o777, True)
    if os.path.exists(unmergeDir) == 0:
        os.makedirs(unmergeDir, 0o777, True)
    cmd = "fastp -i " + R1 + " -I " + R2 + " -o " + os.path.join(unmergeDir, prefix + "_R1.fastq") + " -O " + os.path.join(unmergeDir, prefix + "_R2.fastq") + \
    " --trim_front1 30 --max_len1 750 --trim_front2 30 --max_len2 750 --cut_front --cut_tail --cut_window_size 20 --cut_mean_quality 30" + \
    " --merge --merged_out " + os.path.join(mergeDir, prefix + ".fastq") + " --correction --overlap_len_require 20 " + \
    " --html " + os.path.join(fastpDir, prefix + ".html") + " --json " + os.path.join(fastpDir, prefix + ".json") + " --report_title " + prefix + "-fastq-merge-report"
    subprocess.call(cmd, shell=True)
def RunFastpParallel(R1List, R2List, prefixList, OutDir):
    pool = Pool(processes = 4)
    pool.starmap(RunFastp, zip(R1List, R2List, prefixList, repeat(OutDir)))
    pool.close()
    pool.join()
    pool.terminate()

## fastp single end trimming and merging
def RunFastpSingle(single, OutDir):
    prefix = os.path.split(single)[1].replace(".fastq", "")
    trimDir = os.path.join(OutDir, "trim")
    fastpDir = os.path.join(OutDir, "fastp")
    if os.path.exists(OutDir) == 0:
        os.makedirs(OutDir, 0o777, True)
    if os.path.exists(trimDir) == 0:
        os.makedirs(trimDir, 0o777, True)
    if os.path.exists(fastpDir) == 0:
        os.makedirs(fastpDir, 0o777, True)        
    cmd = "fastp -i " + single + " -o " + os.path.join(trimDir, prefix + ".fastq")  + \
    " --trim_front1 30 --cut_front --cut_tail --cut_window_size 20 --cut_mean_quality 30" + \
    " --html " + os.path.join(fastpDir, prefix + ".html") + " --json " + os.path.join(fastpDir, prefix + ".json") + " --report_title " + prefix + "-fastq-merge-report"
    subprocess.call(cmd, shell=True)
def RunFastpSingleParallel(singleList, OutDir):
    pool = Pool(processes = 4)
    pool.starmap(RunFastpSingle, zip(singleList, repeat(OutDir)))
    pool.close()
    pool.join()
    pool.terminate()

## parse pairs fastq to fasta file
def parseFastqPairs(InDir, OutDir):
    allSeq = []
    fastaFileList = []
    fastaDir = os.path.join(OutDir, "fasta")
    if os.path.exists(fastaDir) == 0:
        os.makedirs(fastaDir, 0o777, True)
    for file in os.listdir(InDir):
        if file.endswith(".fastq") and os.path.getsize(os.path.join(InDir, file)) > 0:
            for seq in SeqIO.parse(os.path.join(InDir, file), "fastq"):
                seq.id = file.replace(".fastq", "")
                seq.name = ""
                seq.description = seq.description.split(" ")[1]
                fastaFile = os.path.join(fastaDir, file.replace(".fastq", ".fasta"))
                fastaFileList.append(fastaFile)
                SeqIO.write(seq, fastaFile, "fasta")
                allSeq.append(seq)
    SeqIO.write(allSeq, os.path.join(OutDir, "pairsMerge.fasta"), "fasta")
    return fastaFileList

## Parse fastq file to fasta file
def parseFastqSingle(InDir, OutDir):
    allSeq = []
    fastaFileList = []
    fastaDir = os.path.join(OutDir, "fasta")
    if os.path.exists(fastaDir) == 0:
        os.makedirs(fastaDir, 0o777, True)
    for file in os.listdir(InDir):
        if file.endswith(".fastq") and os.path.getsize(os.path.join(InDir, file)) > 0:
            for seq in SeqIO.parse(os.path.join(InDir, file), "fastq"):
                seq.id = file.replace(".fastq", "")
                seq.name = ""
                seq.description = ""
                if "_R2" in seq.id:
                    seq.seq = seq.seq.reverse_complement()
                allSeq.append(seq)
                fastaFile = os.path.join(fastaDir, file.replace(".fastq", ".fasta"))
                fastaFileList.append(fastaFile)
                SeqIO.write(seq, fastaFile, "fasta")
    SeqIO.write(allSeq, os.path.join(OutDir, "single.fasta"), "fasta")
    return fastaFileList

## Parse unmerge fastq
def parseFastqUnmerge(filePath, OutDir):
    allSeq = []
    fastaFileList = []
    fastaDir = os.path.join(OutDir, "fasta")
    if os.path.exists(fastaDir) == 0:
        os.makedirs(fastaDir, 0o777, True)
    for file in filePath:
        for seq in SeqIO.parse(file, "fastq"):
            seq.id = os.path.split(file)[1].replace(".fastq", "")
            seq.name = ""
            seq.description = ""
            if "_R2" in seq.id:
                seq.seq = seq.seq.reverse_complement()
            allSeq.append(seq)
            fastaFile = os.path.join(fastaDir, os.path.split(file)[1].replace(".fastq", ".fasta"))
            fastaFileList.append(fastaFile)
            SeqIO.write(seq, fastaFile, "fasta")
    SeqIO.write(allSeq, os.path.join(OutDir, "unmerge.fasta"), "fasta")
    return fastaFileList


def RunBlastnParallel(fastaList, db, jobs, threads, OutDir):
    pool = Pool(processes=jobs)
    pool.starmap(RunBlastn, zip(fastaList, repeat(db), repeat(threads), repeat(OutDir)))
    pool.close()
    pool.join()
    pool.terminate()
def RunBlastn(fasta, db, threads, OutDir):
    blastDir = os.path.join(OutDir, "blast")
    if os.path.exists(blastDir) == 0:
        os.makedirs(blastDir, 0o777, True)
    OutFile = os.path.join(blastDir, os.path.split(fasta)[1].replace(".fasta", "") + "_blast.tsv")
    cmd = "blastn -query " + fasta  + " -out " + OutFile + " -evalue 1.0 -max_target_seqs 5 -outfmt 6 -db " + db + " -num_threads " + str(threads) 
    print(cmd)
    subprocess.call(cmd, shell=True)

def parseIndex(blastDir, OutDir):
    taxaList = []
    df = pd.DataFrame()
    out = pd.DataFrame()
    for file in os.listdir(blastDir):
        #print(file)
        if file.endswith(".tsv") and os.path.getsize(os.path.join(blastDir, file)) > 0:
            df1 = pd.read_table(os.path.join(blastDir, file), header = None)
            df1.columns = ["qseqid", "sseqid", "pident", "length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore"]
            df = df.append(df1[0:1])
    df = df.reset_index()
    df2 = pd.read_table(indexPath)
    index = dict(zip(df2["ID"], df2["Taxa"]))
    for i in range(len(df)):
        taxa = index[df["sseqid"][i]]
        taxaList.append(taxa)
    out["qseqid"] = df["qseqid"]
    out["taxa"] = taxaList
    out["pident"] = df["pident"]
    out["length"] = df["length"]
    return out

parser = argparse.ArgumentParser(description='Run Diamond')
parser.add_argument('-i', '--input', dest='InDir', type=str, required=True,
                    help="the path of the reads")
parser.add_argument('-o', '--output', dest='OutDir', type=str, required=True,
                    help="the output path of reads")
parser.add_argument('-d', '--database', dest='database', type=str,  required=False, default="/mnt/d/Lab/16S-Taxa-Phlyo/database/SILVA_138.1_SSURef_NR99_tax_silva.fasta",
                    help="the reference_reads path")
parser.add_argument('-r', '--index', dest='index', type=str,  required=False, default="/mnt/d/Lab/16S-Taxa-Phlyo/database/silva-138-99-index.tsv",
                    help="the reference_taxonomy path")
parser.add_argument('-j', '--jobs', dest='jobs', type=str,  required=False, default='2',
                    help="the number of jobs run in parallel")
parser.add_argument('-t', '--threads', dest='threads', type=str, required=False, default='4',
                    help="the number of threads run for a job")
args = parser.parse_args()

InDir = os.path.abspath(args.InDir)
OutDir = os.path.abspath(args.OutDir)
db = os.path.abspath(args.database)
indexPath = os.path.abspath(args.index)
jobs = int(args.jobs)
threads = int(args.threads)

if os.path.exists(OutDir) == 0:
    os.makedirs(OutDir, 0o777, True)
singleOutDir = os.path.join(OutDir, "single")
if os.path.exists(singleOutDir) == 0:
    os.makedirs(singleOutDir, 0o777, True)


single, pairs = manifestGen(InDir, OutDir)

pairs.to_csv(os.path.join(OutDir, "pairsTable.csv"))
pairs_q = parsePairs(pairs, OutDir)
## abi to fastq
ABI2FastqParallel(pairs["R1"], pairs_q["R1"])
ABI2FastqParallel(pairs["R2"], pairs_q["R2"])
## fastp
RunFastpParallel(pairs_q["R1"], pairs_q["R2"], pairs_q["ID"], OutDir)
## fastq to fasta
fastaFileList = parseFastqPairs(os.path.join(OutDir, "merge"), OutDir)
## blastn
RunBlastnParallel(fastaFileList, db, jobs, threads, OutDir)
## blast tsv parse
pairsOut = parseIndex(os.path.join(OutDir, "blast"), OutDir)
pairsOut.to_csv(os.path.join(OutDir, "pairsOut.csv"))

# Single
if len(single) > 0:
    singleList = [item for sublist in list(single.values()) for item in sublist]
    single_q_List = parseSingle(single, singleOutDir)
    ABI2FastqParallel(singleList, single_q_List)
    RunFastpSingleParallel(single_q_List, singleOutDir)
    fastaFileList = parseFastqSingle(os.path.join(singleOutDir, "trim"), singleOutDir)
    RunBlastnParallel(fastaFileList, db, jobs, threads, singleOutDir)
    singleOut = parseIndex(os.path.join(singleOutDir, "blast"), singleOutDir)
    singleOut.to_csv(os.path.join(OutDir, "singleOut.csv"))
    pairsOut = pairsOut.append(singleOut)

# Unmerge
unmergeOutDir = os.path.join(OutDir, "unmerge")
unmergeFastqDir = os.path.join(unmergeOutDir, "fastq")
if os.path.exists(unmergeOutDir) == 0:
    os.makedirs(unmergeOutDir, 0o777, True)

unmergeList = []
for file in os.listdir(unmergeFastqDir):
    filePath = os.path.join(unmergeFastqDir, file)
    if file.endswith(".fastq") and os.path.getsize(filePath) > 0:
        unmergeList.append(filePath)
if len(unmergeList) > 0:
    fastaFileList = parseFastqUnmerge(unmergeList, unmergeOutDir)
    RunBlastnParallel(fastaFileList, db, jobs, threads, unmergeOutDir)
    unmergeOut = parseIndex(os.path.join(unmergeOutDir, "blast"), unmergeOutDir)
    unmergeOut.to_csv(os.path.join(OutDir, "unmergeOut.csv"))
    pairsOut = pairsOut.append(unmergeOut)

pairsOut.to_csv(os.path.join(OutDir, "finalOut.csv"))