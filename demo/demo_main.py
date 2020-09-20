from Bio import SeqIO
import pandas as pd
import time
import chunking_algo
import naive_search
import z_algo
import os

def fileselect(fileName,subset_size):
    start_time = time.time()
    records = []
    for record in SeqIO.parse("Sequences/"+fileName, "fasta"):
        records.append(str(record.seq))

    print("--- Fasta load: %s seconds ---" % (time.time() - start_time))
    stringSequence = "".join(records)
    totalLen = sum(len(s) for s in records)
    print("Total length of sequence: ", totalLen)
    if subset_size != 0:
        stringSequence = stringSequence[0:subset_size]
    seqSeries = pd.Series(list(stringSequence))
    #print(seqSeries[0:300])
    print("Total length of subset: ", len(seqSeries))
    print("Genome sequence processing done.\n")
    print("Sample: " + stringSequence[0:200])
    return stringSequence, seqSeries

searchedPattern = "CCTACAA"
subset_size = 100000  # 10000000    # size of genome sequence subset

active = True
while active:
    print("----- Select option (Enter number)-----\n" +
          "(1) Fasta file selection\n" +
          "(2) Change search pattern\n" +
          "(3) Naive Algorithm\n" +
          "(4) Z Algorithm\n" +
          "(5) Chunking Algorithm\n" +
          "(6) Run All Algorithms\n" +
          "(7) Quit")
    c = int(input())
    if c==1:
        availableFasta = []
        for file in os.listdir("Sequences"):
            if not file.endswith(".txt"):
                availableFasta.append(file)
        print("Genome sequences in directory: ", availableFasta)
        print("Please enter file name: ")
        fileName = input()
        print("Enter subset size(0 for complete genome, current is %d): "%subset_size)
        subset_size = int(input())
        stringSequence, seqSeries = fileselect(fileName,subset_size)
    elif c==2:
        print("Enter pattern to search: ") #CCTACAA
        searchedPattern = input()
    elif c==3:
        print('Naive Search Algorithm')
        print(sorted(naive_search.naivePatternSearch(seqSeries, pd.Series([x for x in searchedPattern]))))
    elif c==4:
        print('Z Algorithm')
        print(z_algo.ZAlgo(stringSequence, searchedPattern))
    elif c==5:
        print('Chunking Algorithm')
        print(sorted(chunking_algo.chunkingAlgo(seqSeries, pd.Series(list(searchedPattern)))))
    elif c==6:
        print('Naive Search Algorithm')
        print(sorted(naive_search.naivePatternSearch(seqSeries, pd.Series(list(searchedPattern)))))
        print('\nZ Algorithm')
        print(z_algo.ZAlgo(stringSequence, searchedPattern))
        print('\nChunking Algorithm')
        print(sorted(chunking_algo.chunkingAlgo(seqSeries, pd.Series([x for x in searchedPattern]))))
    elif c==7:
        active = False
    print()



