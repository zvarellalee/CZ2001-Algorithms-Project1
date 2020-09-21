from Bio import SeqIO
from collections import defaultdict
import pandas as pd
import time
import os
# -------- Naive Algorithm --------
def naivePatternSearch(sequence, pattern):
    #Simple linear Pattern Searching
    start_time = time.time()
    indexList = list() # list for storing matched pattern starting index
    patternLen = len(pattern)
    for keyParent, valueParent in sequence.items(): # iter over every char in sequence (split key/value)
        for keyChild, valueChild in pattern.items(): # iter over every char in searched pattern (split key/value)
            try:
                # check if following individual chars matches all chars in searched pattern break internal for loop if a single char does not match.
                # if current index char and following char fulfill searched pattern append to indexlist
                if sequence[keyParent + keyChild] != valueChild: break
                elif (patternLen == (keyChild + 1)): indexList.append(keyParent)
            except: break
    print("--- Processing %s seconds ---" % (time.time() - start_time))
    return set(indexList)

# -------- Z Algorithm --------
def calZArr(myStr, patLen, Z, lookupDict):
    strLen = len(myStr)
    L, R = 0, 0
    count = 0
    for i in range(1, strLen):
        if i > R:
            # reset left and R, prefix should start from i
            L, R = i, i
            # compare S[0...] and S[i...]
            # S[R-L] == S[0], increment R will traverse array until index R
            while (R < strLen) and (myStr[R] == myStr[R - L]):
                R += 1
            Z[i] = R - L
            # deduct 1 as incremented R to compare next (S[i], S[R-L]) but not equal
            R -= 1
        else:
            k = i - L
            # S[i] == S[k]
            # Z[i] = min(Z[k], R - i +1)
            # compute {L, R}
            if Z[k] < R - i + 1:
                # [L, R] stays same
                Z[i] = Z[k]
            else:
                # reset
                L = i
                while (R < strLen) and (myStr[R - L] == myStr[R]):
                    R += 1
                Z[i] = R - L
                R -= 1
        # save in dictionary, Key = Z-Value at the index, Value = index in sequence
        lookupDict[Z[i]].append(i - patLen - 1)
    return lookupDict
    # return Z

def ZAlgo(text, pattern):
    # create string for finding Z array
    concatStr = pattern + "$" + text
    myLen = len(concatStr)
    patLen = len(pattern)

    # construct Z array
    Z = [0] * myLen
    lookupDict = defaultdict(list)
    start_time = time.time()
    zDict = calZArr(concatStr, patLen, Z, lookupDict)
    print("--- Pre-Processing %s seconds ---" % (time.time() - start_time))

    start_time = time.time()
    # retrieve from dictionary any potential matches
    matches = zDict[patLen]
    print("--- Processing %s seconds ---" % (time.time() - start_time))
    return matches

# -------- Chunking Algorithm --------
def buildMap(sequence, chunk):
    #expect m/c time complexity; m - sequence length, c - chunk size
    m = len(sequence)
    ptr = 0
    Map = {}
    while ptr < m :
        string = "".join(sequence[ptr:(ptr + chunk)].to_numpy())
        if string in Map:
            Map[string].append(ptr)
        else:
            Map[string]=[ptr]
        ptr += chunk
    return Map

def buildCombination(pattern, chunk):
    #expect n-(c-1) time complexity; n - pattern length, c - chunk size
    m = len(pattern)
    ptr = 0
    Map = {}
    while ptr + chunk <= m :
        string = "".join(pattern[ptr:(ptr + chunk)].to_numpy())
        if string in Map:
            Map[string].append(ptr)
        else:
            Map[string]=[ptr]
        ptr += 1
    return Map

def checkPattern(sequence, pattern, sequenceLoc, patternLoc, chunk):
    # expect n - c time complexity
    sequenceLength = len(sequence)
    #check front of chunk
    for x in range(0,patternLoc):
        sequenceIndex = sequenceLoc + (x - patternLoc) + sequence.index[0]
        if (0 > sequenceIndex or sequenceIndex >= sequenceLength) : return False
        elif (pattern[x] != sequence[sequenceIndex]) : return False
    #check back of chunk
    for x in range(patternLoc + chunk,len(pattern)):
        sequenceIndex = sequenceLoc + (x - patternLoc) + sequence.index[0]
        if (0 > sequenceIndex or sequenceIndex >= sequenceLength) : return False
        elif (pattern[x] != sequence[sequenceIndex]) : return False
    return True

def chunkingAlgo(sequence, pattern):
    start_time = time.time()

    # chunk const has to be less then half of the length of sequence to ensure pattern will definitely hit in the chucking of sequence.
    chunk = int(len(pattern) / 2)
    sequenceSet = buildMap(sequence, chunk)

    print("--- Pre-Processing %s seconds ---" % (time.time() - start_time))

    start_time = time.time()

    patternSet = buildCombination(pattern, chunk)

    start_time = time.time()

    indexlist = list()
    for keyPattern, valuePattern in patternSet.items():
        if keyPattern in sequenceSet:
            for sequenceElement in sequenceSet[keyPattern]:
                for patternElement in valuePattern:
                    if checkPattern(sequence, pattern, sequenceElement, patternElement, chunk):
                        indexlist.append(sequenceElement - patternElement + sequence.index[0])

    print("--- Processing %s seconds ---" % (time.time() - start_time))

    return set(indexlist)

# -------- Other Functions --------
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

# -------- Main --------
searchedPattern = "CCTACAA"
subset_size = 100000  # 10000000    # size of genome sequence subset
seqSeries = pd.Series(dtype=object)
stringSequence = ''
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
        if len(availableFasta) != 0:
            print("Please enter file name: ")
            fileName = input()
            print("Enter subset size(0 for complete genome, current is %d): "%subset_size)
            subset_size = int(input())
            try:
                stringSequence, seqSeries = fileselect(fileName,subset_size)
            except(FileNotFoundError):
                print("File not found in directory. Please make sure the correct name is inputted.")
        else:
            print("No FASTA files found inside Sequences directory.")
    elif c==2:
        print("Enter pattern to search(Current is %s): "%searchedPattern) #CCTACAA
        searchedPattern = input()
    elif c==3:
        print('Naive Search Algorithm')
        print(sorted(naivePatternSearch(seqSeries, pd.Series([x for x in searchedPattern]))))
    elif c==4:
        print('Z Algorithm')
        print(ZAlgo(stringSequence, searchedPattern))
    elif c==5:
        print('Chunking Algorithm')
        print(sorted(chunkingAlgo(seqSeries, pd.Series(list(searchedPattern)))))
    elif c==6:
        print('Naive Search Algorithm')
        print(sorted(naivePatternSearch(seqSeries, pd.Series([x for x in searchedPattern]))))
        print('\nZ Algorithm')
        print(ZAlgo(stringSequence, searchedPattern))
        print('\nChunking Algorithm')
        print(sorted(chunkingAlgo(seqSeries, pd.Series([x for x in searchedPattern]))))
    elif c==7:
        active = False
    print()