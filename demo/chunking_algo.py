import time

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
    #check front of chunk
    for x in range(0,patternLoc):
        if (pattern[x] != sequence[sequenceLoc + (x - patternLoc)]) : return False # + 10000000
    #check back of chunk
    for x in range(patternLoc + chunk,len(pattern)):
        if (pattern[x] != sequence[sequenceLoc + (x - patternLoc)]) : return False # + 10000000
    return True


def chunkingAlgo(sequence, pattern):
    start_time = time.time()
    # chunk const has to be less then half of the length of sequence to ensure pattern will definitely hit in the chucking of sequence.
    chunk = int(len(pattern) / 2)
    sequenceSet = buildMap(sequence, chunk)
    patternSet = buildCombination(pattern, chunk)

    #print(patternSet)
    #print(sequenceSet)

    print("--- Pre-Processing %s seconds ---" % (time.time() - start_time))

    start_time = time.time()

    indexlist = list()
    for keyPattern, valuePattern in patternSet.items():
        if keyPattern in sequenceSet:
            for sequenceElement in sequenceSet[keyPattern]:
                for patternElement in valuePattern:
                    if checkPattern(sequence, pattern, sequenceElement, patternElement, chunk):
                        indexlist.append(sequenceElement - patternElement)

    print("--- Processing %s seconds ---" % (time.time() - start_time))

    return set(indexlist)