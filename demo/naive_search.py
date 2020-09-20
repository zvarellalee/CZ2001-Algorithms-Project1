import time

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