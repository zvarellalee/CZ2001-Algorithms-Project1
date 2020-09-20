import time
from collections import defaultdict


# use a defaultdict to avoid KeyError

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
