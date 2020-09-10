def RabinKarpSearch(sequence, pattern, p):  # constant prime p for limiting maximum hash value
    m = len(pattern) # length of pattern/window size
    n = len(sequence) # length of sequence
    d = 4 # only 4 letters in DNA alphabet - A, C, G, T
    pattern_hash = 0
    seq_h = 0
    index_list = []

    # compute hash value of pattern and first window of sequence
    for i in range(m):
        pattern_hash += ord(pattern[i]) * pow(d, m-i-1)
        seq_h += ord(sequence[i]) * pow(d, m-i-1)
    pattern_hash %= p
    seq_hash = seq_h % p
    # go through whole sequence and append all indexes where pattern equals to sequence window to index list
    for i in range(n-m+1):
        j = i
        if (pattern_hash == seq_hash):
            for k in range(m):
                if (pattern[k] != sequence[j]):
                    break
                if (k == m-1):
                    index_list.append(i)
                j += 1
        if (i == n-m):
            break
        # move sequence window
        seq_h = d * (seq_h - (ord(sequence[i]) * pow(d, m-1))) + (ord(sequence[i+m]) * pow(d, 0))
        seq_hash = seq_h % p

    if (index_list == []):
        print("empty")
    else:
        print(index_list)

RabinKarpSearch("ACGTACGTAGGTACCCGGGGTACCGTTCACCTT", "ACC", 101)
