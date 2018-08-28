def Count(Motifs):
    # create count_matrix that represents the presence of each nuclotide
    # A,C,G,and T in each column of group of motifs
    count = {}
    k = len(Motifs[0])
    for symbol in "ACGT":
        count[symbol] = []
        for j in range(k):
            count[symbol].append(0)

    t = len(Motifs)
    for i in range(t):
        for j in range(k):
            symbol = Motifs[i][j]
            count[symbol][j] += 1

    return count

def Profile(Motifs):
    profile_matrix = Count(Motifs)
    l = len(Motifs)
    for symbol in "ACGT":
        for j in range(len(profile_matrix[symbol])):
            profile_matrix[symbol][j] = (profile_matrix[symbol][j])/l

    return profile_matrix



motifs = ["AACGTA", "CCCGTT", "CACCTT", "GGATTA", "TTCCGG"]
print(Profile(motifs))
