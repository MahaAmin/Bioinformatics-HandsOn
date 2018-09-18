# Input:  A set of kmers Motifs
# Output: Count(Motifs)
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

# Input:  A list of kmers Motifs
# Output: the profile matrix of Motifs, as a dictionary of lists.
def Profile(Motifs):
    profile_matrix = Count(Motifs)
    l = len(Motifs)
    for symbol in "ACGT":
        for j in range(len(profile_matrix[symbol])):
            profile_matrix[symbol][j] = (profile_matrix[symbol][j])/l

    return profile_matrix


def Consensus(Motifs):
    k = len(Motifs[0])
    count = Count(Motifs)

    consensus = ""
    for j in range(k):
        m = 0
        frequencySymbol = ""
        for symbol in "ACGT":
            if count[symbol][j] > m :
                m = count[symbol][j]
                frequencySymbol = symbol
        consensus += frequencySymbol
    return consensus

# Input:  A set of k-mers Motifs
# Output: The score of these k-mers.
def Score(Motifs):
    consensus = Consensus(Motifs)
    k = len(Motifs)
    score = 0
    for i in range(len(Motifs[0])):
        for j in range(k):
            if(Motifs[j][i] != consensus[i]):
                score += 1
    return score


def Pr(Text, Profile):
    prob = 1
    for i in range(len(Text)):
        prob = prob * Profile[Text[i]][i]
    return prob


def ProfileMostProbablePattern(Text, k, Profile):
    best_p = 0
    best_kmer = Text[0:k]
    for i in range(len(Text)-k+1):
        p = Pr(Text[i:i+k], Profile)
        if(p > best_p):
            best_p = p
            best_kmer = Text[i:i+k]
    return best_kmer

text = "ACGGGGATTACC"
probab = {'A':[0.2, 0.2, 0.3, 0.2, 0.3],
          'C':[0.4, 0.3, 0.1, 0.5, 0.1],
          'G':[0.3, 0.3, 0.5, 0.2, 0.4],
          'T':[0.1, 0.2, 0.1, 0.1, 0.2]}


def GreedyMotifSearch(Dna, k, t):
    BestMotifs = []
    for i in range(0, t):
        BestMotifs.append(Dna[i][0:k])

    n = len(Dna[0])
    for i in range(n-k+1):
        Motifs = []
        Motifs.append(Dna[0][i:i+k])
        for j in range(1, t):
            P = Profile(Motifs[0:j])
            Motifs.append(ProfileMostProbablePattern(Dna[j], k, P))

        if Score(Motifs) < Score(BestMotifs):
            BestMotifs = Motifs

    return BestMotifs


def CountWithPseudocounts(Motifs):
    count = {}
    k = len(Motifs[0])
    for symbol in "ACGT":
        count[symbol] = []
        for j in range(k):
            count[symbol].append(1)

    t = len(Motifs)
    for i in range(t):
        for j in range(k):
            symbol = Motifs[i][j]
            count[symbol][j] += 1
    return count


def ProfileWithPseudocounts(Motifs):
    counts = CountWithPseudocounts(Motifs)
    t = len(Motifs) + 4
    k = len(Motifs[0])
    profile = {}
    for symbol in "ACGT":
        profile[symbol] = []
        for j in range(len(counts[symbol])):
            profile[symbol].append((counts[symbol][j] / t))

    return profile


def GreedyMotifSearchWithPseudocounts(Dna, k, t):
    BestMotifs = []
    for i in range(0, t):
        BestMotifs.append(Dna[i][0:k])

    n = len(Dna[0])
    for i in range(n-k+1):
        Motifs = []
        Motifs.append(Dna[0][i:i+k])
        for j in range(1, t):
            P = ProfileWithPseudocounts(Motifs[0:j])
            Motifs.append(ProfileMostProbablePattern(Dna[j], k, P))

        if Score(Motifs) < Score(BestMotifs):
            BestMotifs = Motifs

    return BestMotifs


def Motifs(Profile, DNA):
    best_kmers = []
    t = len(DNA)

    for i in range(t):
        k = len(DNA[i])
        kmer = DNA[i][0:4]
        best_prob = Pr(kmer, Profile)
        for j in range(1,k-3):
            if Pr(DNA[i][j:j+4], Profile) > best_prob:
                best_prob = Pr(DNA[i][j:j+4], Profile)
                kmer = DNA[i][j:j+4]
        best_kmers.append(kmer)

    return best_kmers


ProfileTestCase0= { 'A': [0.8, 0.0, 0.0, 0.2 ],'C': [ 0.0, 0.6, 0.2, 0.0], 'G': [ 0.2 ,0.2 ,0.8, 0.0], 'T': [ 0.0, 0.2, 0.0, 0.8]}
DnaTC0=['TTACCTTAAC','GATGTCTGTC','ACGGCGTTAG','CCCTAACGAG','CGTCAGAGGT']

print(Motifs(ProfileTestCase0, DnaTC0))
