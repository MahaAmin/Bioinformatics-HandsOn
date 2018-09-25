import random
from operator import itemgetter

def RandomMotifs(DNA, k, t):
    randomKmers = []
    for i in range(t):
        start = random.randint(1, len(DNA[i]) - k - 1)
        kmer = DNA[i][start:start+k]
        randomKmers.append(kmer)
    return randomKmers


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


def Pr(Text, Profile):
    prob = 1
    for i in range(len(Text)):
        prob = prob * Profile[Text[i]][i]
    return prob


def Normalize(Probabilities):
    k = len(Probabilities)
    sumOfKmers = 0
    for key in Probabilities.keys():
        sumOfKmers += Probabilities[key]
    for key in Probabilities.keys():
        Probabilities[key] = Probabilities[key] / sumOfKmers

    return Probabilities


def WeightedDie(Probabilities):
    kmer = ''
    p = random.uniform(0,1)
    prob = 0
    for i,j in sorted(Probabilities.items(), key=itemgetter(1)):
        if prob <= p < j + prob:
            return i
        prob += j


def ProfileGeneratedString(Text, profile, k):
        n = len(Text)
        probabilities = {}
        for i in range(0,n-k+1):
            probabilities[Text[i:i+k]] = Pr(Text[i:i+k], profile)
        probabilities = Normalize(probabilities)
        return WeightedDie(probabilities)


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


def Score(Motifs):
    consensus = Consensus(Motifs)
    k = len(Motifs)
    score = 0
    for i in range(len(Motifs[0])):
        for j in range(k):
            if(Motifs[j][i] != consensus[i]):
                score += 1
    return score


def GibbsSampler(Dna, k, t, N):
    Motifs = RandomMotifs(Dna, k, t)
    BestMotifs = Motifs

    for j in range(1, N):
        i = random.randint(0, t-1)
        motifsExceptOne = []
        for counter in range(t):
            if counter == i :
                ignoredString = Dna[i]
            else:
                motifsExceptOne.append(Motifs[counter])

        profile = ProfileWithPseudocounts(motifsExceptOne)

        Motifs[i] = ProfileGeneratedString(ignoredString, profile, k)

        if Score(Motifs) < Score(BestMotifs):
            BestMotifs = Motifs
    return BestMotifs


DNA = ["CGCCCCTCTCGGGGGTGTTCAGTAAACGGCCA", "GGGCGAGGTATGTGTAAGTGCCAAGGTGCCAG",
"TAGTACCGAGACCGAAAGAAGTATACAGGCGT", "TAGATCAAGTTTCAGGTGCACGTCGGTGAACC", "AATCCACCAGCTCCACGTGCAATGTTGGCCTA"]

print(GibbsSampler(DNA, 8, 5, 100))
