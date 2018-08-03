# Input:  A string Text and an integer k
# Output: A list containing all most frequent k-mers in Text
def FrequentWords(Text, k):
    # your code here
    words = []
    freq = FrequencyMap(Text, k)
    m = max(freq.values())
    for key in freq:
        # add each key to words whose corresponding frequency value is equal to m
        if(freq[key] == m):
            words.append(key)
    return words

# Copy your FrequencyMap() function here.
def FrequencyMap(Text, k):
    # your code here.
    freq = {}
    n = len(Text)
    for i in range(n-k+1):
        Pattern = Text[i:i+k]
        freq[Pattern] = 0
    # hint: your code goes here!
    for i in range(n-k+1):
        Pattern = Text[i:i+k]
        freq[Pattern] += 1
    return freq

# Input:  A string Pattern
# Output: The reverse of Pattern
def Reverse(Pattern):
    rev = ""
    for i in range(len(Pattern)):
        rev = Pattern[i] + rev
    return rev

# Input:  A DNA string Pattern
# Output: The complementary string of Pattern (with every nucleotide replaced by its complement).
def Complement(Pattern):
    # your code here
    comp = ""
    for i in range(len(Pattern)):
        if(Pattern[i] == 'A'):
            comp += 'T'
        elif(Pattern[i] == 'T'):
            comp += 'A'
        elif(Pattern[i] == 'G'):
            comp += 'C'
        elif(Pattern[i] == 'C'):
            comp += 'G'
    return comp


# Input:  A DNA string Pattern
# Output: The reverse complement of Pattern
def ReverseComplement(Pattern):
    # your code here
    reverse = Reverse(Pattern)
    reverseComplement = Complement(reverse)
    return reverseComplement


# fill in your PatternMatching() function along with any subroutines that you need.
def PatternMatching(Genome, Pattern):
    positions = [] # output variable
    # your code here
    for i in range(len(Genome) - len(Pattern) +1):
        if (Genome[i:i+len(Pattern)] == Pattern):
            positions.append(i)
    return positions

def PatternCount(Text, Pattern):
    count = 0
    for i in range(len(Text)-len(Pattern)+1):
        if Text[i:i+len(Pattern)] == Pattern:
            count = count+1
    return count

def symbolArray(genome, symbol):
array = {}
    # your code here
    n = len(genome)
    extendedGenome = genome + genome[0:n//2]
    # look at the first half of Genome to compute first array value
    array[0] = PatternCount(symbol, genome[0:n//2])

    for i in range(1, n):
        # start by setting the current array value equal to the previous array value
        array[i] = array[i-1]

        # the current array value can differ from the previous array value by at most 1
        if extendedGenome[i-1] == symbol:
            array[i] = array[i]-1
        if extendedGenome[i+(n//2)-1] == symbol:
            array[i] = array[i]+1
    return array

genome = input()
print(symbolArray(genome, 'C'))
