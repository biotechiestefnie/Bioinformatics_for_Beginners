
# Problem: finding the most conserved motif matrix - we score 1 for every non - conserved,
    # lower case value from the matrix.We arrange each k - mer then consider position
    # 0 to every other k - mer position 0, tally the number of unmatched bases.

# example input:
# ATTGC
# ATCGA
# Motifs: ACCGC
# ATCGG
# GGAGC

# This makes a 5x5 matrix, with the vertical tally, or output, being 1, 2, 2, 0, 2

# Score(Motifs)= 1 + 2 + 2 + 0 + 2 = 7

# Count(Motifs):
# A: 4, 0, 1, 0, 1
# T: 0, 3, 1, 0, 0
# C: 0, 1, 3, 0, 3
# G: 1, 1, 0, 4, 1

# Profile(Motifs)= each of the count motifs divided by the number of k-mers, so
# Profile(Motifs)= A: .8, 0, .2, 0, .2
T: 0, .6, .2, 0, 0
C: 0, .2, .6, 0, .6
G: .2, .2, 0, .8, .2


# Note that the total of each vertical column=1**** column 3 error?

# Finally, we make a cobsensus motif to determine the candidate regulatory motif for these
# sequences:

# Consensus(Motif): A, T, C, G, C
# Input:  A set of kmers Motifs
# Output: Count(Motifs)

def Count(Motifs):
    count = {}
    # initializing the count dictionary

    k = len(Motifs[0])
    for symbol in "ACGT":
        count[symbol] = []
        for j in range(k):
            count[symbol].append(0)

    # range over all nucleotides symbol and create a list of zeros corresponding to
    # count[symbol]
    t = len(Motifs)
    for i in range(t):
        for j in range(k):
            symbol = Motifs[i][j]
            count[symbol][j] += 1
    return count


# range over all elements symbol=Motifs[i][j] of the count matrix and add 1 to
# count[symbol][j]

def Profile(Motifs):
    t = len(Motifs)
    k = len(Motifs[0])
    profile = {}
    for symbol in "ATCG":
        profile[symbol] = []
        for j in range(k):
            profile[symbol].append(0)

    for i in range(t):
        for j in range(k):
            symbol = Motifs[i][j]
        profile[symbol][j] += 1 / k
    return profile


# all the same as Count(Motifs), but divide final values by k

# and finally,
def Consensus(Motifs):
    k = len(Motifs[0])
    count = Count(Motifs)
    consensus = ""
    for j in range(k):
        m = 0
        frequentSymbol = ""
        for symbol in "ACGT":
            if count[symbol][j] > m:
                m = count[symbol][j]
                frequentSymbol = symbol
        consensus += frequentSymbol
    return consensus


def Count(Motifs):
    count = {}
    # initializing the count dictionary

    # Set k equal to the length of Motifs[0]
    k = len(Motifs[0])

    # Iterate over all nucleotides symbol and create a list of zeroes corresponding
    # to count[symbol]
    for symbol in "ACGT":
        count[symbol] = []
        for j in range(k):
            count[symbol].append(0)

    # Iterate over all elements symbol = Motifs[i][j] of the count matrix and add
    # 1 to count[symbol][j]
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
        frequentSymbol = ""
        for symbol in "ACGT":
            if count[symbol][j] > m:
                m = count[symbol][j]
                frequentSymbol = symbol
        consensus += frequentSymbol
    return consensus
    # Insert your Count(Motifs) function here.


def Count(Motifs):
    count = {}  # initializing the count dictionary

    # Set k equal to the length of Motifs[0]
    k = len(Motifs[0])

    # Iterate over all nucleotides symbol and create a list of zeroes corresponding to count[symbol]
    for symbol in "ACGT":
        count[symbol] = []
        for j in range(k):
            count[symbol].append(0)

    # Iterate over all elements symbol = Motifs[i][j] of the count matrix and add 1 to count[symbol][j]
    t = len(Motifs)
    for i in range(t):
        for j in range(k):
            symbol = Motifs[i][j]
            count[symbol][j] += 1

    return count


# Input:  A set of k-mers Motifs
# Output: The score of these k-mers.

def Score(Motifs):
    consensus = Consensus(Motifs)

    # Get the consensus string using Consensus function
    k = len(consensus)
    t = len(Motifs)
    score = 0

    # Sum the number of symbols in the j-th column of Motifs that do not match the symbol
    # in position j of the consensus string
    for j in range(k):
        for i in range(t):
            if Motifs[i][j] != consensus[j]:
                score += 1

    return score


def Score(Motifs):
    consensus = Consensus(Motifs)
    k = len(consensus)
    t = len(Motifs)
    score = 0

    for j in range(k):
        for i in range(t):
            if Motifs[i][j] != consensus[j]:
                score += 1

    return score


def Count(Motifs):
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
    count = Count(Motifs)
    profile = {}

    k = len(Motifs[0])
    t = len(Motifs)

    for symbol in "ACGT":
        profile[symbol] = []
        for j in range(k):
            profile[symbol].append(count[symbol][j] / t)

    return profile


def Consensus(Motifs):
    profile = Profile(Motifs)
    consensus = ""

    k = len(Motifs[0])
    for j in range(k):
        m = 0
        frequentSymbol = ""
        for symbol in "ACGT":
            if profile[symbol][j] > m:
                m = profile[symbol][j]
                frequentSymbol = symbol
        consensus += frequentSymbol

    return consensus


def ProfileMostProbableKmer(Text, k, Profile):
    max_prob = -1
    most_probable = ""

    for i in range(len(Text) - k + 1):
        pattern = Text[i:i + k]
        prob = Pr(pattern, Profile)
        if prob > max_prob:
            max_prob = prob
            most_probable = pattern

    return most_probable


def Pr(Text, Profile):
    p = 1.0

    for i in range(len(Text)):
        symbol = Text[i]
        column = i

        p *= Profile[symbol][column]

    return p


def GreedyMotifSearch(Dna, k, t):
    n = len(Dna[0])
    best_motifs = [Dna[i][0:k] for i in range(t)]

    for i in range(n - k + 1):
        motifs = [Dna[0][i:i + k]]
        for j in range(1, t):
            profile = Profile(motifs[0:j])
            motifs.append(ProfileMostProbableKmer(Dna[j], k, profile))

        if Score(motifs) < Score(best_motifs):
            best_motifs = motifs

    return best_motifs


 if __name__ == "__main__":
     k, t = map(int, input().split())
     dna_input = [input().strip() for _ in range(t)]
     result = GreedyMotifSearch(dna_input, k, t)
     print("\n".join(result))


# For LaPlace's Law, Write a function CountWithPseudocounts(Motifs)
# that takes a list of strings Motifs as input and returns the count matrix
# of Motifs with pseudocounts as a dictionary of lists:
# Input:  A set of kmers Motifs
# Input:  A set of kmers Motifs
# Output: CountWithPseudocounts(Motifs)
def CountWithPseudocounts(Motifs):
    t = len(Motifs)
    # number of k-mers
    k = len(Motifs[0])
    # length of each k-mer
    # Initialize the count matrix with pseudocounts
    counts = {'A': [1] * k, 'C': [1] * k, 'G': [1] * k, 'T': [1] * k}

    # Iterate through each position in the motifs
    for i in range(k):
        for motif in Motifs:
            counts[motif[i]][i] += 1

    return counts


# write a function ProfileWithPseudocounts(Motifs) that takes a list of strings Motifs
# as input and returns the profile matrix of Motifs with pseudocounts as a dictionary
# of lists.
def ProfileWithPseudocounts(Motifs):
    t = len(Motifs)
    # Number of k-mers
    k = len(Motifs[0])
    # Length of each k-mer
    profile = {}
    # Initialize an empty dictionary for the profile

    # Compute the count matrix with pseudocounts
    count = CountWithPseudocounts(Motifs)

    # Calculate the probabilities for each base at each position
    for base in "ACGT":
        profile[base] = [count[base][i] / (t + 4) for i in range(k)]

    return profile


def CountWithPseudocounts(Motifs):
    t = len(Motifs)
    # Number of k-mers
    k = len(Motifs[0])
    # Length of each k-mer
    count = {}
    # Initialize an empty dictionary for counts

    # Initialize the count matrix with pseudocounts
    for base in "ACGT":
        count[base] = [1] * k

    # Iterate through each position in the k-mers
    for i in range(k):
        for j in range(t):
            base = Motifs[j][i]
            count[base][i] += 1

    return count


# Write a function GreedyMotifSearchWithPseudocounts(Dna, k, t)
# that takes as input a list of strings Dna followed by integers
# k and t and returns the result of running GreedyMotifSearch(),
# where each profile matrix is generated with pseudocounts:

# To solve this code challenge, we need to modify the original GreedyMotifSearch()
# algorithm to incorporate pseudocounts. This means we will use CountWithPseudocounts(Motifs)
# and ProfileWithPseudocounts(Motifs) instead of the regular Count(Motifs) and Profile(Motifs).

# Steps to Approach:

# 1. GreedyMotifSearchWithPseudocounts(Dna, k, t):

# This function will run the GreedyMotifSearch() algorithm but with pseudocounts applied
# to the profile matrix.

# We will update the motif finding part by using the new versions of the Count() and
# Profile() functions.

# 2. Use of Pseudocounts:

# Instead of using the regular count and profile matrices, we will use CountWithPseudocounts
# to ensure that there are no zero probabilities for motifs that haven't been observed yet.
# This will help in calculating more accurate probabilities, especially for rare or unseen
# motifs.

# 3. Implementing the GreedyMotifSearch with Pseudocounts:

# Follow the original GreedyMotifSearch() algorithm, but modify it to use CountWithPseudocounts
# and ProfileWithPseudocounts in each iteration.


# Here is the code for the function:

# CountWithPseudocounts(Motifs): Count function using pseudocounts (Laplace's Rule of Succession)
def CountWithPseudocounts(Motifs):
    t = len(Motifs)
    # number of strings
    k = len(Motifs[0])
    # length of each motif (k-mer)

    # Initialize count matrix with pseudocounts
    count = {'A': [1] * k, 'C': [1] * k, 'G': [1] * k, 'T': [1] * k}

    for motif in Motifs:
        for i in range(k):
            count[motif[i]][i] += 1

    return count


# ProfileWithPseudocounts(Motifs): Profile matrix using pseudocounts
def ProfileWithPseudocounts(Motifs):
    t = len(Motifs)
    count = CountWithPseudocounts(Motifs)
    # Get the count matrix using pseudocounts

    # Calculate profile by dividing each count by the total (with pseudocounts)
    profile = {}
    for symbol in count:
        profile[symbol] = [count[symbol][i] / float(t + 4) for i in range(len(count[symbol]))]  # +4 for pseudocounts

    return profile


# GreedyMotifSearchWithPseudocounts(Dna, k, t): Greedy Motif Search using pseudocounts
def GreedyMotifSearchWithPseudocounts(Dna, k, t):
    BestMotifs = []

    # Step 1: Initialize BestMotifs with the first k-mer from each string
    BestMotifs = [Dna[i][:k] for i in range(t)]
    # Take the first k-mer from each DNA string

    # Step 2: Greedy search for the best motif starting with the first k-mer
    for i in range(t):
        Motifs = [BestMotifs[0]]
        # Start with the first motif
        for j in range(1, t):
            profile = ProfileWithPseudocounts(Motifs)
            # Calculate profile with pseudocounts
            best_motif = MostProbableKmer(Dna[j], k, profile)
            # Find the most probable k-mer in the j-th string
            Motifs.append(best_motif)
        # Add this motif to the Motifs list

        # Calculate the score of the current Motifs
        current_score = Score(Motifs)
        # Score function is usually defined by sum of column-wise entropy or other scoring method
        best_score = Score(BestMotifs)

        if current_score < best_score:
            BestMotifs = Motifs
        # Update BestMotifs if the current score is better (lower score)

    return BestMotifs


# Function to find the most probable k-mer in a string based on a profile matrix
def MostProbableKmer(Dna, k, profile):
    max_prob = -1
    best_kmer = None
    for i in range(len(Dna) - k + 1):
        kmer = Dna[i:i + k]
        prob = 1.0
        for j in range(k):
            prob *= profile[kmer[j]][j]
        if prob > max_prob:
            max_prob = prob
            best_kmer = kmer
    return best_kmer


# Function to compute the score of a set of motifs (sum of column-wise entropy)
def Score(Motifs):
    k = len(Motifs[0])
    t = len(Motifs)
    score = 0
    for j in range(k):
        counts = {'A': 0, 'C': 0, 'G': 0, 'T': 0}
        for i in range(t):
            counts[Motifs[i][j]] += 1
        most_common = max(counts.values())
        score += (t - most_common)
    return score


# Explanation of Functions:

# 1. CountWithPseudocounts(Motifs):

# Counts the occurrences of each base at each position across all motifs,
# with pseudocounts (i.e., we add 1 to each count to avoid zero probabilities).

# 2. ProfileWithPseudocounts(Motifs):

# Uses the count matrix to calculate a profile matrix by normalizing the counts,
# accounting for pseudocounts.

# 3. GreedyMotifSearchWithPseudocounts(Dna, k, t):

# Starts with the first k-mer from each DNA string.

# Iteratively builds the motifs by selecting the most probable k-mer
# at each position based on the profile matrix with pseudocounts.

# Updates the best motifs based on the motif score.

# 4. MostProbableKmer(Dna, k, profile):

# Given a DNA string and a profile matrix, finds the most probable k-mer
# by calculating the probability of each possible k-mer.

# 5. Score(Motifs):

# Computes a score for a set of motifs. A lower score indicates better motifs.


# Example:

# For the sample input:

Dna = ["GGCGTTCAGGCA", "AAGAATCAGTCA", "CAAGGAGTTCGC", "CACGTCAATCAC", "CAATAATATTCG"]
k = 3
t = 5

# The output would be:

TTC
ATC
TTC
ATC
TTC


# How It Works:

# The function starts by selecting the first k-mer from each DNA string as the initial
# motifs.

# Then, for each subsequent DNA string, it calculates the most probable k-mer based on
# the profile matrix (with pseudocounts).

# It iteratively refines the motif set, always choosing the k-mer with the highest
# probability in the next DNA string.

# The best motifs are chosen based on the motif set with the lowest score, which is
# determined by the motif's alignment and how close it is to the consensus string.


# This approach allows GreedyMotifSearch() to work better, especially with rare motifs,
# by ensuring that even unseen motifs are assigned a small probability due to the
# pseudocounts.


# Apply GreedyMotifSearchWithPseudocounts() to find motifs in the DosR dataset with k
# equal to 15


def Score(Motifs):
    consensus = Consensus(Motifs)
    k = len(consensus)
    t = len(Motifs)
    score = 0

    for j in range(k):
        for i in range(t):
            if Motifs[i][j] != consensus[j]:
                score += 1

    return score


def Count(Motifs):
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
    count = Count(Motifs)
    profile = {}

    k = len(Motifs[0])
    t = len(Motifs)

    for symbol in "ACGT":
        profile[symbol] = []
        for j in range(k):
            profile[symbol].append(count[symbol][j] / t)

    return profile


def Consensus(Motifs):
    profile = Profile(Motifs)
    consensus = ""

    k = len(Motifs[0])
    for j in range(k):
        m = 0
        frequentSymbol = ""
        for symbol in "ACGT":
            if profile[symbol][j] > m:
                m = profile[symbol][j]
                frequentSymbol = symbol
        consensus += frequentSymbol

    return consensus


def CountWithPseudocounts(Motifs):
    count = {symbol: [1] * len(Motifs[0]) for symbol in "ACGT"}
    for i in range(len(Motifs)):
        for j in range(len(Motifs[i])):
            symbol = Motifs[i][j]
            count[symbol][j] += 1
    return count


def ProfileWithPseudocounts(Motifs):
    t = len(Motifs) + 4
    count = CountWithPseudocounts(Motifs)
    profile = {symbol: [float(c) / t for c in count[symbol]] for symbol in "ACGT"}
    return profile


def ScoreWithPseudocounts(Motifs):
    count = CountWithPseudocounts(Motifs)
    score = 0
    for j in range(len(Motifs[0])):
        max_count = max(count[symbol][j] for symbol in "ACGT")
        score += sum(count[symbol][j] for symbol in "ACGT") - max_count
    return score


def ProfileMostProbablePattern(Text, k, profile):
    max_prob = -1
    most_probable = ""
    for i in range(len(Text) - k + 1):
        pattern = Text[i:i + k]
        prob = 1
        for j in range(k):
            prob *= profile[pattern[j]][j]
        if prob > max_prob:
            max_prob = prob
            most_probable = pattern
    return most_probable


def GreedyMotifSearchWithPseudocounts(Dna, k, t):
    BestMotifs = [string[:k] for string in Dna]

    for i in range(len(Dna[0]) - k + 1):
        Motifs = [Dna[0][i:i + k]]

        for j in range(1, t):
            profile = ProfileWithPseudocounts(Motifs)
            motif = ProfileMostProbablePattern(Dna[j], k, profile)
            Motifs.append(motif)

        if ScoreWithPseudocounts(Motifs) < ScoreWithPseudocounts(BestMotifs):
            BestMotifs = Motifs

    return BestMotifs


# The ten strings occurring in the hyperlinked DosR dataset below.
Dna = [
    "GCGCCCCGCCCGGACAGCCATGCGCTAACCCTGGCTTCGATGGCGCCGGCTCAGTTAGGGCCGGAAGTCCCCAATGTGGCAGACCTTTCGCCCCTGGCGGACGAATGACCCCAGTGGCCGGGACTTCAGGCCCTATCGGAGGGCTCCGGCGCGGTGGTCGGATTTGTCTGTGGAGGTTACACCCCAATCGCAAGGATGCATTATGACCAGCGAGCTGAGCCTGGTCGCCACTGGAAAGGGGAGCAACATC",
    "CCGATCGGCATCACTATCGGTCCTGCGGCCGCCCATAGCGCTATATCCGGCTGGTGAAATCAATTGACAACCTTCGACTTTGAGGTGGCCTACGGCGAGGACAAGCCAGGCAAGCCAGCTGCCTCAACGCGCGCCAGTACGGGTCCATCGACCCGCGGCCCACGGGTCAAACGACCCTAGTGTTCGCTACGACGTGGTCGTACCTTCGGCAGCAGATCAGCAATAGCACCCCGACTCGAGGAGGATCCCG",
    "ACCGTCGATGTGCCCGGTCGCGCCGCGTCCACCTCGGTCATCGACCCCACGATGAGGACGCCATCGGCCGCGACCAAGCCCCGTGAAACTCTGACGGCGTGCTGGCCGGGCTGCGGCACCTGATCACCTTAGGGCACTTGGGCCACCACAACGGGCCGCCGGTCTCGACAGTGGCCACCACCACACAGGTGACTTCCGGCGGGACGTAAGTCCCTAACGCGTCGTTCCGCACGCGGTTAGCTTTGCTGCC",
    "GGGTCAGGTATATTTATCGCACACTTGGGCACATGACACACAAGCGCCAGAATCCCGGACCGAACCGAGCACCGTGGGTGGGCAGCCTCCATACAGCGATGACCTGATCGATCATCGGCCAGGGCGCCGGGCTTCCAACCGTGGCCGTCTCAGTACCCAGCCTCATTGACCCTTCGACGCATCCACTGCGCGTAAGTCGGCTCAACCCTTTCAAACCGCTGGATTACCGACCGCAGAAAGGGGGCAGGAC",
    "GTAGGTCAAACCGGGTGTACATACCCGCTCAATCGCCCAGCACTTCGGGCAGATCACCGGGTTTCCCCGGTATCACCAATACTGCCACCAAACACAGCAGGCGGGAAGGGGCGAAAGTCCCTTATCCGACAATAAAACTTCGCTTGTTCGACGCCCGGTTCACCCGATATGCACGGCGCCCAGCCATTCGTGACCGACGTCCCCAGCCCCAAGGCCGAACGACCCTAGGAGCCACGAGCAATTCACAGCG",
    "CCGCTGGCGACGCTGTTCGCCGGCAGCGTGCGTGACGACTTCGAGCTGCCCGACTACACCTGGTGACCACCGCCGACGGGCACCTCTCCGCCAGGTAGGCACGGTTTGTCGCCGGCAATGTGACCTTTGGGCGCGGTCTTGAGGACCTTCGGCCCCACCCACGAGGCCGCCGCCGGCCGATCGTATGACGTGCAATGTACGCCATAGGGTGCGTGTTACGGCGATTACCTGAAGGCGGCGGTGGTCCGGA",
    "GGCCAACTGCACCGCGCTCTTGATGACATCGGTGGTCACCATGGTGTCCGGCATGATCAACCTCCGCTGTTCGATATCACCCCGATCTTTCTGAACGGCGGTTGGCAGACAACAGGGTCAATGGTCCCCAAGTGGATCACCGACGGGCGCGGACAAATGGCCCGCGCTTCGGGGACTTCTGTCCCTAGCCCTGGCCACGATGGGCTGGTCGGATCAAAGGCATCCGTTTCCATCGATTAGGAGGCATCAA",
    "GTACATGTCCAGAGCGAGCCTCAGCTTCTGCGCAGCGACGGAAACTGCCACACTCAAAGCCTACTGGGCGCACGTGTGGCAACGAGTCGATCCACACGAAATGCCGCCGTTGGGCCGCGGACTAGCCGAATTTTCCGGGTGGTGACACAGCCCACATTTGGCATGGGACTTTCGGCCCTGTCCGCGTCCGTGTCGGCCAGACAAGCTTTGGGCATTGGCCACAATCGGGCCACAATCGAAAGCCGAGCAG",
    "GGCAGCTGTCGGCAACTGTAAGCCATTTCTGGGACTTTGCTGTGAAAAGCTGGGCGATGGTTGTGGACCTGGACGAGCCACCCGTGCGATAGGTGAGATTCATTCTCGCCCTGACGGGTTGCGTCTGTCATCGGTCGATAAGGACTAACGGCCCTCAGGTGGGGACCAACGCCCCTGGGAGATAGCGGTCCCCGCCAGTAACGTACCGCTGAACCGACGGGATGTATCCGCCCCAGCGAAGGAGACGGCG",
    "TCAGCACCATGACCGCCTGGCCACCAATCGCCCGTAACAAGCGGGACGTCCGCGACGACGCGTGCGCTAGCGCCGTGGCGGTGACAACGACCAGATATGGTCCGAGCACGCGGGCGAACCTCGTGTTCTGGCCTCGGCCAGTTGTGTAGAGCTCATCGCTGTCATCGAGCGATATCCGACCACTGATCCAAGTCGGGGGCTCTGGGGACCGAAGTCCCCGGGCTCGGAGCTATCGGACCTCACGATCACC"]

# set t equal to the number of strings in Dna and k equal to 15
t = len(Dna)
k = 15

# Call GreedyMotifSearchWithPseudocounts(Dna, k, t) and store the output in a variable
# called Motifs
Motifs = GreedyMotifSearchWithPseudocounts(Dna, k, t)

# Print the Motifs variable
print(Motifs)
# Print Score(Motifs)
print(Score(Motifs))

------------------------------------------------------------------------------
# You have now seen the power of pseudocounts illustrated on a small example.
# Running GreedyMotifSearchWithPseudocounts() to solve the Subtle Motif Problem
# returns a collection of 15-mers Motifs with Score(Motifs) equal to 41 and
# Consensus(Motifs) equal to "AAAAAtAgaGGGGtt". Thus, Laplace’s Rule of Succession
# has provided a significant improvement over the original GreedyMotifSearch(),
# which returned the consensus string "gttAAAtAgaGatGtG" with Score(Motifs) equal to 58.
-----------------------------------------------------------------------------------------


# A function, Motifs(Profile, Dna) that takes as input a profile matrix Profile
# corresponding to a list of strings Dna and that returns a list of the profile-most
# probable k-mers in each string from Dna, using Profile:

# Input: A profile matrix Profile and a list of strings Dna
# Output: Motifs(Profile, Dna)

def ProfileMostProbablePattern(Text, k, profile):
    max_prob = -1
    most_probable = ""
    for i in range(len(Text) - k + 1):
        pattern = Text[i:i + k]
        prob = 1
        for j in range(k):
            prob *= profile[pattern[j]][j]
        if prob > max_prob:
            max_prob = prob
            most_probable = pattern
    return most_probable


def Motifs(Profile, Dna):
    k = len(Profile["A"])
    most_probable_kmers = []
    for text in Dna:
        most_probable = ProfileMostProbablePattern(text, k, Profile)
        most_probable_kmers.append(most_probable)
    return most_probable_kmers


# Test input
profile_input = [
    [0.8, 0.0, 0.0, 0.2],
    [0.0, 0.6, 0.2, 0.0],
    [0.2, 0.2, 0.8, 0.0],
    [0.0, 0.2, 0.0, 0.8]
]

Dna_input = [
    "TTACCTTAAC",
    "GATGTCTGTC",
    "ACGGCGTTAG",
    "CCCTAACGAG",
    "CGTCAGAGGT"
]

profile = {"A": profile_input[0], "C": profile_input[1], "G": profile_input[2], "T": profile_input[3]}

result = Motifs(profile, Dna_input)

-------------------------------------------------------------------------------------
# Simulating the process of generating a random integer is difficult and takes a lot of math.
# Python provides a module called random for generating them. You can think of a module as
# a “bundle” of related functions. To use the random module, we place the following statement
# at the top of our file.

import random

# Inside of the random module is a built-in function called randint(1, M) that
# generates a random integer between 1 and M, inclusively. To call this function,
# we use

random.randint(1, M)
------------------------------------------------------------------------------------------------

# A function, RandomMotifs(Dna, k, t) that uses random.randint() to choose a random
# k-mer from each of t different strings Dna, and returns a list of t strings:

# We start by importing the random module:
import random


# Input:  A list of strings Dna, and integers k and t
# Output: RandomMotifs(Dna, k, t)
def RandomMotifs(Dna, k, t):
    motifs = []
    for i in range(t):
        string = Dna[i]
        random_position = random.randint(0, len(string) - k)
        random_kmer = string[random_position: random_position + k]
        motifs.append(random_kmer)
    return motifs

# Now, we can develop RandomizedMotifSearch().
# The code below stops running as soon as the score of the motifs that we
# generate stops improving. It uses the loop “while True”, which iterates until
# it encounters a return statement. It can be dangerous to use such a loop, since
# it could lead to an infinite loop in which a program never terminates. However,
# in this particular case, the motif score must eventually stop improving, so that
# RandomizedMotifSearch() must eventually terminate.

