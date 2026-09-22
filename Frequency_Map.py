def FrequencyMap(Text, k):
    """
    This code gives a dictionary of each key (k-mer) along with how many times it occurs
    in a sequence
    param text (string): Text to be processed
    param k (int): number of kmers
    returns: (dict): Frequency of each key (k-mer) along with how many times it occurs
    """
    freq = {}
    n = len(Text)

    # REMOVE these — they overwrite user input
    # text = "ATCTAGGGATCTAATCGGAGATATATTAGTGACTCATACTGACAT"
    # k = 3

    for i in range(n - k + 1):
        pattern = Text[i:i + k]
        freq[pattern] = 0

    for i in range(n - k + 1):
        pattern = Text[i:i + k]
        freq[pattern] += 1

    return freq


def frequent_words(text, k):
    """
    This code is used to determine the most frequent kmers for a given sequence
    and frequency of each one
    :param text: (string): Text to be processed
    :param k: (int): Number of kmers
    :return: (list): Frequency of each kmer
    """

    list = []
    freq = FrequencyMap(text, k)
    m = max(freq.values())

    for key, value in freq.items():
        if value == m:
            list.append(key)

    return list


def PatternMatching(pattern, genome):
    """
    Finds all occurrences of a pattern in sequence
    input: strings pattern and genome
    output: all starting positions in genome where pattern appears as substring
    '''

    positions = []

    pattern_length = len(pattern)
    genome_length = len(genome)

    for i in range(genome_length - pattern_length + 1):
        if genome[i: i + pattern_length] == pattern:
            positions.append(i)
    return positions


def pattern_count(symbol, extended_genome):
    count = 0
    # output variable
    for i in range(len(extended_genome) - len(symbol) + 1):
        if extended_genome[i:i + len(symbol)] == symbol:
            count += 1
    return count


def faster_symbol_array(genome, symbol):
    array = {}
    n = len(genome)
    extended_genome = genome + genome[0:n // 2]

    # look at the first half of Genome to compute first array value
    array[0] = pattern_count(symbol, genome[0:n // 2])

    for i in range(1, n):
        # start by setting the current array value equal to the previous array value
        array[i] = array[i - 1]

        # the current array value can differ from the previous array value by at most 1
        if extended_genome[i - 1] == symbol:
            array[i] = array[i] - 1
        if extended_genome[i + (n // 2) - 1] == symbol:
            array[i] = array[i] + 1
    return array


import sys

lines = sys.stdin.read().splitlines()
print(faster_symbol_array(lines[0], lines[1]))


# ****running this code with the E. coli data found at
# https://bioinformaticsalgorithyms.com/data/realdatasets/Replication/E_coli.txt with
# symbol="C" results in the file being interpreted as bytes, not characters, so Python
# compares symbol "C" with b"A", b"G", b"T", & b"C", & because there is no match, it
# doesn't increment count in fxn PatternCount()

# recall in original PatternCount() we have:
# if patterns match, increment counter
# if (current_pattern == pattern):
#  count = count + 1

# This is corrected by changing the code to:
# if patterns match. increment counter
# current_Pattern = current_Pattern.decode()
# if current_Pattern == Pattern:
#  count = count + 1
