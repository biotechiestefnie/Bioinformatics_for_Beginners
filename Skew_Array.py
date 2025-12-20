# The function skew_array(Genome) below takes a sequence as input and returns the skew
# array of Genome in the form of a list whose i-th element is skew[i]. We use this to
# determine the count of #G-#C at position i in Genome:

def skew_array(genome):
    skew = [0]
    n = len(genome)
    for i in range(n):
        if genome == "C":
            skew.append(skew[i] - 1)
        elif genome[i] == "G":
            skew.append(skew[i] + 1)
        else:
            skew.append(skew[i])
    return skew

# This next function returns the minimum values in the skew diagram:
# Input: A DNA string Genome
# Output: All integer(s) i minimizing skew[i] among all values of i (from 0 to len(Genome)

def minimum_skew(genome):
    positions = []
    # list to store positions with minimum skew
    skew = skew_array(genome)
    # Call the SkewArray function to calculate skew
    m = min(skew)
    # Find the min skew value(s) in the list
    for i in range(len(genome)):
        # Iterate through the skew list
        if skew[i] == m:
            # Check if current skew equals min
            positions.append(i)
    # Add the position to the results
    return positions
