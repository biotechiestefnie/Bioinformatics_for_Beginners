def reverse_complement(pattern):
    """
    This code gives us the reverse complement of a given pattern
    :param pattern: (string) the pattern to be reversed
    :return: (string) the reverse complement
    """
    rev = reverse(pattern)
    comp = complement(rev)
    pattern =
    return comp


def Complement(rev):
    comp = ""
    for char in pattern:
        if char == "A":
            comp += "T"
        elif char == "T":
            comp += "A"
        elif char == "C":
            comp += "G"
        elif char == "G":
            comp += "C"
    return comp

