def PatternCount(Text, Pattern):
    """
    This code gives the number of times a Pattern (k-mer) occurs in a string Text,
    or given sequence
    :param Text:
    :param Pattern:
    :return: Count of occurences of Pattern in Text
    """

    count = 0
    for i in range(len(Text) - len(Pattern) + 1):
        if Text[i:i + len(Pattern)] == Pattern:
            count = count + 1
    return count



