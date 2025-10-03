def relativeFrequencies(s):
    frequencies = {}
    total = len(s)
    
    for letter in s:
        count = s.count(letter)
        frequencies[letter] = count / total * 100  #percent
    
    return frequencies

S="ATTTCGCCGATA"

for ch, freq in relativeFrequencies(S).items():
    print(f"{ch}: {freq:.2f}%")