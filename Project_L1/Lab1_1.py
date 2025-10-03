#design app to find alphabet is sequence s = "ATTTCGCCGATA"
alphabet =[]
S="ATTTCGCCGATA"
def findAlphabet(s):
    for char in s:
        charExists = False
        for letter in alphabet:
            if char == letter:
                charExists = True
                break
        if not charExists:
            alphabet.append(char)
    return alphabet

print(findAlphabet(S))