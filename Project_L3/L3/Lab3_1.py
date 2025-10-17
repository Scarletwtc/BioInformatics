def temperature_formula1(G,C,A,T):
    temp = 4*(G+C) + 2*(A+T)
    return temp

S="ATTTCGCCGATA"

G = S.count("G")
C = S.count("C")
A = S.count("A")
T = S.count("T")

print("Simple formula: " + str(temperature_formula1(G,C,A,T))+ " °C")

Na = 0.001
import math 

def temperature_formula2(G,C):
    temp = 81.5 + 16.6*math.log10(Na)+ 0.41*(((G+C)/len(S))*100)-600/len(S)
    return temp

print("Salt adjusted formula: " + str(temperature_formula2(G,C)) + " °C")

