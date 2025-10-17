import math
import matplotlib.pyplot as plt

def read_fasta(filename):
    with open(filename, 'r') as file:
        lines = file.readlines()
    sequence = ''.join(line.strip() for line in lines if not line.startswith('>'))
    return sequence.upper()

def temperature_formula1(G, C, A, T):
    temp = 4 * (G + C) + 2 * (A + T)
    return temp


def temperature_formula2(G, C, seq_length, Na=0.001):
    temp = -(81.5 + 16.6 * math.log10(Na) + 0.41 * (((G + C) / seq_length) * 100) - (600 / seq_length))
    return temp

def sliding_window_temperature1(sequence, window_size=9):
    P = []
    for i in range(len(sequence) - window_size + 1):
        window = sequence[i:i + window_size]
        G = window.count("G")
        C = window.count("C")
        A = window.count("A")
        T = window.count("T")
        temp = temperature_formula1(G, C, A, T)   
        P.append(temp)
    return P

def sliding_window_temperature2(sequence, window_size=9):
    P = []
    for i in range(len(sequence) - window_size + 1):
        window = sequence[i:i + window_size]
        G = window.count("G")
        C = window.count("C")
        temp = temperature_formula2(G, C, len(window))
        P.append(temp)
    return P


def plot_vectors(P1, P2):
    plt.figure(figsize=(10, 5))
    plt.plot(range(1, len(P1) + 1), P1, color='darkred', linewidth=1.5, label='Formula 1 (Simple)')
    plt.plot(range(1, len(P2) + 1), P2, color='darkblue', linewidth=1.5, label='Formula 2 (Salt-adjusted)')
    plt.title("Sliding Window DNA Melting Temperature (Window = 9)")
    plt.xlabel("Sequence Position")
    plt.ylabel("Temperature (°C)")
    plt.grid(True, linestyle='--', alpha=0.6)
    plt.legend()
    plt.tight_layout()
    plt.show()


def main():
    fasta_file = "C:\LaylaYear4\BioInformatics\BioInformatics\Project_L3\L3\dna.fasta"   
    sequence = read_fasta(fasta_file)
    P1 = sliding_window_temperature1(sequence, window_size=9)
    P2 = sliding_window_temperature2(sequence, window_size=9)
    plot_vectors(P1, P2)


if __name__ == "__main__":
    main()
