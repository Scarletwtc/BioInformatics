#make an artificial DNA sequence of 200-400 bases in lenght in which to simulate
#3-4 transposable element
# start= ATGCA end = TACGT
import random

# Function to generate a random DNA sequence of given length
def generate_dna_sequence(length):
    return ''.join(random.choice('ATCG') for _ in range(length))

# Function to generate a transposable element
def generate_transposable_element(min_len=5, max_len=20):
    internal_length = random.randint(min_len, max_len)
    internal_seq = ''.join(random.choice('ATCG') for _ in range(internal_length))
    return 'ATGCA' + internal_seq + 'TACGT'

# Function to insert TEs at random positions with * markers
def insert_transposable_elements(dna_seq, num_tes):
    dna_list = list(dna_seq)
    for _ in range(num_tes):
        te = generate_transposable_element()
        te_marked = '*' + te + '*'  # mark the TE
        pos = random.randint(0, len(dna_list))
        dna_list.insert(pos, te_marked)
    return ''.join(dna_list)

# Generate a DNA sequence of 200-400 bases
seq_length = random.randint(200, 400)
dna_sequence = generate_dna_sequence(seq_length)

# Insert 3-4 transposable elements
num_tes_to_insert = random.randint(3, 4)
artificial_sequence = insert_transposable_elements(dna_sequence, num_tes_to_insert)

# Output
print("Original DNA length:", seq_length)
print("Number of TEs inserted:", num_tes_to_insert)
print("\nArtificial DNA sequence (TEs marked with *):\n")
print(artificial_sequence)
print("\nFinal length:", len(artificial_sequence))
