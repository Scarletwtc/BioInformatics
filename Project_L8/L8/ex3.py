from pathlib import Path
import matplotlib.pyplot as plt

# -------------------------------------------------------------------
# Basic FASTA reader
# -------------------------------------------------------------------
def read_fasta(filepath):
    """
    Read a FASTA file and return a single concatenated uppercase sequence.
    Ignores header lines starting with '>'.
    """
    seq_parts = []
    with open(filepath, 'r') as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith(">"):
                continue
            seq_parts.append(line)
    return "".join(seq_parts).upper()


# -------------------------------------------------------------------
# Reverse complement (for inverted repeats)
# -------------------------------------------------------------------
def reverse_complement(seq):
    """
    Return the reverse complement of a DNA sequence.
    If you ONLY want the complement (like ATGCA -> TACGT),
    remove the [::-1] at the end.
    """
    comp_table = str.maketrans("ACGTNacgtn", "TGCANtgcan")
    return seq.translate(comp_table)[::-1]


# -------------------------------------------------------------------
# Main function to find inverted repeats
# -------------------------------------------------------------------
def find_inverted_repeats(seq, min_len=4, max_len=6, max_spacer=None):
    """
    Find pairs of inverted repeats in sequence `seq`.

    - min_len, max_len: length range for the repeat (e.g. 4..6).
    - max_spacer: maximum allowed distance (in bp) between the left and right repeats.
      If None, no upper bound is enforced.

    Returns a list of dicts describing each pair.
    """
    seq = seq.upper()
    n = len(seq)
    results = []

    for i in range(n - min_len + 1):
        # Sliding window for the left repeat
        for L in range(min_len, max_len + 1):
            if i + L > n:
                break

            left = seq[i:i + L]

            # Skip ambiguous repeats (you can remove this if you want Ns included)
            if any(b not in "ACGT" for b in left):
                continue

            # "Invert" the left window: reverse complement
            right_pattern = reverse_complement(left)

            # Second sliding window: search for the inverted pattern to the right
            start_search = i + L   # don’t allow overlapping left & right repeats
            j = seq.find(right_pattern, start_search)

            while j != -1:
                spacer_len = j - (i + L)

                if max_spacer is None or spacer_len <= max_spacer:
                    results.append({
                        "left_start_0": i,
                        "left_end_0": i + L - 1,
                        "right_start_0": j,
                        "right_end_0": j + L - 1,
                        "left_start_1": i + 1,           # 1-based coords
                        "left_end_1": i + L,
                        "right_start_1": j + 1,
                        "right_end_1": j + L,
                        "repeat_len": L,
                        "spacer_len": spacer_len,
                        "left_seq": left,
                        "right_seq": seq[j:j + L],
                        "putative_transposon_len": (j + L) - i  # total span
                    })

                # Look for next occurrence of this right_pattern
                j = seq.find(right_pattern, j + 1)

    return results


# -------------------------------------------------------------------
# Main script
# -------------------------------------------------------------------
if __name__ == "__main__":
    genome_files = [
        r"C:\LaylaYear4\BioInformatics\BioInformatics\Project_L8\L8\corona.fasta",
        r"C:\LaylaYear4\BioInformatics\BioInformatics\Project_L8\L8\Nile.fasta",
        r"C:\LaylaYear4\BioInformatics\BioInformatics\Project_L8\L8\acute.fasta"
    ]

    # You can tweak these if needed
    MIN_LEN = 4
    MAX_LEN = 6
    MAX_SPACER = None  # or e.g. 10000 if you want to limit transposon size

    for filepath in genome_files:
        path = Path(filepath)
        print("=" * 80)
        print(f"Processing genome: {path.name}")
        seq = read_fasta(path)
        print(f"Sequence length: {len(seq)} bp")

        hits = find_inverted_repeats(
            seq,
            min_len=MIN_LEN,
            max_len=MAX_LEN,
            max_spacer=MAX_SPACER
        )
        print(f"Found {len(hits)} inverted-repeat pairs (potential transposons)\n")

        # Show first N hits as example
        N = 20
        for idx, hit in enumerate(hits[:N], start=1):
            print(f"Hit {idx}:")
            print(f"  Left  repeat: {hit['left_seq']}  "
                  f"(pos {hit['left_start_1']}–{hit['left_end_1']})")
            print(f"  Right repeat: {hit['right_seq']} "
                  f"(pos {hit['right_start_1']}–{hit['right_end_1']})")
            print(f"  Repeat length: {hit['repeat_len']} bp")
            print(f"  Spacer length: {hit['spacer_len']} bp")
            print(f"  Putative transposon length (full span): {hit['putative_transposon_len']} bp")
            print()

        if not hits:
            print("  No inverted repeats detected in the specified length range.")

        # ------------------------------------------------------------------
        # PLOTTING
        # ------------------------------------------------------------------
        if hits:
            # Collect spacer lengths and transposon lengths
            spacer_lengths = [h["spacer_len"] for h in hits]
            transposon_lengths = [h["putative_transposon_len"] for h in hits]

            # Histogram of spacer lengths
            plt.figure()
            plt.hist(spacer_lengths, bins=50)
            plt.title(f"Spacer length distribution in {path.name}")
            plt.xlabel("Spacer length (bp)")
            plt.ylabel("Count")
            plt.tight_layout()
            plt.show()

            # Histogram of putative transposon lengths
            plt.figure()
            plt.hist(transposon_lengths, bins=50)
            plt.title(f"Putative transposon length distribution in {path.name}")
            plt.xlabel("Transposon length (bp)")
            plt.ylabel("Count")
            plt.tight_layout()
            plt.show()

            # Distribution of repeat lengths (4, 5, 6)
            repeat_len_counts = {}
            for h in hits:
                L = h["repeat_len"]
                repeat_len_counts[L] = repeat_len_counts.get(L, 0) + 1

            lengths = sorted(repeat_len_counts.keys())
            counts = [repeat_len_counts[L] for L in lengths]

            plt.figure()
            plt.bar(lengths, counts)
            plt.title(f"Repeat length counts in {path.name}")
            plt.xlabel("Repeat length (bp)")
            plt.ylabel("Number of inverted repeat pairs")
            plt.xticks(lengths)
            plt.tight_layout()
            plt.show()
        else:
            print("No hits found, so no plots for this genome.")
