import re
from collections import Counter
from Bio import SeqIO
import matplotlib.pyplot as plt

# Your FASTA file paths
files = [
    r"C:\LaylaYear4\BioInformatics\BioInformatics\Project_L7\L7\corona.fasta",
    r"C:\LaylaYear4\BioInformatics\BioInformatics\Project_L7\L7\Boebeisivirus.fasta",
    r"C:\LaylaYear4\BioInformatics\BioInformatics\Project_L7\L7\acute.fasta",
    r"C:\LaylaYear4\BioInformatics\BioInformatics\Project_L7\L7\HepatitisB.fasta",
    r"C:\LaylaYear4\BioInformatics\BioInformatics\Project_L7\L7\HepatitisC.fasta",
    r"C:\LaylaYear4\BioInformatics\BioInformatics\Project_L7\L7\papillomavirus.fasta",
    r"C:\LaylaYear4\BioInformatics\BioInformatics\Project_L7\L7\Nile.fasta",
    r"C:\LaylaYear4\BioInformatics\BioInformatics\Project_L7\L7\Measles.fasta",
    r"C:\LaylaYear4\BioInformatics\BioInformatics\Project_L7\L7\Zika.fasta",
    r"C:\LaylaYear4\BioInformatics\BioInformatics\Project_L7\L7\Enterobacteria.fasta"
]


def read_sequence(filepath):
    """Reads the first DNA sequence from a FASTA file."""
    record = next(SeqIO.parse(filepath, "fasta"))
    seq = str(record.seq).upper()
    return seq[:min(len(seq), 3000)]


def find_repetitions(seq, min_len=2, max_len=6):
    """Finds patterns of length 6–10 that repeat consecutively."""
    repeats = Counter()
    for l in range(min_len, max_len + 1):
        pattern = re.compile(rf'((?:[ACGT]{{{l}}}))(\1+)')
        for match in pattern.finditer(seq):
            unit = match.group(1)
            repeat_count = len(match.group(0)) // len(unit)
            if repeat_count > 1:
                repeats[unit] += 1
    return repeats


# Prepare figure for all subplots (2 rows × 5 columns = 10 plots)
fig, axes = plt.subplots(2, 5, figsize=(20, 8))
axes = axes.flatten()

for i, filepath in enumerate(files):
    seq = read_sequence(filepath)
    repetitions = find_repetitions(seq)
    title = filepath.split("\\")[-1].replace(".fasta", "")

    ax = axes[i]
    if repetitions:
        patterns = list(repetitions.keys())
        counts = list(repetitions.values())
        ax.bar(range(len(patterns)), counts)
        ax.set_xticks(range(len(patterns)))
        ax.set_xticklabels(patterns, rotation=90, fontsize=6)
        ax.set_ylabel("Repeats", fontsize=8)
        ax.set_xlabel("Pattern", fontsize=8)
    else:
        ax.text(0.5, 0.5, "No repetitions", ha='center', va='center', fontsize=9)

    ax.set_title(title, fontsize=10)

# Adjust layout
plt.tight_layout()
plt.suptitle("Repetition Frequencies in 10 Influenza Genomes", fontsize=14, y=1.02)
plt.show()
