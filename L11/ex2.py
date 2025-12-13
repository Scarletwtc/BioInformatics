import os
import re
import itertools
from typing import Dict, Tuple, List

import numpy as np
import matplotlib.pyplot as plt


# =========================
# SETTINGS (change if needed)
# =========================
MATCH_SCORE = 1
MISMATCH_SCORE = -1
GAP_SCORE = -1

# Chunked NW strategy for long genomes
CHUNK_SIZE = 350
OVERLAP = 70
SEARCH_BACK = 150
SEARCH_FORWARD = 250

# Console alignment printing
WRAP_WIDTH = 100
PRINT_HEAD_BLOCKS = 3     # show first N blocks
PRINT_TAIL_BLOCKS = 1     # show last N blocks (0 to disable)

# Fragmented "teacher-style" blocks (alignment-axis)
WIN_COLS = 250            # window size in alignment columns
STEP_COLS = 50            # sliding step in columns
IDENTITY_THRESH = 90    # increase (96-98) for MORE fragmentation
MIN_COMPARED = 80         # min non-gap comparisons inside window


# =========================
# Loading sequences (FASTA / raw / GenBank)
# =========================
def load_sequence(path: str) -> str:
    """
    Loads sequence from FASTA, raw text, or GenBank.
    Returns uppercase string containing only A/C/G/T/N.
    """
    with open(path, "r", encoding="utf-8", errors="ignore") as f:
        text = f.read()

    # GenBank: ORIGIN ... //
    if "ORIGIN" in text and "//" in text:
        origin = text.split("ORIGIN", 1)[1].split("//", 1)[0]
        seq = re.sub(r"[^ACGTNacgtn]", "", origin).upper()
        return seq

    lines = [ln.strip() for ln in text.splitlines() if ln.strip()]
    if not lines:
        raise ValueError(f"{path} is empty.")

    # FASTA
    if lines[0].startswith(">"):
        seq = "".join(ln for ln in lines if not ln.startswith(">"))
    else:
        seq = "".join(lines)

    seq = re.sub(r"[^ACGTNacgtn]", "", seq).upper()
    return seq


def discover_sequence_files(folder: str) -> List[str]:
    """
    Finds sequence1..sequence10 files (any extension).
    e.g. sequence1, sequence1.fasta, sequence1.txt, sequence1.gb
    """
    files = []
    for name in os.listdir(folder):
        if re.match(r"^sequence(10|[1-9])(\..+)?$", name, flags=re.IGNORECASE):
            files.append(os.path.join(folder, name))

    def idx(p: str) -> int:
        base = os.path.basename(p)
        m = re.match(r"^sequence(10|[1-9])", base, flags=re.IGNORECASE)
        return int(m.group(1)) if m else 999

    files.sort(key=idx)
    return files


# =========================
# Needleman–Wunsch (SMALL sequences)
# =========================
def nw_global_align(a: str, b: str, match: int, mismatch: int, gap: int) -> Tuple[str, str]:
    n, m = len(a), len(b)
    score = np.zeros((n + 1, m + 1), dtype=int)
    tb = np.zeros((n + 1, m + 1), dtype=np.int8)  # 0 diag, 1 up, 2 left

    for i in range(1, n + 1):
        score[i, 0] = score[i - 1, 0] + gap
        tb[i, 0] = 1
    for j in range(1, m + 1):
        score[0, j] = score[0, j - 1] + gap
        tb[0, j] = 2

    for i in range(1, n + 1):
        ai = a[i - 1]
        for j in range(1, m + 1):
            diag = score[i - 1, j - 1] + (match if ai == b[j - 1] else mismatch)
            up = score[i - 1, j] + gap
            left = score[i, j - 1] + gap

            best = diag
            direction = 0
            if up > best:
                best = up
                direction = 1
            if left > best:
                best = left
                direction = 2

            score[i, j] = best
            tb[i, j] = direction

    # Traceback
    i, j = n, m
    out_a, out_b = [], []
    while i > 0 or j > 0:
        if i > 0 and j > 0 and tb[i, j] == 0:
            out_a.append(a[i - 1])
            out_b.append(b[j - 1])
            i -= 1
            j -= 1
        elif i > 0 and (j == 0 or tb[i, j] == 1):
            out_a.append(a[i - 1])
            out_b.append("-")
            i -= 1
        else:
            out_a.append("-")
            out_b.append(b[j - 1])
            j -= 1

    return "".join(reversed(out_a)), "".join(reversed(out_b))


# =========================
# Chunked strategy (BIG genomes)
# =========================
def trim_overlap(aln_a: str, aln_b: str, overlap_letters: int) -> Tuple[str, str, int]:
    """
    Trim columns from LEFT until overlap_letters non-gap letters removed from aln_a.
    Returns (trimmed_a, trimmed_b, removed_b_letters)
    """
    if overlap_letters <= 0:
        return aln_a, aln_b, 0

    removed_a = 0
    removed_b = 0
    cut = 0

    while cut < len(aln_a) and removed_a < overlap_letters:
        if aln_a[cut] != "-":
            removed_a += 1
        if aln_b[cut] != "-":
            removed_b += 1
        cut += 1

    return aln_a[cut:], aln_b[cut:], removed_b


def identity_percent(aln_a: str, aln_b: str) -> float:
    matches = 0
    compared = 0
    for x, y in zip(aln_a, aln_b):
        if x != "-" and y != "-":
            compared += 1
            if x == y:
                matches += 1
    return (matches / compared * 100.0) if compared else 0.0


def chunked_pairwise_alignment(seq1: str, seq2: str) -> Tuple[str, str, float]:
    """
    Step-by-step alignment:
    chunk seq1, align chunk vs window in seq2, stitch with overlap trimming.
    Returns aligned_seq1, aligned_seq2, identity%.
    """
    p1 = 0
    p2 = 0
    first = True

    out1, out2 = [], []

    while p1 < len(seq1):
        chunk1 = seq1[p1: p1 + CHUNK_SIZE]
        if not chunk1:
            break

        w_start = max(0, p2 - SEARCH_BACK)
        w_end = min(len(seq2), p2 + len(chunk1) + SEARCH_FORWARD)
        window2 = seq2[w_start:w_end]

        aln1, aln2 = nw_global_align(chunk1, window2, MATCH_SCORE, MISMATCH_SCORE, GAP_SCORE)

        removed_b = 0
        if not first and OVERLAP > 0:
            aln1, aln2, removed_b = trim_overlap(aln1, aln2, OVERLAP)

        out1.append(aln1)
        out2.append(aln2)

        # Advance pointers
        step1 = CHUNK_SIZE if first else (CHUNK_SIZE - OVERLAP)
        p1 += step1

        kept_b = sum(1 for c in aln2 if c != "-")
        p2 = w_start + removed_b + kept_b

        first = False
        if step1 <= 0:
            break

    aligned1 = "".join(out1)
    aligned2 = "".join(out2)
    sim = identity_percent(aligned1, aligned2)
    return aligned1, aligned2, sim


# =========================
# Console printing
# =========================
def match_line(a: str, b: str) -> str:
    return "".join("|" if (x == y and x != "-" and y != "-") else " " for x, y in zip(a, b))


def print_alignment_console(name_a: str, name_b: str, aln_a: str, aln_b: str, similarity: float):
    print("\n" + "=" * 100)
    print(f"ALIGNMENT: {name_a} vs {name_b}   | Identity: {similarity:.2f}%")
    print("=" * 100)

    mline = match_line(aln_a, aln_b)
    width = WRAP_WIDTH
    total_blocks = (len(aln_a) + width - 1) // width

    head = min(PRINT_HEAD_BLOCKS, total_blocks)
    tail = min(PRINT_TAIL_BLOCKS, max(0, total_blocks - head))

    # Print head blocks
    for blk in range(head):
        s = blk * width
        e = min(len(aln_a), (blk + 1) * width)
        print(aln_a[s:e])
        print(mline[s:e])
        print(aln_b[s:e])
        print()

    # Print tail blocks
    if tail > 0:
        if total_blocks > head + tail:
            print("... (middle omitted) ...\n")
        start_blk = total_blocks - tail
        for blk in range(start_blk, total_blocks):
            s = blk * width
            e = min(len(aln_a), (blk + 1) * width)
            print(aln_a[s:e])
            print(mline[s:e])
            print(aln_b[s:e])
            print()


# =========================
# Teacher-style fragmented blocks (alignment-axis)
# =========================
def alignment_blocks(aln_a: str, aln_b: str,
                     win_cols: int, step_cols: int,
                     id_thresh: float, min_compared: int) -> List[Tuple[int, int]]:
    """
    Compute high-identity blocks along the ALIGNMENT axis (columns).
    Returns list of (start_col, end_col).
    """
    L = len(aln_a)
    good = np.zeros(L, dtype=bool)

    for s in range(0, max(1, L - win_cols + 1), step_cols):
        e = min(L, s + win_cols)

        matches = 0
        compared = 0
        for x, y in zip(aln_a[s:e], aln_b[s:e]):
            if x != "-" and y != "-":
                compared += 1
                if x == y:
                    matches += 1

        if compared >= min_compared:
            ident = matches / compared * 100.0
            if ident >= id_thresh:
                good[s:e] = True

    segs = []
    i = 0
    while i < L:
        if not good[i]:
            i += 1
            continue
        j = i
        while j < L and good[j]:
            j += 1
        segs.append((i, j))
        i = j

    return segs


def plot_blocks_like_teacher(name_a: str, name_b: str,
                            aln_a: str, aln_b: str,
                            out_path: str):
    """
    Draw two tracks (G1, G2). Similar regions are rectangular bars at the SAME x-positions
    (alignment-axis), producing the "long bar - gap - short bar ..." look.
    """
    segs = alignment_blocks(
        aln_a, aln_b,
        win_cols=WIN_COLS,
        step_cols=STEP_COLS,
        id_thresh=IDENTITY_THRESH,
        min_compared=MIN_COMPARED
    )

    fig, ax = plt.subplots(figsize=(12, 2.4))

    y1, y2, h = 12, 4, 5

    # Faint baseline tracks (full alignment length)
    ax.broken_barh([(0, len(aln_a))], (y1, h), alpha=0.12)
    ax.broken_barh([(0, len(aln_a))], (y2, h), alpha=0.12)

    # Bars at same x positions for both genomes
    if segs:
        blocks = [(s, e - s) for s, e in segs]
        ax.broken_barh(blocks, (y1, h), alpha=0.9)
        ax.broken_barh(blocks, (y2, h), alpha=0.9)

    ax.set_yticks([y1 + h / 2, y2 + h / 2])
    ax.set_yticklabels([name_a, name_b])
    ax.set_xlabel("Alignment position (columns)")
    ax.set_title(f"Similar blocks (≥ {IDENTITY_THRESH:.0f}% windows) — {name_a} vs {name_b}")
    ax.set_xlim(0, len(aln_a))
    ax.set_ylim(0, 22)
    ax.grid(True, axis="x", alpha=0.25)

    plt.tight_layout()
    fig.savefig(out_path, dpi=150)
    plt.close(fig)


# =========================
# Main (all pairwise)
# =========================
def main():
    folder = "."
    out_dir = "results"
    plots_dir = os.path.join(out_dir, "plots")
    os.makedirs(out_dir, exist_ok=True)
    os.makedirs(plots_dir, exist_ok=True)

    files = discover_sequence_files(folder)
    if len(files) < 10:
        raise RuntimeError(f"Found {len(files)} files. Need sequence1..sequence10 in {os.path.abspath(folder)}")

    # Load sequences
    seqs: Dict[str, str] = {}
    for f in files[:10]:
        name = os.path.splitext(os.path.basename(f))[0]  # "sequence1"
        seqs[name] = load_sequence(f)
        print(f"Loaded {name}: {len(seqs[name])} bp from {os.path.basename(f)}")

    names = sorted(seqs.keys(), key=lambda x: int(re.findall(r"\d+", x)[0]))

    # Similarity matrix
    sim = {a: {b: (100.0 if a == b else None) for b in names} for a in names}

    pairs = list(itertools.combinations(names, 2))
    pairwise_rows = []   # will hold one row per alignment

    print(f"\nRunning {len(pairs)} pairwise alignments...\n")

    for idx, (a, b) in enumerate(pairs, 1):
        print(f"[{idx}/{len(pairs)}] Aligning {a} vs {b} ...")

        aln_a, aln_b, similarity = chunked_pairwise_alignment(seqs[a], seqs[b])

        sim[a][b] = similarity
        sim[b][a] = similarity

        # 1) Show alignment in console (head+tail)
        print_alignment_console(a, b, aln_a, aln_b, similarity)

        pairwise_rows.append({
            "Genome_A": a,
            "Genome_B": b,
            "Len_A": len(seqs[a]),
            "Len_B": len(seqs[b]),
            "Identity_%": round(similarity, 2)
        })


        # Save alignment text
        txt_path = os.path.join(out_dir, f"alignment_{a}_vs_{b}.txt")
        with open(txt_path, "w", encoding="utf-8") as f:
            f.write(f"{a} length: {len(seqs[a])} bp\n")
            f.write(f"{b} length: {len(seqs[b])} bp\n")
            f.write(f"Identity: {similarity:.2f} %\n\n")
            f.write(aln_a + "\n")
            f.write(match_line(aln_a, aln_b) + "\n")
            f.write(aln_b + "\n")

        # 2) Teacher-style fragmented bars (alignment-axis)
        plot_path = os.path.join(plots_dir, f"blocks_{a}_vs_{b}.png")
        plot_blocks_like_teacher(a, b, aln_a, aln_b, plot_path)

        print(f"    -> Identity: {similarity:.0f}% | saved: {txt_path}")
        print(f"    -> Blocks plot: {plot_path}\n")

    # Save similarity matrix CSV
    csv_path = os.path.join(out_dir, "similarity_matrix.csv")
    with open(csv_path, "w", encoding="utf-8") as f:
        f.write("," + ",".join(names) + "\n")
        for r in names:
            row = [r] + [f"{sim[r][c]:.2f}" for c in names]
            f.write(",".join(row) + "\n")
    
    # Save pairwise similarity table
    table_path = os.path.join(out_dir, "pairwise_similarity_table.csv")
    with open(table_path, "w", encoding="utf-8") as f:
        f.write("Genome_A,Genome_B,Len_A,Len_B,Identity_%\n")
        for row in pairwise_rows:
            f.write(f"{row['Genome_A']},{row['Genome_B']},{row['Len_A']},{row['Len_B']},{row['Identity_%']}\n")

    print(f"\nPairwise similarity table saved to: {table_path}")


    print("Done.")
    print(f"- Alignments: {os.path.abspath(out_dir)}")
    print(f"- Plots:      {os.path.abspath(plots_dir)}")
    print(f"- Matrix:     {csv_path}")
    print(f"- Pairwise table: {table_path}")


    print("\nPAIRWISE SIMILARITY SCORES")
    print("-" * 60)
    for row in pairwise_rows:
        print(f"{row['Genome_A']:>10} vs {row['Genome_B']:<10}  = {row['Identity_%']:6.2f}%")



if __name__ == "__main__":
    main()
