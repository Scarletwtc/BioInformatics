import tkinter as tk
from tkinter import ttk
import numpy as np
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from matplotlib.figure import Figure
from matplotlib.colors import ListedColormap


class NeedlemanWunschAligner:
    """Needleman–Wunsch global alignment (score matrix + traceback)."""

    # Traceback codes
    DIAG = 0
    UP = 1
    LEFT = 2

    def __init__(self, seq1, seq2, match_score=1, mismatch_score=-1, gap_score=0):
        self.seq1 = seq1.upper().strip()
        self.seq2 = seq2.upper().strip()
        self.match_score = match_score
        self.mismatch_score = mismatch_score
        self.gap_score = gap_score

        self.score_matrix = None
        self.traceback_matrix = None

    def align(self):
        n = len(self.seq1)
        m = len(self.seq2)

        self.score_matrix = np.zeros((n + 1, m + 1), dtype=float)
        self.traceback_matrix = np.zeros((n + 1, m + 1), dtype=int)

        # Initialize first row/col
        for i in range(1, n + 1):
            self.score_matrix[i, 0] = i * self.gap_score
            self.traceback_matrix[i, 0] = self.UP

        for j in range(1, m + 1):
            self.score_matrix[0, j] = j * self.gap_score
            self.traceback_matrix[0, j] = self.LEFT

        # Fill matrices
        for i in range(1, n + 1):
            for j in range(1, m + 1):
                diag = self.score_matrix[i - 1, j - 1] + (
                    self.match_score if self.seq1[i - 1] == self.seq2[j - 1] else self.mismatch_score
                )
                up = self.score_matrix[i - 1, j] + self.gap_score
                left = self.score_matrix[i, j - 1] + self.gap_score

                best = max(diag, up, left)
                self.score_matrix[i, j] = best

                # Tie-breaking: prefer DIAG, then UP, then LEFT (common/simple choice)
                if best == diag:
                    self.traceback_matrix[i, j] = self.DIAG
                elif best == up:
                    self.traceback_matrix[i, j] = self.UP
                else:
                    self.traceback_matrix[i, j] = self.LEFT

        return self.traceback()

    def traceback(self):
        aligned1 = []
        aligned2 = []
        match_line = []

        i, j = len(self.seq1), len(self.seq2)
        path = [(i, j)]

        while i > 0 or j > 0:
            if i > 0 and j > 0 and self.traceback_matrix[i, j] == self.DIAG:
                a = self.seq1[i - 1]
                b = self.seq2[j - 1]
                aligned1.append(a)
                aligned2.append(b)
                match_line.append("|" if a == b else " ")
                i -= 1
                j -= 1

            elif i > 0 and (j == 0 or self.traceback_matrix[i, j] == self.UP):
                a = self.seq1[i - 1]
                aligned1.append(a)
                aligned2.append("-")
                match_line.append(" ")
                i -= 1

            else:
                b = self.seq2[j - 1]
                aligned1.append("-")
                aligned2.append(b)
                match_line.append(" ")
                j -= 1

            path.append((i, j))

        aligned1 = "".join(reversed(aligned1))
        aligned2 = "".join(reversed(aligned2))
        match_line = "".join(reversed(match_line))
        path = list(reversed(path))

        matches = sum(1 for a, b in zip(aligned1, aligned2) if a == b and a != "-")
        length = len(aligned1)
        similarity = (matches / length * 100) if length else 0.0

        return {
            "seq1": aligned1,
            "seq2": aligned2,
            "match_line": match_line,
            "matches": matches,
            "length": length,
            "similarity": similarity,
            "traceback_path": path,
        }


class NeedlemanWunschGUI:
    """GUI for Needleman–Wunsch alignment with matrix + traceback deviation views."""

    def __init__(self, root):
        self.root = root
        self.root.title("Needleman-Wunsch Sequence Aligner")
        self.root.geometry("1400x820")

        self.seq1_default = "ACCGTGAAGCCAATAC"
        self.seq2_default = "AGCGTGCAGCCAATAC"

        self.plot_traceback = tk.BooleanVar(value=True)
        self.plot_grid = tk.BooleanVar(value=True)

        self.create_widgets()

    def create_widgets(self):
        # Main layout: left controls, right visuals+output
        left_frame = ttk.Frame(self.root, padding="10")
        left_frame.grid(row=0, column=0, sticky="nsw")

        right_frame = ttk.Frame(self.root, padding="10")
        right_frame.grid(row=0, column=1, sticky="nsew")

        self.root.columnconfigure(1, weight=1)
        self.root.rowconfigure(0, weight=1)
        right_frame.columnconfigure(0, weight=1)
        right_frame.rowconfigure(0, weight=1)

        # === Left: sequences ===
        seq_frame = ttk.LabelFrame(left_frame, text="Sequences", padding="10")
        seq_frame.grid(row=0, column=0, sticky="ew", pady=6)

        ttk.Label(seq_frame, text="Sq 1 =").grid(row=0, column=0, sticky="w", pady=2)
        self.seq1_entry = ttk.Entry(seq_frame, width=28)
        self.seq1_entry.grid(row=0, column=1, pady=2)
        self.seq1_entry.insert(0, self.seq1_default)

        ttk.Label(seq_frame, text="Sq 2 =").grid(row=1, column=0, sticky="w", pady=2)
        self.seq2_entry = ttk.Entry(seq_frame, width=28)
        self.seq2_entry.grid(row=1, column=1, pady=2)
        self.seq2_entry.insert(0, self.seq2_default)

        # === Left: parameters ===
        param_frame = ttk.LabelFrame(left_frame, text="Parameters", padding="10")
        param_frame.grid(row=1, column=0, sticky="ew", pady=6)

        ttk.Label(param_frame, text="Gap =").grid(row=0, column=0, sticky="w")
        self.gap_entry = ttk.Entry(param_frame, width=10)
        self.gap_entry.grid(row=0, column=1, padx=5, pady=2)
        self.gap_entry.insert(0, "0")

        ttk.Label(param_frame, text="Mach =").grid(row=1, column=0, sticky="w")
        self.match_entry = ttk.Entry(param_frame, width=10)
        self.match_entry.grid(row=1, column=1, padx=5, pady=2)
        self.match_entry.insert(0, "1")

        ttk.Label(param_frame, text="MMach =").grid(row=2, column=0, sticky="w")
        self.mismatch_entry = ttk.Entry(param_frame, width=10)
        self.mismatch_entry.grid(row=2, column=1, padx=5, pady=2)
        self.mismatch_entry.insert(0, "-1")

        # === Left: options ===
        options_frame = ttk.LabelFrame(left_frame, text="Options", padding="10")
        options_frame.grid(row=2, column=0, sticky="ew", pady=6)

        ttk.Checkbutton(options_frame, text="Plot TraceBack", variable=self.plot_traceback)\
            .grid(row=0, column=0, sticky="w")
        ttk.Checkbutton(options_frame, text="Plot grid", variable=self.plot_grid)\
            .grid(row=1, column=0, sticky="w")

        ttk.Button(left_frame, text="Align", command=self.perform_alignment)\
            .grid(row=3, column=0, pady=14, ipadx=24, ipady=6)

        # === Right: plots side-by-side ===
        plots_outer = ttk.Frame(right_frame)
        plots_outer.grid(row=0, column=0, sticky="nsew")
        plots_outer.columnconfigure(0, weight=1)
        plots_outer.columnconfigure(1, weight=1)
        plots_outer.rowconfigure(0, weight=1)

        self.matrix_frame = ttk.Frame(plots_outer)
        self.matrix_frame.grid(row=0, column=0, sticky="nsew", padx=(0, 10))

        self.traceback_frame = ttk.Frame(plots_outer)
        self.traceback_frame.grid(row=0, column=1, sticky="nsew")

        # === Right: results text ===
        result_frame = ttk.LabelFrame(right_frame, text="Show Alignment:", padding="10")
        result_frame.grid(row=1, column=0, sticky="ew", pady=(10, 0))

        self.result_text = tk.Text(result_frame, height=16, width=90, font=("Courier", 10))
        self.result_text.grid(row=0, column=0, sticky="ew")

        scrollbar = ttk.Scrollbar(result_frame, orient="vertical", command=self.result_text.yview)
        scrollbar.grid(row=0, column=1, sticky="ns")
        self.result_text.configure(yscrollcommand=scrollbar.set)

        result_frame.columnconfigure(0, weight=1)

        # Do one initial alignment so it looks “ready”
        self.perform_alignment()

    def perform_alignment(self):
        seq1 = self.seq1_entry.get().strip().upper()
        seq2 = self.seq2_entry.get().strip().upper()

        try:
            gap_score = int(self.gap_entry.get())
            match_score = int(self.match_entry.get())
            mismatch_score = int(self.mismatch_entry.get())
        except ValueError:
            self.result_text.delete("1.0", tk.END)
            self.result_text.insert("1.0", "Error: Invalid parameter values!")
            return

        if not seq1 or not seq2:
            self.result_text.delete("1.0", tk.END)
            self.result_text.insert("1.0", "Error: Please enter both sequences!")
            return

        aligner = NeedlemanWunschAligner(seq1, seq2, match_score, mismatch_score, gap_score)
        result = aligner.align()

        self.display_results(result, len(seq1), len(seq2))
        self.visualize_matrix(aligner, result)
        self.visualize_traceback_deviation(aligner, result)

    def display_results(self, result, n, m):
        self.result_text.delete("1.0", tk.END)

        # Teacher screenshot shows 77% (truncation), not 78% (rounding).
        similarity_int = int(result["similarity"])

        output = ""
        output += f"{result['seq1']}\n"
        output += f"{result['match_line']}\n"
        output += f"{result['seq2']}\n\n"
        output += f"Matches = {result['matches']}\n"
        output += f"Lenght = {result['length']}\n"
        output += f"\nSimilarity = {similarity_int} %\n\n"
        output += f"Tracing back: M[{n},{m}]\n"

        self.result_text.insert("1.0", output)

    def visualize_matrix(self, aligner, result):
        for w in self.matrix_frame.winfo_children():
            w.destroy()

        fig = Figure(figsize=(6.6, 5.2))
        ax = fig.add_subplot(111)

        # Similar look to your teacher screenshot
        im = ax.imshow(
            aligner.score_matrix,
            cmap="magma",
            origin="upper",
            aspect="equal"
        )

        ax.set_title("Graphic representation of the alignment matrix", fontsize=11)

        # Optional grid overlay
        if self.plot_grid.get():
            ax.set_xticks(np.arange(-0.5, aligner.score_matrix.shape[1], 1), minor=True)
            ax.set_yticks(np.arange(-0.5, aligner.score_matrix.shape[0], 1), minor=True)
            ax.grid(which="minor", linewidth=0.5, alpha=0.25)
            ax.tick_params(which="minor", bottom=False, left=False)

        # Traceback path line (optional)
        if self.plot_traceback.get() and result.get("traceback_path"):
            path = result["traceback_path"]
            ys = [p[0] for p in path]
            xs = [p[1] for p in path]
            ax.plot(xs, ys, linewidth=2, alpha=0.9)

        # Tick labels like the sample (letters along borders)
        ax.set_xticks(range(len(aligner.seq2) + 1))
        ax.set_yticks(range(len(aligner.seq1) + 1))
        ax.set_xticklabels(["-"] + list(aligner.seq2), fontsize=8)
        ax.set_yticklabels(["-"] + list(aligner.seq1), fontsize=8)

        fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
        fig.tight_layout()

        canvas = FigureCanvasTkAgg(fig, master=self.matrix_frame)
        canvas.draw()
        canvas.get_tk_widget().pack(fill=tk.BOTH, expand=True)

    def visualize_traceback_deviation(self, aligner, result):
        for w in self.traceback_frame.winfo_children():
            w.destroy()

        n = len(aligner.seq1)
        m = len(aligner.seq2)

        fig = Figure(figsize=(6.6, 5.2))
        ax = fig.add_subplot(111)

        deviation = np.zeros((n + 1, m + 1), dtype=int)
        if result.get("traceback_path"):
            for i, j in result["traceback_path"]:
                deviation[i, j] = 1

        # Yellow background + red path (like your teacher screenshot)
        cmap = ListedColormap(["#FFF7B3", "#C8382B"])
        ax.imshow(deviation, cmap=cmap, origin="upper", aspect="equal")

        ax.set_title("Traceback path deviation from optimal alignment (diagonal)", fontsize=11)

        # Strong grid like the screenshot
        ax.set_xticks(np.arange(-0.5, m + 1, 1), minor=True)
        ax.set_yticks(np.arange(-0.5, n + 1, 1), minor=True)
        ax.grid(which="minor", linewidth=1.0)
        ax.tick_params(which="minor", bottom=False, left=False)

        # Labels along the edges
        ax.set_xticks(range(m + 1))
        ax.set_yticks(range(n + 1))
        ax.set_xticklabels(["-"] + list(aligner.seq2), fontsize=8)
        ax.set_yticklabels(["-"] + list(aligner.seq1), fontsize=8)

        fig.tight_layout()

        canvas = FigureCanvasTkAgg(fig, master=self.traceback_frame)
        canvas.draw()
        canvas.get_tk_widget().pack(fill=tk.BOTH, expand=True)


def main():
    root = tk.Tk()
    app = NeedlemanWunschGUI(root)
    root.mainloop()


if __name__ == "__main__":
    main()
