#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
AI-Designed FASTA Letter Percentage Viewer
------------------------------------------
- Choose a FASTA file (DNA/RNA/Protein) — supports very large files (>100 MB).
- Interprets FASTA: first line is header (info), remaining lines are sequence (often 80 chars/line).
- Streams the file in chunks, so memory stays low and the GUI remains responsive.
- Outputs percentage for ALL letters (A–Z) found in the sequence (case-insensitive).

Run:

"""

import os
import string
import tkinter as tk
from tkinter import filedialog, messagebox, ttk
from collections import Counter

CHUNK_SIZE = 4 * 1024 * 1024  # 4 MB chunks


class FastaLetterApp:
    def __init__(self, root):
        self.root = root
        self.root.title("FASTA Letter Percentage Viewer (AI)")
        self.root.geometry("760x560")
        self.root.minsize(680, 520)

        # Streaming state
        self.fh = None
        self.path = ""
        self.size = 0
        self.bytes_read = 0
        self.header = ""
        self.header_done = False

        self.counter = Counter()
        self.total_letters = 0

        self._build_ui()

    def _build_ui(self):
        wrapper = ttk.Frame(self.root, padding=12)
        wrapper.pack(fill="both", expand=True)

        title = ttk.Label(wrapper, text="FASTA Letter Percentage Viewer", font=("Segoe UI", 16, "bold"))
        title.pack(anchor="w", pady=(0, 8))

        # Buttons
        row = ttk.Frame(wrapper)
        row.pack(fill="x", pady=(0, 8))
        self.btn_open = ttk.Button(row, text="📂 Choose FASTA...", command=self.choose_file)
        self.btn_open.pack(side="left")
        self.btn_reset = ttk.Button(row, text="Reset", command=self.reset, state="disabled")
        self.btn_reset.pack(side="left", padx=(8, 0))

        # Progress
        self.progress = ttk.Progressbar(wrapper, mode="determinate")
        self.progress.pack(fill="x")

        # Info
        info = ttk.LabelFrame(wrapper, text="FASTA Info")
        info.pack(fill="x", pady=(12, 8))
        self.var_path = tk.StringVar(value="File: —")
        self.var_header = tk.StringVar(value="Header: —")
        self.var_len = tk.StringVar(value="Letters counted: —")

        ttk.Label(info, textvariable=self.var_path).pack(anchor="w", padx=8, pady=(6, 2))
        ttk.Label(info, textvariable=self.var_header, wraplength=720, justify="left").pack(anchor="w", padx=8, pady=2)
        ttk.Label(info, textvariable=self.var_len).pack(anchor="w", padx=8, pady=(2, 8))

        # Table
        out = ttk.LabelFrame(wrapper, text="Percentages for Letters (A–Z)")
        out.pack(fill="both", expand=True)

        self.tree = ttk.Treeview(out, columns=("letter", "count", "percent"), show="headings", height=14)
        self.tree.heading("letter", text="Letter")
        self.tree.heading("count", text="Count")
        self.tree.heading("percent", text="Percent")
        self.tree.column("letter", width=100, anchor="center")
        self.tree.column("count", width=140, anchor="e")
        self.tree.column("percent", width=140, anchor="e")
        self.tree.pack(fill="both", expand=True, padx=6, pady=6)

        # Status hint
        ttk.Label(wrapper, text="Tip: Works with very large FASTA files by streaming; only letters A–Z are counted.", foreground="#555").pack(anchor="w", pady=(2, 0))

        # Use a modern theme if available
        style = ttk.Style(self.root)
        if "clam" in style.theme_names():
            style.theme_use("clam")

    def choose_file(self):
        path = filedialog.askopenfilename(
            title="Select FASTA file",
            filetypes=(("FASTA files", "*.fa *.fasta *.fna *.faa *.fas"), ("All files", "*.*")),
        )
        if not path:
            return

        try:
            self.reset()
            self.path = path
            self.size = os.path.getsize(path)
            self.var_path.set(f"File: {path}  ({self.size:,} bytes)")

            # Open file text-mode with universal newlines and error replacement
            self.fh = open(path, "r", encoding="utf-8", errors="replace", newline=None)
            self.progress["value"] = 0
            self.progress["maximum"] = max(1, self.size)
            self.btn_open.config(state="disabled")
            self.btn_reset.config(state="disabled")
            self.root.after(0, self._process_next_chunk)
        except Exception as e:
            messagebox.showerror("Error", f"Could not open file:\n{e}")
            self.reset()

    def reset(self):
        try:
            if self.fh and not self.fh.closed:
                self.fh.close()
        except Exception:
            pass

        self.fh = None
        self.path = ""
        self.size = 0
        self.bytes_read = 0
        self.header = ""
        self.header_done = False

        self.counter.clear()
        self.total_letters = 0

        self.var_path.set("File: —")
        self.var_header.set("Header: —")
        self.var_len.set("Letters counted: —")
        self.progress["value"] = 0
        for iid in self.tree.get_children():
            self.tree.delete(iid)

        self.btn_open.config(state="normal")
        self.btn_reset.config(state="disabled")

    def _process_next_chunk(self):
        if not self.fh:
            return

        try:
            # Header handling (first line only)
            if not self.header_done:
                first = self.fh.readline()
                self.bytes_read = self.fh.tell()
                if not first:
                    self._finish()
                    return
                if first.startswith(">"):
                    self.header = first.strip()
                else:
                    # Some files may lack '>' header — handle gracefully
                    self.header = "(no header line)"
                    self._accumulate_letters(first)
                self.var_header.set(f"Header: {self.header}")
                self.header_done = True
                self.progress["value"] = min(self.bytes_read, self.size)

            # Stream sequence
            chunk = self.fh.read(CHUNK_SIZE)
            if not chunk:
                self._finish()
                return

            self.bytes_read = self.fh.tell()
            self._accumulate_letters(chunk)
            self.progress["value"] = min(self.bytes_read, self.size)

            # Schedule the next chunk so GUI stays responsive
            self.root.after(1, self._process_next_chunk)
        except Exception as e:
            messagebox.showerror("Read Error", str(e))
            self._finish(error=True)

    def _accumulate_letters(self, text):
        # Remove whitespace/newlines, uppercase, keep only letters A–Z
        if not text:
            return
        text = text.replace("\n", "").replace("\r", "").replace(" ", "").upper()
        if not text:
            return
        allowed = set(string.ascii_uppercase)
        # Fast path: build a small Counter from filtered chars
        filtered = [ch for ch in text if ch in allowed]
        if not filtered:
            return
        self.counter.update(filtered)
        self.total_letters += len(filtered)

    def _finish(self, error=False):
        try:
            if self.fh and not self.fh.closed:
                self.fh.close()
        finally:
            self.fh = None

        # Populate table
        for iid in self.tree.get_children():
            self.tree.delete(iid)

        if self.total_letters > 0:
            self.var_len.set(f"Letters counted: {self.total_letters:,}")
            # Show letters found in alphabetical order with percent
            for letter in sorted(self.counter):
                cnt = self.counter[letter]
                pct = (cnt / self.total_letters) * 100.0
                self.tree.insert("", "end", values=(letter, f"{cnt:,}", f"{pct:.2f}%"))
        else:
            self.var_len.set("Letters counted: 0")

        self.btn_reset.config(state="normal")
        if not error:
            messagebox.showinfo("Done", "Finished reading FASTA file and computing letter percentages.")

def main():
    root = tk.Tk()
    FastaLetterApp(root)
    root.mainloop()

if __name__ == "__main__":
    main()
