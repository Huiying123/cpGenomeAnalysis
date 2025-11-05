import sys
from Bio import SeqIO

# --- User Configuration ---
# Your input multi-fasta alignment file
FASTA_FILE = "all.62.alignment.trimmed.fasta"

# Sliding window parameters
WINDOW_SIZE = 600
STEP_SIZE = 200

# Output file for the results table
OUTPUT_FILE = "sliding_window_pi2.tsv"
# -------------------------

def calculate_pi_pairwise(sequences):
    """Calculates nucleotide diversity (pi) by directly comparing all pairs of sequences."""
    n_sequences = len(sequences)
    if n_sequences < 2:
        return 0.0
    
    total_diffs = 0
    total_comparisons = 0
    seq_length = len(sequences[0])

    if seq_length == 0:
        return 0.0

    # Iterate through all possible pairs of sequences (i, j)
    for i in range(n_sequences):
        for j in range(i + 1, n_sequences):
            # Compare the pair site by site
            for k in range(seq_length):
                base1 = sequences[i][k]
                base2 = sequences[j][k]
                
                # Only compare valid bases, ignore gaps and Ns
                if base1 in ('A', 'T', 'G', 'C') and base2 in ('A', 'T', 'G', 'C'):
                    total_comparisons += 1
                    if base1 != base2:
                        total_diffs += 1
    
    if total_comparisons == 0:
        return 0.0
    
    pi = total_diffs / total_comparisons
    return pi

# --- Main Script ---
print(f"Reading alignment from: {FASTA_FILE}")
try:
    # Read all records into memory as strings for speed
    records = [str(record.seq).upper() for record in SeqIO.parse(FASTA_FILE, "fasta")]
    alignment_length = len(records[0])
    print(f"Found {len(records)} sequences of length {alignment_length} bp.")
except FileNotFoundError:
    print(f"ERROR: File not found at '{FASTA_FILE}'")
    sys.exit(1)

print(f"Calculating Pi in windows (size={WINDOW_SIZE}, step={STEP_SIZE})...")

with open(OUTPUT_FILE, "w") as f_out:
    f_out.write("Window_Midpoint\tPi\n")

    for i in range(0, alignment_length - WINDOW_SIZE + 1, STEP_SIZE):
        start = i
        end = i + WINDOW_SIZE
        midpoint = start + (WINDOW_SIZE / 2)

        # Extract the sequences for the current window
        window_sequences = [rec[start:end] for rec in records]
        
        pi_for_window = calculate_pi_pairwise(window_sequences)
        
        f_out.write(f"{midpoint}\t{pi_for_window:.6f}\n")
        
        if i % (STEP_SIZE * 10) == 0:
             print(f"Window {start}-{end}... Pi = {pi_for_window:.6f}")

print(f"\nAnalysis complete! Results saved to {OUTPUT_FILE}")
