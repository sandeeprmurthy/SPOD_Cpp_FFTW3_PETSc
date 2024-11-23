import os
import re
import matplotlib.pyplot as plt
import numpy as np

# Enable LaTeX formatting
plt.rcParams.update({
    "text.usetex": True,
    "font.family": "serif",
    "axes.labelsize": 12,
    "font.size": 10,
    "legend.fontsize": 10
})

# Directory containing eigenvalue files
eigenvalues_dir = "eigenvalues"

# Regular expression to extract frequency index
freq_pattern = re.compile(r"eigenvalues_freq_(\d+)\.txt")

# Dictionary to store eigenvalues from all files
eigenvalues_by_order = {}

# Read all eigenvalue files
for filename in sorted(os.listdir(eigenvalues_dir), key=lambda x: int(freq_pattern.search(x).group(1))):
    if filename.startswith("eigenvalues_freq_") and filename.endswith(".txt"):
        freq = int(freq_pattern.search(filename).group(1))
        filepath = os.path.join(eigenvalues_dir, filename)
        
        with open(filepath, "r") as file:
            lines = file.readlines()
            for idx, line in enumerate(lines):
                eigenvalue = float(line.strip()[1:-1].split(",")[0])  # Extract real part of eigenvalue
                if idx not in eigenvalues_by_order:
                    eigenvalues_by_order[idx] = []
                eigenvalues_by_order[idx].append((freq, eigenvalue))

# Plotting
plt.figure(figsize=(10, 6))
for order, values in eigenvalues_by_order.items():
    freqs, eigenvalues = zip(*values)
    plt.plot(freqs, eigenvalues, label=f"Eigenvalue {order + 1}")

# Configure plot
plt.xscale("log")
plt.yscale("log")
plt.xlabel(r"Frequency (log scale)")
plt.ylabel(r"Eigenvalue Magnitude (log scale)")
plt.title(r"Eigenvalues Across Frequencies")
plt.legend()
plt.grid(True, which="both", linestyle="--", linewidth=0.5)

# Save and display plot
plt.savefig("eigenvalues_plot.pdf", bbox_inches="tight")
plt.show()

