import pandas as pd
import matplotlib.pyplot as plt

# Path to your .narrowPeak file
file = input("Enter the path to your .narrowPeak file: ").strip()
# file = "your_peaks.narrowPeak.bed"

# Read the tab-separated file (no header, BED-like)
df = pd.read_csv(file, sep="\t", header=None)

# BED: columns 1=start, 2=end
df["length"] = df[2] - df[1]


# --- Define bins (last bin is >5000 bp) ---
bins = [0, 50, 100,300, 500,700, 1000, 2000, 3000, 4000, 5000, float("inf")]
labels = ["0-50","50-100", "100-300","300-500", "500-700","700-1000", "1000-2000", "2000-3000", "3000-4000", "4000-5000", ">5000"]

# --- Compute frequencies ---
df["bin"] = pd.cut(df["length"], bins=bins, labels=labels, right=False)
freq = df["bin"].value_counts().sort_index()

# --- Print frequency table ---
print("\nPeak Length Distribution:")
print(freq)

# --- Plot histogram ---
plt.figure(figsize=(8, 5))
plt.bar(freq.index.astype(str), freq.values, color="skyblue", edgecolor="black")
plt.title(f"Peak Length Distribution ({file})")
plt.xlabel("Peak length (bp)")
plt.ylabel("Frequency")
plt.xticks(rotation=45)
plt.tight_layout()
plt.show()