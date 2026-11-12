import os
headers = []
sequences = []
current = ""

with open(os.path.join(os.path.dirname(__file__), "globin_input.fasta"), "r") as f:
    for line in f:
        line = line.strip()

        if not line:
            continue

        if line.startswith(">"):
            if current:
                sequences.append(current)
            headers.append(line)
            current = ""
        else:
            current += line

if current:
    sequences.append(current)

print("Number of records:", len(headers))

for h, s in zip(headers, sequences):
    print(h, "| Length:", len(s))
    from collections import Counter

print("\nAMINO ACID COMPOSITION ANALYSIS")

for h, s in zip(headers, sequences):
    print("\nProtein:", h[:50])

    amino_acids = Counter(s)

    for aa, count in sorted(amino_acids.items()):
        percentage = (count / len(s)) * 100
        print(aa, ":", count, "|", round(percentage, 2), "%")
from collections import Counter

print("\nCOMPARATIVE PROTEIN ANALYSIS")

for h, s in zip(headers, sequences):
    amino_acids = Counter(s)
    aa, count = amino_acids.most_common(1)[0]
    percentage = (count / len(s)) * 100

    print("\nProtein:", h)
    print("Total length:", len(s))
    print("Most abundant amino acid:", aa)
    print("Count:", count)
    print("Percentage:", round(percentage, 2), "%")
    print("\nSEQUENCE LENGTH STATISTICS")

lengths = [len(s) for s in sequences]

print("Total proteins:", len(lengths))
print("Shortest sequence:", min(lengths))
print("Longest sequence:", max(lengths))
print("Average length:", round(sum(lengths)/len(lengths), 2))
print("\nFASTA DATA QUALITY CHECK")

valid_amino_acids = set("ACDEFGHIKLMNPQRSTVWY")

print("Checking protein sequences...")

for i, (h, s) in enumerate(zip(headers, sequences), 1):
    invalid = set(s.upper()) - valid_amino_acids

    print("\nProtein:", i)
    print("Empty sequence:", len(s) == 0)
    print("Invalid amino acids:", invalid if invalid else "None")
    print("Sequence length:", len(s))

print("\nTotal headers:", len(headers))
print("Total sequences:", len(sequences))
print("Header-sequence count matches:", len(headers) == len(sequences))
import os

# Create results folder
base_folder = os.path.dirname(__file__)
results_folder = os.path.join(base_folder, "..", "results")
os.makedirs(results_folder, exist_ok=True)

# Save sequence summary
report_path = os.path.join(results_folder, "globin_analysis_report.txt")

with open(report_path, "w") as report:
    report.write("SWYNEX TASK 2 - GLOBIN ANALYSIS\n")
    report.write("=" * 40 + "\n\n")

    report.write("Total protein sequences: " + str(len(sequences)) + "\n")

    for i, (h, s) in enumerate(zip(headers, sequences), 1):
        report.write("\nProtein " + str(i) + "\n")
        report.write("Header: " + h + "\n")
        report.write("Sequence length: " + str(len(s)) + "\n")

print("\nAnalysis report saved successfully!")
print("Saved at:", report_path)
import matplotlib.pyplot as plt

protein_names = [
    "Horse alpha",
    "Human epsilon",
    "Mouse epsilon-Y2",
    "Chicken beta",
    "Zebrafish beta-1"
]

sequence_lengths = [len(s) for s in sequences]

plt.figure(figsize=(9, 5))
plt.bar(protein_names, sequence_lengths)

plt.title("Globin Protein Sequence Length Comparison")
plt.xlabel("Protein")
plt.ylabel("Sequence Length (Amino Acids)")
plt.xticks(rotation=25)
plt.tight_layout()

# Save graph in results folder
graph_path = os.path.join(results_folder, "sequence_length_graph.png")
plt.savefig(graph_path, dpi=300)
plt.show()

print("Graph saved successfully:", graph_path)

