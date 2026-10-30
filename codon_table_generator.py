import pandas as pd

codon_table = {
    "TTT": "Phenylalanine", "TTC": "Phenylalanine", "TTA": "Leucine", "TTG": "Leucine",
    "TCT": "Serine", "TCC": "Serine", "TCA": "Serine", "TCG": "Serine",
    "TAT": "Tyrosine", "TAC": "Tyrosine", "TAA": "Stop", "TAG": "Stop",
    "TGT": "Cysteine", "TGC": "Cysteine", "TGA": "Stop", "TGG": "Tryptophan",
    "CTT": "Leucine", "CTC": "Leucine", "CTA": "Leucine", "CTG": "Leucine",
    "CCT": "Proline", "CCC": "Proline", "CCA": "Proline", "CCG": "Proline",
    "CAT": "Histidine", "CAC": "Histidine", "CAA": "Glutamine", "CAG": "Glutamine",
    "CGT": "Arginine", "CGC": "Arginine", "CGA": "Arginine", "CGG": "Arginine",
    "ATT": "Isoleucine", "ATC": "Isoleucine", "ATA": "Isoleucine", "ATG": "Methionine",
    "ACT": "Threonine", "ACC": "Threonine", "ACA": "Threonine", "ACG": "Threonine",
    "AAT": "Asparagine", "AAC": "Asparagine", "AAA": "Lysine", "AAG": "Lysine",
    "AGT": "Serine", "AGC": "Serine", "AGA": "Arginine", "AGG": "Arginine",
    "GTT": "Valine", "GTC": "Valine", "GTA": "Valine", "GTG": "Valine",
    "GCT": "Alanine", "GCC": "Alanine", "GCA": "Alanine", "GCG": "Alanine",
    "GAT": "Aspartic acid", "GAC": "Aspartic acid", "GAA": "Glutamic acid", "GAG": "Glutamic acid",
    "GGT": "Glycine", "GGC": "Glycine", "GGA": "Glycine", "GGG": "Glycine",
}

start_codons = ["ATG"]
stop_codons = ["TAA", "TAG", "TGA"]

print("Total codons:", len(codon_table))


bases = ["T", "C", "A", "G"]
grid = []
for first in bases:
    for third in bases:
        row = {}
        for second in bases:
            codon = first + second + third
            row[second] = codon + " - " + codon_table[codon]
        grid.append(row)

grid_df = pd.DataFrame(grid)
grid_df.columns = ["2nd base T", "2nd base C", "2nd base A", "2nd base G"]
grid_df.index = [f + "_" + t for f in bases for t in bases]
grid_df.index.name = "1st _ 3rd"
print(grid_df)


def clean_sequence(seq):
    seq = seq.upper().replace(" ", "").replace("\n", "")
    return seq


def is_valid_dna(seq):
    for base in seq:
        if base not in "ATGC":
            return False
    return True


def split_codons(seq):
    codons = []
    for i in range(0, len(seq) - len(seq) % 3, 3):
        codons.append(seq[i:i + 3])
    leftover = seq[len(seq) - len(seq) % 3:]
    return codons, leftover


def translate(seq):
    seq = clean_sequence(seq)
    if seq == "":
        print("Error: the sequence is empty.")
        return None
    if not is_valid_dna(seq):
        bad = []
        for base in seq:
            if base not in "ATGC" and base not in bad:
                bad.append(base)
        print("Error: invalid characters found:", ", ".join(bad))
        print("A DNA sequence can only have A, T, G and C.")
        return None

    codons, leftover = split_codons(seq)
    rows = []
    for codon in codons:
        amino = codon_table[codon]
        note = ""
        if codon in start_codons:
            note = "Start codon"
        elif codon in stop_codons:
            note = "Stop codon"
        rows.append({"Codon": codon, "Amino acid": amino, "Note": note})

    starts = [c for c in codons if c in start_codons]
    stops = [c for c in codons if c in stop_codons]

    print("Sequence:", seq)
    print("Codons:", " - ".join(codons))
    print("Amino acids:", " - ".join(r["Amino acid"] for r in rows))
    if starts:
        print("Start codon:", ", ".join(starts))
    else:
        print("Start codon: none found")
    if stops:
        print("Stop codon:", ", ".join(stops))
    else:
        print("Stop codon: none found")
    if leftover:
        print("Incomplete codon at the end:", leftover, "(" + str(len(leftover)) + " leftover base(s), not translated)")

    return pd.DataFrame(rows)


translate("ATGGCTTTTTAA")


translate("atg gct ttt taa gc")


translate("ATGXYZ")


import re
import urllib.request

github_url = "https://raw.githubusercontent.com/biopython/biopython/master/Bio/Data/CodonTable.py"

one_letter = {
    "A": "Alanine", "R": "Arginine", "N": "Asparagine", "D": "Aspartic acid",
    "C": "Cysteine", "Q": "Glutamine", "E": "Glutamic acid", "G": "Glycine",
    "H": "Histidine", "I": "Isoleucine", "L": "Leucine", "K": "Lysine",
    "M": "Methionine", "F": "Phenylalanine", "P": "Proline", "S": "Serine",
    "T": "Threonine", "W": "Tryptophan", "Y": "Tyrosine", "V": "Valine",
}

github_table = {}
try:
    text = urllib.request.urlopen(github_url, timeout=20).read().decode()
    start = text.index('name="Standard"')
    end = text.index("register_ncbi_table(", start)
    block = text[start:end]
    for codon, letter in re.findall(r'"([ACGT]{3})": "([A-Z])"', block):
        github_table[codon] = one_letter[letter]
    stop_part = block[block.index("stop_codons"):]
    for codon in re.findall(r'"([ACGT]{3})"', stop_part.split("]")[0]):
        github_table[codon] = "Stop"
    print("Loaded", len(github_table), "codons from GitHub (Biopython)")
except Exception as e:
    print("Could not download the GitHub table:", e)


if github_table:
    compare_rows = []
    matches = 0
    for codon in codon_table:
        mine = codon_table[codon]
        theirs = github_table.get(codon, "missing")
        same = mine == theirs
        if same:
            matches += 1
        compare_rows.append({"Codon": codon, "My table": mine, "GitHub table": theirs, "Match": same})
    compare_df = pd.DataFrame(compare_rows)
    print("Matches:", matches, "out of", len(codon_table))
    print("Mismatches:", len(codon_table) - matches)
    print(compare_df[compare_df["Match"] == False] if matches < 64 else compare_df.head(10))
else:
    print("Skipped comparison because the GitHub table was not loaded.")


rna_table = {}
for codon in codon_table:
    rna_table[codon.replace("T", "U")] = codon_table[codon]

rna_df = pd.DataFrame(list(rna_table.items()), columns=["RNA codon", "Amino acid"])
print(rna_df.to_string())


def dna_to_rna(seq):
    return clean_sequence(seq).replace("T", "U")


def translate_rna(seq):
    seq = clean_sequence(seq)
    if "T" in seq and "U" in seq:
        print("Error: the sequence has both T and U.")
        return None
    seq = seq.replace("T", "U")
    for base in seq:
        if base not in "AUGC":
            print("Error: invalid character", base, "- an RNA sequence can only have A, U, G and C.")
            return None

    codons, leftover = split_codons(seq)
    rows = []
    for codon in codons:
        amino = rna_table[codon]
        note = ""
        if codon == "AUG":
            note = "Start codon"
        elif codon in ["UAA", "UAG", "UGA"]:
            note = "Stop codon"
        rows.append({"Codon": codon, "Amino acid": amino, "Note": note})

    print("RNA sequence:", seq)
    print("Codons:", " - ".join(codons))
    print("Amino acids:", " - ".join(r["Amino acid"] for r in rows))
    if leftover:
        print("Incomplete codon at the end:", leftover, "(not translated)")
    return pd.DataFrame(rows)


print(dna_to_rna("ATGGCTTTTTAA"))
translate_rna("ATGGCTTTTTAA")


translate_rna("AUGGCUUUUUAA")


try:
    user_seq = input("Enter a DNA or RNA sequence (press Enter for ATGGCTTTTTAA): ")
except Exception:
    user_seq = ""

if user_seq.strip() == "":
    user_seq = "ATGGCTTTTTAA"

if "U" in user_seq.upper():
    result = translate_rna(user_seq)
else:
    result = translate(user_seq)

print(result)
