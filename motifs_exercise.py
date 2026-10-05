lecture_dna = [
    "TGACGTATAAGTTGCGATGGACGAGATAGCAGAGAATAGGCAACGAGAGATAAGCAG",
    "GACGGTAGCAGATAGACAGATGAAGAGTATGAATTGCACAGATAGCAGATAGCAGAT",
    "GGAGTGTGACGTAGCAGAGACGAAAGACGTAGAGTAGCAGTAGCAGATAGAGGGAGT",
    "TAGACAGTATAGAGACAGCGAGTCGGATAGCACCCAGTATGACGATAGCAATGACAG",
    "GCAGTAGAGCAGATTAGCATTGACAGATAGACGATTGGAGAGATGTGTGGATGACGA",
    "GGCAGGTAGCACACTGGGTCGATAAAGAGTAGCATAGAGACATAGACATATTTTAGC",
]

def count_matrix(motifs):
    l = len(motifs[0])

    counts = {base: [0] * l for base in "ACGT"}
    for motif in motifs:
        for column, base in enumerate(motif):
            if base in motif:
                counts[base][column] += 1
    return counts  

def score(motifs):
    counts = count_matrix(motifs)
    l = len(motifs[0])
    score = 0

    for column in range(l):
        column_counts = [counts[base][column] for base in "ACGT"]
        max_count = max(column_counts)
        score += max_count
    return score

def consensus(motifs):
    counts = count_matrix(motifs)
    l = len(motifs[0])
    consensus = ""

    for column in range(l):
        column_counts = {base: counts[base][column] for base in "ACGT"}
        max_base = max(column_counts, key=column_counts.get)
        consensus += max_base
    return consensus

def hamming_distance(seq1, seq2):
    return sum(base1 != base2 for base1, base2 in zip(seq1, seq2))

def total_distance(pattern, sequences):
    l = len(pattern)
    total_distance = 0

    for seq in sequences:
        min_distance = float('inf')
        for i in range(len(seq) - l + 1):
            window = seq[i:i+l]
            distance = hamming_distance(pattern, window)
            if distance < min_distance:
                min_distance = distance
        total_distance += min_distance
    return total_distance

red = ["TAAGTT", "TGAATT", "GGAGTG", "CGAGTC", "TGTGTG", "TGGGTC"]  # slide 19
best = ["AGATAG", "AGATAG", "AGATAG", "AGACAG", "AGATAG", "AGGTAG"]

print(score(red))                                # 26
print(consensus(best), score(best))              # AGATAG 34
print(hamming_distance("TAAGTT", "TGAATT"))      # 2
print(total_distance("TGCGTT", lecture_dna))     # 13