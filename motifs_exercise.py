from itertools import count


lecture_dna = [
    "TGACGTATAAGTTGCGATGGACGAGATAGCAGAGAATAGGCAACGAGAGATAAGCAG",
    "GACGGTAGCAGATAGACAGATGAAGAGTATGAATTGCACAGATAGCAGATAGCAGAT",
    "GGAGTGTGACGTAGCAGAGACGAAAGACGTAGAGTAGCAGTAGCAGATAGAGGGAGT",
    "TAGACAGTATAGAGACAGCGAGTCGGATAGCACCCAGTATGACGATAGCAATGACAG",
    "GCAGTAGAGCAGATTAGCATTGACAGATAGACGATTGGAGAGATGTGTGGATGACGA",
    "GGCAGGTAGCACACTGGGTCGATAAAGAGTAGCATAGAGACATAGACATATTTTAGC",
]

#TASK 1
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

#print(score(red))                                # 26
#print(consensus(best), score(best))              # AGATAG 34
#print(hamming_distance("TAAGTT", "TGAATT"))      # 2
#print(total_distance("TGCGTT", lecture_dna))     # 13

#TASK 2
class MotifProfile:
    def __init__(self, motifs, pseudocount=1):
        self.motifs = motifs
        self.pseudocount = pseudocount
        self.l = len(motifs[0])
        t = len(motifs)
        counts = count_matrix(motifs)
        self.ppm = {base: [(count + pseudocount) / (t + 4 * pseudocount) for count in column] for base, column in counts.items()}

#profile = MotifProfile(["ATCCGTA", "GTGCATA", "AAGCGTA", "ATGCGTG"])
#print(profile.l)         # 7
#print(profile.ppm["A"])    # [0.5, 0.25, 0.125, 0.125, 0.25, 0.125, 0.5]

    def lmer_probability(self, lmer):
        prob = 1
        for i, base in enumerate(lmer):
            prob *= self.ppm[base][i]
        return prob
    
    def most_probable_lmer(self, sequence):
        best_lmer = None
        best_prob = -1
        for i in range(len(sequence) - self.l + 1):
            window = sequence[i:i+self.l]
            lmer = self.lmer_probability(window)
            if lmer > best_prob:
                best_prob = lmer
                best_lmer = window
        return best_lmer

    def consensus(self):
        consensus = ""
        for column in range(self.l):
            column_counts = {base: self.ppm[base][column] for base in "ACGT"}
            max_base = max(column_counts, key=column_counts.get)
            consensus += max_base
        return consensus

#profile = MotifProfile(["ATCCGTA", "GTGCATA", "AAGCGTA", "ATGCGTG"])
#print(profile.consensus())                       # ATGCGTA
#print(round(profile.lmer_probability("ATGCGTA"), 4))  # 0.0122

#two = MotifProfile(["GTAC", "TTAA"])
#print(two.most_probable_lmer("ACTGGATGACCC"))    # TGAC
#print(round(two.lmer_probability("TGAC"), 4))         # 0.0093

#from Bio import motifs
#from Bio.Seq import Seq

#bio = motifs.create([Seq(site) for site in ["ATCCGTA", "GTGCATA", "AAGCGTA", "ATGCGTG"]])
#bio.pseudocounts = 1
#print(bio.consensus)     # ATGCGTA
#print(bio.pwm["A"])      # the same numbers as your profile.ppm["A"]

