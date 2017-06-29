# sample count Matrix
count = {"A": [2, 2, 0, 0, 0, 0, 9, 1, 1, 1, 3, 0],
         "C": [1, 6, 0, 0, 0, 0, 0, 4, 1, 2, 4, 6],
         "G": [0, 0,10,10, 9, 9, 1, 0, 0, 0, 0, 0],
         "T": [7, 2, 0, 0, 1, 1, 0, 5, 8, 7, 3, 4]
        } 

#Input: A set of kmers Motifs
#Output: count matrix of Motifs (as a dictionary of lists) 
def Count(Motifs): 
    count = {}
    k = len(Motifs[0])
    for symbol in "ACGT":
        count[symbol] = []
        for j in range(k):
            count[symbol].append(0)
    t = len(Motifs)
    for i in range(t):
        for j in range(k):
            symbol = Motifs[i][j]
            count[symbol][j] += 1 
    return count 

#Input: A list of kmers Motifs
#Output: profile matrix of Motifs (as a dictionary of lists)
def Profile(Motifs):
    t = len(Motifs)
    k = len(Motifs[0])
    profile = Count(Motifs)
    for x in profile: 
        for y in range(k):
            profile[x][y] /= float(t)
    return profile 

#Input: A set of kmers Motifs 
#Output: A consensus string of Motifs (a.k.a the most popular string) 
def Consensus(Motifs):
    k = len(Motifs[0])
    count = Count(Motifs)
    consensus = ""
    for j in range(k):
        m = 0
        frequentSymbol = ""
        for symbol in "ACGT":
            if count[symbol][j] > m:
                m = count[symbol][j]
                frequentSymbol = symbol
        consensus += frequentSymbol 
    return consensus

#Input: A set of k-mers Motifs
#Output: The score of these kmers
def Score(Motifs):
    ConsensusSequence = Consensus(Motifs)
    t = len(Motifs)
    Score = 0
    for n in range(t):
        for m in range(len(ConsensusSequence)):
            if Motifs[n][m]: != ConsensusSequence[m]:
                    Score += 1
    return Score 

#Input: String Text & profile matrix Profile
#Output: The probability of the kmer (which is used to determine if it is close to being the consensus string) 
def Pr(Text, Profile):
    p = 1
    for i in range(len(Text)):
        key = Text[i]
        d = Profile[key][i]
        p = p * d
    return p 

#Input: String Text, an integer k, and profile matrix Profile
#Output: The k-mer that is most likely to be the consensus string based on Profile 
def ProfileMostProbablePattern(Text, k, Profile):
    best_prob = -1
    best_kmer = ""
    for i in range(len(Text)-k+1):
        kmer = Text[i:i+k]
        prob = Pr(kmer, Profile)
        if prob > best_prob:
            best_prob = prob 
            best_kmer = kmer 
    return best_kmer 

#Input: A list of kmers Dna, and integers k and t (t is the no. of kmers in Dna)
#Output: GreedyMotifSearch(Dna, k, t)
def GreedyMotifSearch(Dna, k, t):
    BestMotifs = []
    for i in range(0, t):
        BestMotifs.append(Dna[i][0:k])
    n = len(Dna[0])
    for i in range(n-k+1):
        Motifs = []
        Motifs.append(Dna[0][i:i+k])
        for j in range(1, t):
            P = Profile(Motifs[0:j])
            Motif.append(ProfileMostProbablePattern(Dna[j], k, P))
        if Score(Motifs) < Score(BestMotifs):
            BestMotifs = Motifs
    return BestMotifs

#
            
