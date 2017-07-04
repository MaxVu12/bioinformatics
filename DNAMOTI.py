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
            if Motifs[n][m] != ConsensusSequence[m]:
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

#Input: A set of kmers Motifs
#Output: Add 1 to every element in Motifs
def CountWithPseudocounts(Motifs): 
    count = {}
    t = len(Motifs)
    k = len(Motifs[0])
    for symbol in "ACGT": 
        count[symbol] = []
        for j in range(k):
            count[symbol].append(1)
    for i in range(t):
		for j in range(k):
            symbol = Motifs[i][j]
            count[symbol][j] += 1 
	return count 

#Input: A set of kmers Motifs 
#Output: The Profile matrix adjusted with Laplace's Rule of Succession
def ProfileWithPseudocounts(Motifs):
	t = len(Motifs)
	k = len(Motifs[0])
	profile = CountWithPseudocounts(Motifs)
	for x in profile:
		for y in range(k):
			profile[x][y] /= float(t+4) 
	return profile

#Input: A list of kmers Dna, and integers k and t (where t is the number of kmers in Dna)
#Output: The consensus string from each string of Dna using greedy search
def GreedyMotifSearchWithPseudocounts(Dna, k, t):
	BestMotifs = []
	for i in range(0, t):
		BestMotifs.append(Dna[i][0:k])
	n = len(Dna[0])
	for i in range(n-k+1):
		Motifs = []
		Motif.append(Dna[0][i:i+k])
		for j in range(1, t):
			P = ProfileWithPseudocounts(Motifs[0:j])
			Motif.append(ProfileMostProbablePattern(Dna[j], k, P))
		if Score(Motifs) < Score(BestMotifs):
	return BestMotifs

#Input: A profile matrix Profile and a list of strings Dna
#Output: A motif from Dna using the Profile
def Motif(Profile, Dna):
	motiflist = []
	k = len(Profile['A'])
	for i in range(len(Dna)):
		mostprobable = ProfileMostProbablePattern(Dna[i], k, Profile)
		motiflist.append(mostprobable)
	return motiflist

#Input:a list of strings Dna and integers k and t
#Output: A random motif
#t is essentially len(Dna). Also, remember to import random
def RandomMotifs(Dna, k, t):
	motiflist = []
	for i in range(t):
		s = random.randint(0, t-1) 
		motiflist.append(Dna[i][s:s+k])
	return motiflist

#Input: Positive integers k and t, followed by a list of strings Dna
#Output: A random motif (the one with the best score)
def RandomizedMotifSearch(Dna, k, t):
	M = RandomMotifs(Dna, k, t)
	BestMotifs = M 
	while True:
		Profile = ProfileWithPseudocounts(M)
		M = Motifs(Profile, Dna)
		if Score(M) < Score(BestMotifs):
			BestMotifs = M
		else: 
			return BestMotifs

#Input: A dictionary Probabilites, where keys are k-mers and values are the probabilities of these k-mers (which do not necessary sum up to 1)
#Output: A normalized dictionary where the probability of each k-mer was divided by the sum of all k-mers' probabilities
def Normalize(Probabilities):
	total = 0
	for vals in Probabilities.values():
		total += vals
	for k in Probabilities.keys():
		Probabilities[k] = Probabilities[k]/total
	return Probabilities 

#Input: A dictionary Probabilities whose keys are k-mers and whose values are the probabilities of these kmers
#Output: A randomly chosen k-mer with respect to the values in Probabilities 
#Remember to import random
def WeightedDie(Probabilities):
	kmer = ''
	r = random.uniform(0, 1)
	prob = 0.0
	for key, va in Probabilities.items():
		prob += float(va)
		if float(r) <= float(prob):
			kmer = key
			return kmer 

#Input: A string Text, a profile matrix Profile, and an integer k 
#Output: A randomly generated kmer from Text whose probabilities are generated from profile 
def ProfileGeneratedString(Text, profile, k):
	n = len(Text)
	probabilities = {}
	for i in range(0, n-k+1):
		probabilities[Text[i:i+k]] = Pr(Text[i:i+k], profile)
	probabilities = Normalize(probabilities)
	return WeightedDie(probabilities)

# To-do
#     GibbsSampler(Dna, k, t, N)
#         randomly select k-mers Motifs = (Motif1, …, Motift) in each string from Dna 
#         ﻿BestMotifs ← Motifs
#         for j ← 1 to N
#             i ← randomly generated integer between 1 and t
#             Profile ← profile matrix formed from all strings in Motifs except for Motifi
#             Motifi ← Profile-randomly generated k-mer in the i-th string
#             if Score(Motifs) < Score(BestMotifs) 
#                 BestMotifs ← Motifs
#         return BestMotifs

#Input: A list of strings Dna, an integer k that represents the size of k-mer, an integer t that will be used to generate a random number and N, the number of time this algorithm is runned
#Output: The best possible Motifs from Dna 
def GibbsSampler(Dna, k, t, N):
	motif = RandomMotifs(Dna, k, len(Dna))
	BestMotifs = motif  
	for j in range(N):
		i = random.randint(0, t-1)
		new_motif = motif[:i] + motif[i+1:]
		profile = ProfileWithPseudocounts(new_motif)
		motif[i] = ProfileGeneratedString(Dna[i], profile, k)
		if Score(motif) < Score(BestMotifs):
			BestMotifs = motif
	return BestMotifs

#For repeating GibbsSampler 
def RepeatedGibbsSampler(Dna, k, t, N):
    BestScore = float('inf')
    BestMotifs = []
    for i in range(100):
        Motifs = GibbsSampler(Dna, k, t, N)
        CurrScore = Score(Motifs)
        if CurrScore < BestScore:
            BestScore = CurrScore
            BestMotifs = Motifs
    return BestMotifs

		


            
