# Input: Strings Pattern & Text 
# Output: No. of times Pattern appears in Text 

def PatternCount(Pattern, Text): 
    count = 0 
    for i in range(len(Text)-len(Pattern)+1):
        if Text[i:i+len(Pattern)] == Pattern: 
      		count = count + 1
    return count 

# Input: A string Text & an integer k 
# Output: Count dictionary for k-mer in Text
def CountDict(Text, k): 
    Count = {}
    for i in range(len(Text)-k+1):
        Pattern = Text[i:i+k]
        Count[i] = PatternCount(Pattern, Text)
    return Count 

# Input: A string Text & an integer k
# Output: A list containing all most frequent k-mers in Text 
def FrequentWords(Text, k):
    FrequentPatterns = []
    Count = CountDict(Text, k):
    m = max(Count.values())
    for i in Count:
        if Count[i] == m:
            FrequentPatterns.append(Text[i:i+k])
    FrequentPatternNoDuplicates = remove_duplicates(FrequentPatterns)
	return FrequentPatternsNoDuplicates 

#Input: A list Items
#Output: A list containing all objects from Items without duplicates 
def remove_duplicates(Items):
	ItemsNoDuplicates = []
	for x in Items:
		if x not in ItemsNoDuplicates:
			ItemNoDuplicates.append(x)
	return ItemNoDuplicates

#Input: A string Nucleotide
#Output: the complement of Nucleotide
def complement(Nucleotide):
	comp = {'A': 'T','T': 'A', 'G': 'C', 'C': 'G'}
	bases = list(Nucleotide)
	bases = [comp[base] for base in bases]
	return ''.join(bases)

#Input: A string Nucleotide
#Output: The reverse of Nucleotide
def reverse(Nucleotide):
	rev = ""
	for i in Nucleotide:
		rev = i + rev 
	return rev

#Input: A DNA string Pattern 
#Output: The reverse complement of Pattern
def ReverseComplement(Pattern):
	revComp = reverse(complement(Pattern))
	return revComp 

#Input: Two Sstrings, Pattern and Genome
#Output: A list containing all starting positions where Pattern appears as a substring of Genome
def PatternMatching(Pattern, Genome):
	positions = []
	for i in range(len(Genome)-len(Pattern)+1):
		if Genome[i:i+len(Pattern)] == Pattern:
			positions.append(i)
	return positions

#Input: Strings Genome and symbol
#Output: The symbol array of Genome corresponding to symbol
def SymbolArray(Genome, symbol):
	array = {}
	n = len(Genome)
	ExtendedGenome = Genome + Genome[0:n//2]
	for i in range(n):
		array[i] = PatternCount(symbol, ExtendedGenome[i:i(n//2)])
	return array 


#Input: Strings Genome and symbol
#Output: The symbol array of Genome corresponding to symbol
def FasterSymbolArray(Genome, symbol):
	array = {}
	n = len(Genome)
	ExtendedGenome = Genome + Genome[0:n//2]
	array[0] = PatternCount(symbol, Genome[0:n//2])
	for i in range(1, n):
		array[i] = array[i-1]
		if ExtendedGenome[i-1] == symbol:
			array[i] = array[i] - 1
		if ExtendedGenome[i+(n//2)-1] == symbol:
			array[i] = array[i] + 1
	return array

#Input: A string Genome
#Output: the skew array of Genome in the form of a dictionary mapping the i-th symbol of Genome to Skew[i]
def Skew(Genome):
	skew = {}
	n = len(Genome)
	skew[0] = 0 
	for i in range(1, n+1):
		skew[i] = skew[i-1]
		if Genome[i-1] == 'C':
			skew[i] = skew[i] - 1
		if Genome[i-1] == 'G':
			skew[i] = skew[i] + 1
		if Genome[i-1] == 'A' or 'T':
			skew[i] = skew[i]
	return skew

#Input: A DNA string Genome
#Output: A list containing all integers i minimizing Skew(Prefix_i(Text)) over all values of i (from 0 to |Genome|)
def MinimumSkew(Genome):
	positions = []
	askew = Skew(Genome)
	mini = min(askew.values())
	for i in askew:
		if askew[i] == mini:
			pos = i
			position.append(pos)
	return positions

#Input: 2 strings p & q
#Output: An integer value representing the Hamming Distance between p & q (which is the number of mismatches at position i)
def HammingDistance(p, q):
	count = 0
	for i in range(0, len(q)):
		if p[i] != q[i]:
			count += 1
	return count 

#Input: Strings Patterns and Text along with an integer d 
#Output: A list containing all starting positions where Pattern appears
def ApproximatePatternMatching(Pattern, Text, d):
	positions = [] 
	for i in range(0, len(Text) - len(Pattern)+1):
		ham = HammingDistance(Pattern, Text[i:i+len(Pattern)])
		if ham <= d:
			positions.append(i)
	return positions

#Input: Strings Patterns and Text, and an integer d
#Output: The number of times Pattern appears in Text with at most d mismatches
def ApproximatePatternCount(Pattern, Text, d):
	count = 0
	for i in range(0, len(Text)-len(Pattern)+1):
		ham = HammingDistance(Pattern, Text[i:i+len(Pattern)])
		if ham <= d:
			count += 1 
	return count 






    
    
   
	
