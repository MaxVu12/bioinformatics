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
    return FrequentPatterns



    
    
   
	
