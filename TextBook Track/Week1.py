
# coding: utf-8

# In[1]:

from Bio import Seq,SeqIO
import time
import itertools


# In[35]:

def pattern_count(text,pattern):
    count = 0
    for i in range(len(text)-len(pattern)+1):
        if text[i:i+len(pattern)] == pattern:
            count += 1
    return count


# In[36]:

def freq_word(text,k):
    '''
    Solve the Frequent Words Problem.
    Input: A string Text and an integer k.
    Output: All most frequent k-mers in Text.
    '''
    freq_table = {}
    for i in range(len(text)-k+1):
        pattern = text[i:i+k]
        if pattern in freq_table:
            freq_table[pattern]+=1
        else:
            freq_table[pattern] = 1
    max_count = max(freq_table.values())
    res = ''
    for pattern in freq_table:
        if freq_table[pattern] == max_count:
            res += pattern+" "
    return res


# In[55]:

def rev_complement(pattern):
    rec = Seq.Seq(pattern)
    return str(rec.reverse_complement())


# In[61]:

def pattern_matching(pattern,text):
    '''
    Find all occurrences of a pattern in a string.
    Input: Two strings, Pattern and Genome.
    Output: All starting positions where Pattern appears as a substring of Genome.
    '''
    res = []
    for i in range(len(text)-len(pattern)+1):
        if text[i:i+len(pattern)] == pattern:
            res.append(i)
    return res


# In[79]:

pattern = 'CTTGATCAT'


# In[82]:

### fast Biopython package fasta process 
seq_rec = list(SeqIO.parse("/Users/richard/Downloads/Vibrio_cholerae.txt",'fasta'))[0]
seq = seq_rec.seq
genome = str(seq)


# In[83]:

result = pattern_matching(pattern,genome)
' '.join([str(num) for num in result])


# In[ ]:

###
###
### The Clump Finding Problem ###


# In[46]:

def in_range(alist,L,t): 
    if len(alist) < t: return False
    max_range = alist[-1] - alist[0]
    if max_range < L:
        return True
    else:
        return in_range(alist[:-1],L,t)        


# In[49]:

def clump_finding(genome,k, L, t):
    freq_table = {}
    loc_table = {}
    for i in range(len(genome)-k+1):
        pattern = genome[i:i+k]
        if pattern in freq_table:
            freq_table[pattern]+=1
            loc_table[pattern].append(i)
        else:
            freq_table[pattern] = 1
            loc_table[pattern] = [i]
    filtered_table = {pattern:loc_table[pattern] for pattern in loc_table if freq_table[pattern] >= t}
    res = []
    for pattern in filtered_table:
        if in_range(filtered_table[pattern],L,t):
            res.append(pattern)
    return res


# In[70]:

# ##test
genome = 'CCGACAGGCTAGTCTATAATCCTGAGGCGTTACCCCAATACCGTTTACCGTGGGATTTGCTACTACAACTCCTGAGCGCTACATGTACGAAACCATGTTATGTAT'
k = 4
L = 30
t = 3


# In[71]:

### E-coli
### fast Biopython package fasta process 
# seq_rec = list(SeqIO.parse("/Users/richard/Downloads/E-coli.txt",'fasta'))[0]
# seq = seq_rec.seq
# genome = str(seq)
len(genome)


# In[72]:

result = clump_finding(genome,k, L, t)
print result
' '.join(result)


# In[40]:

### Minimum Skew Problem ###
def skewness(seq):
    res = [0]
    count_C = 0
    count_G = 0
    for i in range(len(seq)):
        if seq[i] == 'G':
            count_G += 1
        if seq[i] == 'C':
            count_C += 1
        res.append(count_G - count_C)
    return res


# In[44]:

# seq_rec = list(SeqIO.parse("/Users/richard/Downloads/dataset_7_6.txt",'fasta'))[0]
# seq = seq_rec.seq
# genome = str(seq)


# In[60]:

skew_score = skewness('GATACACTTCCCAGTAGGTACTG')
min_val = min(skew_score)
res = [i for i in range(len(skew_score)) if skew_score[i] == min_val]
res


# In[69]:

skew_score = skewness('CATTCCAGTACTTCATGATGGCGTGAAGA')
min_val = max(skew_score)
res = [i for i in range(len(skew_score)) if skew_score[i] == min_val]
res


# In[37]:

####
#### Approximate Pattern Matching Problem #####
def hamming_dist(str1,str2):
    '''
    Compute the Hamming distance between two strings.
    Input: Two strings of equal length.
    Output: The Hamming distance between these strings.
    '''
    score = 0
    for i in range(len(str1)):
        if str1[i] != str2[i]:
            score += 1
    return score


# In[68]:

hamming_dist('CTTGAAGTGGACCTCTAGTTCCTCTACAAAGAACAGGTTGACCTGTCGCGAAG','ATGCCTTACCTAGATGCAATGACGGACGTATTCCTTTTGCCTCAACGGCTCCT')


# In[54]:

def approx_match(pattern, text, d):
    '''
    Find all approximate occurrences of a pattern in a string.
    Input: Strings Pattern and Text along with an integer d.
    Output: All starting positions where Pattern appears as a substring of Text with at most d mismatches.
    '''
    res = []
    count = 0
    for i in range(len(text)-len(pattern)+1):
        if hamming_dist(text[i:i+len(pattern)],pattern) <= d:
            res.append(i)
            count+=1
    return res,count


# In[143]:

' '.join([str(i) for i in approx_match(pattern, text, d)[0]])


# In[55]:

def count(text, pattern,d):
    '''
    Input: Strings Text and Pattern as well as an integer d.
    Output: Countd(Text, Pattern)
    '''
    return approx_match(pattern, text, d)[1]


# In[73]:

pattern = 'CCC'
text = 'CATGCCATTCGCATTGTCCCAGTGA'
d = 2


# In[74]:

count(text, pattern,d)


# In[150]:

##########
#######
#####  Frequent Words with Mismatches Problem  ###### 


# In[154]:

def freq_word_mismatch(text,k,d):
    patterns = list(make_candidate_fast(text,k,d))
    counts = []
    for pattern in patterns: counts.append(approx_match_count(pattern,text,d))
    max_count = max(counts)
    res = [patterns[i] for i in range(len(patterns)) if counts[i] == max_count]
    return res


# In[155]:

def make_candidate(k):
    return map(''.join, itertools.product('ATGC',repeat = k))


# In[156]:

def make_candidate_fast(text,k,d):
    res = set()
    base_patterns = [text[i:k+i] for i in range(len(text)-k+1)]
    for pattern in base_patterns: res |= set(d_neighbor(pattern,k,d))
    return res


# In[157]:

def d_neighbor(pattern,k,d):
    '''
    all substitutes with at most k mismatches
    a much better version of helper function to make_candidate
    '''
    result = set([pattern])
    template = [letter for letter in pattern]
    for i in range(1,1+d):
        sub_locs = list(itertools.combinations(range(k),i))
        sub_symbols = list(itertools.product('ATGC',repeat = i))        
        for loc in sub_locs:
            temp = template[:]
            for symbol in sub_symbols:
                for j in range(i):
                    temp[loc[j]] = symbol[j]
                result.add(''.join(temp))
    return result


# In[158]:

def hamming_dist_opt(str1,str2,d):
    '''
    Compute the min of Hamming distance between two strings and d.
    Input: Two strings of equal length and max dist allowed
    Output: The min of Hamming distance between these strings and d.
    '''
    score = 0
    for i in range(len(str1)):
        if str1[i] != str2[i]:
            score += 1
            if score > d:
                return d+1
    return score


# In[159]:

def approx_match_count(pattern, text, d):
    count = 0
    for i in range(len(text)-len(pattern)+1):
        if hamming_dist_opt(text[i:i+len(pattern)],pattern,d) <= d:
            count+=1
    return count


# In[ ]:

seq_rec = list(SeqIO.parse("/Users/richard/Downloads/rosalind_1g.txt",'fasta'))[0]
seq = seq_rec.seq
text = str(seq)
k = 9
d = 2


# In[168]:

start = time.time()
res = freq_word_mismatch(text,k,d)
end = time.time()
print ' '.join(res)
print end-start


# In[ ]:

##########
#######
#####  Frequent Words with Mismatches and Reverse Complements Problem  ###### 


# In[180]:

def freq_word_mismatch_both_strand(text,rev_comp_text,k,d):
    patterns = list(make_candidate_fast(text,k,d))
    counts = []
    for pattern in patterns: counts.append(approx_match_count_with_rev_comp(pattern,text,d))
    max_count = max(counts)
    res = [patterns[i] for i in range(len(patterns)) if counts[i] == max_count]
    return res


# In[181]:

def approx_match_count_with_rev_comp(pattern, text, d):
    count = 0
    count_rev = 0
    rev_pattern = str(Seq.Seq(pattern).reverse_complement())
    for i in range(len(text)-len(pattern)+1):
        if hamming_dist_opt(text[i:i+len(pattern)],pattern,d) <= d:
            count+=1
    for i in range(len(text)-len(pattern)+1):
        if hamming_dist_opt(text[i:i+len(pattern)],rev_pattern,d) <= d:
            count_rev+=1        
    return count+count_rev


# In[186]:

seq_rec = list(SeqIO.parse("/Users/richard/Downloads/dataset_9_8.txt",'fasta'))[0]
seq = seq_rec.seq
seq_rev = seq.reverse_complement()
text = str(seq)
text_rev = str(seq_rev)
k = 10
d = 3


# In[187]:

start = time.time()
res = freq_word_mismatch_both_strand(text,text_rev,k,d)
end = time.time()
print ' '.join(res)
print end-start

