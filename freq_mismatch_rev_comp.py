from Bio import Seq,SeqIO
import time
import itertools

def freq_word_mismatch_both_strand(text,rev_comp_text,k,d):
    patterns = list(make_candidate_fast(text,k,d))
    counts = []
    for pattern in patterns: counts.append(approx_match_count_with_rev_comp(pattern,text,d))
    max_count = max(counts)
    res = [patterns[i] for i in range(len(patterns)) if counts[i] == max_count]
    return res

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

def make_candidate_fast(text,k,d):
    res = set()
    base_patterns = [text[i:k+i] for i in range(len(text)-k+1)]
    for pattern in base_patterns: res |= set(d_neighbor(pattern,k,d))
    return res

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
   
######################################
######################################
############ Running #################

seq_rec = list(SeqIO.parse("/Users/richard/Downloads/text.txt",'fasta'))[0]
seq = seq_rec.seq
seq_rev = seq.reverse_complement()
text = str(seq)
text_rev = str(seq_rev)
k = 8
d = 3

start = time.time()
res = freq_word_mismatch_both_strand(text,text_rev,k,d)
end = time.time()
print ' '.join(res)
print end-start