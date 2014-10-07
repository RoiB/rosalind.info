
# coding: utf-8

# In[16]:

from Bio import SeqIO
from Bio import Seq


# In[136]:

seq_records = [_ for _ in SeqIO.parse("/Users/richard/Downloads/rosalind_long.txt",'fasta')]


# In[137]:

seq_strings = [str(seq.seq) for seq in seq_records]


# In[138]:

def right_merge(result, seq_string):
    if seq_string in result: return result
    mid =  len(seq_string)//2
    if seq_string[:mid] in result:
        while seq_string[:mid] in result: mid += 1
        return result + seq_string[mid-1:]
    else:
        return False


# In[139]:

def left_merge(result,seq_string):
    if seq_string in result: return result
    mid =  len(seq_string)//2
    if seq_string[mid:] in result:
        while seq_string[mid:] in result: mid -= 1
        return seq_string[:mid+1]+result
    else:
        return False


# In[140]:

def merge(result,seq_string):
    left_merged = left_merge(result,seq_string)
    if not left_merged:
        return right_merge(result, seq_string)
    else:
        right_merged = right_merge(result,seq_string)
        if right_merged: 
            return min(right_merged,left_merged,key = len)
        else:
            return left_merged


# In[142]:

result, other_strings = seq_strings[0], seq_strings[1:]
while other_strings:
    merged_indices = []
    n = len(other_strings)
    for i in range(n):
        temp =  merge(result,other_strings[i])
        if temp:
            result = temp
            merged_indices.append(i)
    while merged_indices:
        other_strings.pop(merged_indices.pop())


# In[144]:

with open("/Users/richard/Downloads/output","w") as file:
    for i in result: file.write(i)

