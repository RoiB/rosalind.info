
# coding: utf-8

# In[41]:

from Bio import Seq, SeqIO


# In[53]:

seq = list(SeqIO.parse("/Users/richard/Downloads/rosalind_4a.txt","fasta"))[0]


# In[44]:

def sub_str(seq,k):
    seq = str(seq.seq)
    res = [seq[i:i+k] for i in range(len(seq)-k+1)]
    return sorted(res)


# In[54]:

res = sub_str(seq,100)


# In[55]:

with open("/Users/richard/Downloads/output","w") as file:
    for item in res: file.write(item+"\n")

