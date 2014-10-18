
# coding: utf-8

# In[1]:

from Bio import Seq,SeqIO


# In[55]:

def score(seq1,seq2):
    diff = 0
    m = len(seq1)
    for i in range(m):
        if seq1[i]!=seq2[i]:
            diff +=1
    return 1.0 * diff / m


# In[62]:

seqs_rec = SeqIO.parse("/Users/richard/Downloads/rosalind_pdst.txt","fasta")


# In[63]:

seqs = [str(seq.seq) for seq in seqs_rec]


# In[64]:

res = [[score(seq_i,seq_j) for seq_j in seqs] for seq_i in seqs]


# In[66]:

with open("/Users/richard/Downloads/output","w") as file:
    for row in res:
        for col in row:
            file.write(str(col)+" ")
        file.write("\n")

