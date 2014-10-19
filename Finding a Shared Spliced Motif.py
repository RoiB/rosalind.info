
# coding: utf-8

# In[32]:

import seq_align
from Bio import SeqIO,Seq


# In[22]:

score_mat = seq_align.build_scoring_matrix(set(["A","T","G","C"]),10,-10,0)


# In[33]:

seq_recs = SeqIO.parse("/Users/richard/Downloads/rosalind_lcsq.txt","fasta")
seq_x,seq_y = [str(seq.seq) for seq in seq_recs]


# In[34]:

align_mat = seq_align.compute_alignment_matrix(seq_x, seq_y, score_mat, True)


# In[35]:

_,align1,align2 = seq_align.compute_global_alignment(seq_x,seq_y,score_mat,align_mat)
print align1
print align2


# In[36]:

res = ''.join([align1[i] for i in range(len(align1)) if align1[i]!= '-' and align2[i]!='-'])


# In[37]:

with open ("/Users/richard/Downloads/output",'w') as file:
    file.write(res)

