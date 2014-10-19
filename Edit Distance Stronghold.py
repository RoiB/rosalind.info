
# coding: utf-8

# In[69]:

'''
implement edit distance
'''
def edit_dist(seq1,seq2):
    
    #initialize table
    m = len(seq1)
    n = len(seq2)
#     print m,n
#     table = [[float('inf') for j in range(n+1)] for i in range(m+1)]
    table = [[0 for j in range(n+1)] for i in range(m+1)]
    table[0][0] = 0
    for i in range(1,m+1):
        table[i][0] = i
    for j in range(1,n+1):
        table[0][j] = j
    
    #implement algorithms
    def diff(i,j):
        if seq1[i-1] == seq2[j-1]:
            return 0
        else:
            return 1
    def OPT(i,j): return min(table[i-1][j]+1, table[i][j-1]+1,diff(i,j)+ table[i-1][j-1])
    
    #fill up the whole table
    for i in range(1,m+1):
        for j in range(1,n+1):
            table[i][j] = OPT(i,j)
#     for row in table: print row
    return table[m][n]


# In[83]:

from Bio import Seq,SeqIO


# In[86]:

seq_recs = SeqIO.parse('/Users/richard/Downloads/rosalind_edit (1).txt','fasta')
seq1,seq2 = [str(seq.seq) for seq in seq_recs]


# In[87]:

edit_dist(seq1,seq2)

