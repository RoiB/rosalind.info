
# coding: utf-8

# In[2]:

from Bio import Seq,SeqIO


# In[124]:

'''
implement edit distance
'''
def edit_dist(seq1,seq2):
    
    #initialize table
    m = len(seq1)
    n = len(seq2)

    table = [[0 for j in range(n+1)] for i in range(m+1)]
    parent = [[(0,0) for j in range(n+1)] for i in range(m+1)]
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
    def OPT(i,j): 
        min_val = min(table[i-1][j]+1, table[i][j-1]+1,diff(i,j)+ table[i-1][j-1])
        
        if min_val == table[i-1][j]+1: 
            parent[i][j] = (i-1,j)
        elif min_val == table[i][j-1]+1:
            parent[i][j] = (i,j-1)
        else:
            parent[i][j] = (i-1,j-1)
        return min_val
    
    #fill up the whole table
    for i in range(1,m+1):
        for j in range(1,n+1):
            table[i][j] = OPT(i,j)
#     for row in table: print row
    return table[m][n],parent


# In[134]:

seq_recs = SeqIO.parse('/Users/richard/Downloads/rosalind_edta.txt','fasta')
seq1,seq2 = [str(seq.seq) for seq in seq_recs]
score,parent = edit_dist(seq1,seq2)


# In[135]:

def reconstruct(seq_x,seq_y,parent):
    align_x = align_y = ''
    i, j = len(seq_x), len(seq_y)
    score_mat = parent
    
    while i != 0 and j != 0:
        if score_mat[i][j] == (i-1,j-1):
            align_x = seq_x[i-1]+align_x
            align_y = seq_y[j-1]+align_y
            i -= 1; j-= 1
        elif score_mat[i][j] == (i-1,j):
            align_x = seq_x[i-1]+align_x
            align_y = '-'+align_y
            i-=1
        else: 
            align_x = '-'+align_x
            align_y = seq_y[j-1]+align_y
            j-= 1
        
    while i != 0:
        align_x = seq_x[i-1]+align_x
        align_y = '-'+align_y
        i-=1
        
    while j != 0:
        align_x = '-'+align_x
        align_y = seq_y[j-1]+align_y
        j-= 1
    return align_x,align_y


# In[136]:

align_x,align_y = reconstruct(seq1,seq2,parent)
print align_x
print align_y


# In[137]:

with open("/Users/richard/Downloads/output",'w') as file:
    file.write(str(score)+"\n")
    file.write(align_x+"\n")
    file.write(align_y+"\n")

