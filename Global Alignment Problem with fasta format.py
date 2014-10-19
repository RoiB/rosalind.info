
# coding: utf-8

# In[1]:

import BLOSUM62
from Bio import SeqIO,Seq


# In[2]:

def build_scoring_matrix(mat,indel):
    new_mat = {amino1:{amino2:mat[amino1][amino2] for amino2 in mat} for amino1 in mat}
    for amino in new_mat:
        new_mat[amino]['-'] = indel
    new_mat['-'] = {amino:indel for amino in mat}
    return new_mat


# In[3]:

def compute_alignment_matrix(seq_x, seq_y, scoring_matrix):
    '''
    The function computes and returns the alignment matrix for seq_x and seq_y as described in the Homework.
    '''
    len_x = len(seq_x)
    len_y = len(seq_y)
    score = [[None for _ in range(len_y+1)] for _ in range(len_x+1)]
    score[0][0] = 0
    
    for id_i in range(1,len_x+1):
        score[id_i][0] =  score[id_i-1][0] + scoring_matrix[seq_x[id_i-2]]['-']
            
    for id_j in range(1,len_y+1):
        score[0][id_j] =  score[0][id_j-1] + scoring_matrix['-'][seq_y[id_j-2]]
        
    for id_i in range(1,len_x+1):
        for id_j in range(1,len_y+1):
            score[id_i][id_j] =  max(score[id_i-1][id_j]+scoring_matrix[seq_x[id_i-1]]['-'],                                     score[id_i][id_j-1]+scoring_matrix['-'][seq_y[id_j-1]],                                     score[id_i-1][id_j-1] + scoring_matrix[seq_x[id_i-1]][seq_y[id_j-1]])
    return score


# In[4]:

def readfile(filename):
    with open(filename) as file:
        return file.readline().rstrip(), file.readline().rstrip()


# In[5]:

def compute_global_alignment(seq_x, seq_y, scoring_matrix, alignment_matrix):
    '''
    Computes a global alignment of seq_x and seq_y using the global alignment matrix alignment_matrix.
    '''
    id_i = len(seq_x)
    id_j = len(seq_y)
    align_x = align_y = ''
    while id_i != 0 and id_j != 0:
        if alignment_matrix[id_i][id_j] == alignment_matrix[id_i-1][id_j-1] + scoring_matrix[seq_x[id_i-1]][seq_y[id_j-1]]:
            align_x = seq_x[id_i-1] + align_x
            align_y = seq_y[id_j-1] + align_y
            id_i -= 1
            id_j -= 1
        else:
            if alignment_matrix[id_i][id_j] == alignment_matrix[id_i-1][id_j] + scoring_matrix[seq_x[id_i-1]]['-']:
                align_x = seq_x[id_i-1] + align_x
                align_y = '-' + align_y
                id_i -= 1
            else:
                align_x = '-' + align_x
                align_y = seq_y[id_j-1] + align_y
                id_j -= 1
    while id_i != 0:
        align_x = seq_x[id_i-1] + align_x
        align_y = '-'+align_y
        id_i -= 1
    
    while id_j != 0:
        align_y = seq_y[id_j-1] + align_y
        align_x = '-'+align_x
        id_j -= 1
    return alignment_matrix[len(seq_x)][len(seq_y)],align_x,align_y


# In[18]:

mat = build_scoring_matrix(BLOSUM62.BLOSUM62,-5)
seq1,seq2 = SeqIO.parse("/Users/richard/Downloads/rosalind_glob-2.txt",'fasta')
align_mat = compute_alignment_matrix(seq1,seq2,mat)
res = compute_global_alignment(seq1, seq2, mat, align_mat)


# In[19]:

with open('/Users/richard/Downloads/output','w') as file:
        file.write(str(res[0])+"\n")

