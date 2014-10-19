
# coding: utf-8

# In[1]:

import PAM250


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
        score[id_i][0] =  0
            
    for id_j in range(1,len_y+1):
        score[0][id_j] =  0
        
    for id_i in range(1,len_x+1):
        for id_j in range(1,len_y+1):
            score[id_i][id_j] =  max(score[id_i-1][id_j]+scoring_matrix[seq_x[id_i-1]]['-'],                                     score[id_i][id_j-1]+scoring_matrix['-'][seq_y[id_j-1]],                                     score[id_i-1][id_j-1] + scoring_matrix[seq_x[id_i-1]][seq_y[id_j-1]],0)
    return score


# In[4]:

def compute_local_alignment(seq_x, seq_y, scoring_matrix, alignment_matrix):
    '''
    Computes a local alignment of seq_x and seq_y using the local alignment matrix alignment_matrix.
    '''
    id_i, id_j = find_max_idx(alignment_matrix)
    max_i, max_j = id_i, id_j
    
    align_x = align_y = ''
    while id_i != 0 and id_j != 0:
        if alignment_matrix[id_i][id_j] == 0:
            break
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
    return alignment_matrix[max_i][max_j],align_x,align_y


# In[5]:

def find_max_idx(mat):
    '''
    helper function, help to find max value/location in alignment matrix
    '''
    max_val = -1
    max_idx = (-1,-1)
    for id_i in range(len(mat)):
        for id_j in range(len(mat[0])):
            if mat[id_i][id_j] > max_val:
                max_val = mat[id_i][id_j]
                max_idx = (id_i,id_j)    
    return max_idx


# In[6]:

# def readfile(filename):
#     with open(filename) as file:
#         return file.readline().rstrip(), file.readline().rstrip()


# In[7]:

# seq1,seq2 = readfile("/Users/richard/Downloads/rosalind_5f.txt")
# align_mat = compute_alignment_matrix(seq1,seq2,score_mat)


# In[11]:

'''
fasta format version
'''
from Bio import SeqIO,Seq
seq1,seq2 = SeqIO.parse("/Users/richard/Downloads/rosalind_loca-2.txt",'fasta')
score_mat = build_scoring_matrix(PAM250.PAM250,-5)
align_mat = compute_alignment_matrix(seq1,seq2,score_mat)


# In[12]:

score,align_x,align_y = compute_local_alignment(seq1, seq2, score_mat, align_mat)
# print align_x
# print align_y


# In[13]:

align_x = ''.join([letter for letter in align_x if letter is not '-'])
print align_x
align_y = ''.join([letter for letter in align_y if letter is not '-'])
print align_y


# In[14]:

res = score,align_x,align_y
with open('/Users/richard/Downloads/output','w') as file:
    for item in res:
        file.write(str(item)+"\n")


# In[ ]:



