
# coding: utf-8

# In[ ]:

'''
implementation of sequence alignment
Author: Richard Ding
'''


# In[2]:

def build_scoring_matrix(alphabet, diag_score, off_diag_score, dash_score):
    '''
    Takes as input a set of characters alphabet and three scores diag_score, off_diag_score, and dash_score. 
    The function returns a dictionary of dictionaries whose entries are indexed by pairs of characters in alphabet plus '-'. 
    The score for any entry indexed by one or more dashes is dash_score. 
    The score for the remaining diagonal entries is diag_score. 
    Finally, the score for the remaining off-diagonal entries is off_diag_score.
    '''
    alphabet = set(alphabet)
    alphabet.add('-')
    mat = {elem:{letter: -1 for letter in alphabet} for elem in alphabet}
    
    for elem in mat:
        for letter in mat[elem]:
            
            if letter == elem: 
                if letter != '-':
                    mat[elem][letter] = diag_score
                else:
                    mat[elem][letter] = dash_score
            elif letter == '-' or elem == '-':
                mat[elem][letter] = dash_score
            else:
                mat[elem][letter] = off_diag_score
    return mat


# In[3]:

def compute_alignment_matrix(seq_x, seq_y, scoring_matrix, global_flag):
    '''
    The function computes and returns the alignment matrix for seq_x and seq_y as described in the Homework.
    '''
    len_x = len(seq_x)
    len_y = len(seq_y)
#     print len_x,len_y
    score = [[-1 for _ in range(len_y+1)] for _ in range(len_x+1)]
#     print len(score),len(score[0])
#     print score
    score[0][0] = 0
    
    if global_flag:
        for id_i in range(1,len_x+1):
            score[id_i][0] =  score[id_i-1][0] + scoring_matrix[seq_x[id_i-2]]['-']
            
        for id_j in range(1,len_y+1):
#             print id_j
            score[0][id_j] =  score[0][id_j-1] + scoring_matrix['-'][seq_y[id_j-2]]
        
        for id_i in range(1,len_x+1):
            for id_j in range(1,len_y+1):
                score[id_i][id_j] =  max(score[id_i-1][id_j]+scoring_matrix[seq_x[id_i-1]]['-'],                                         score[id_i][id_j-1]+scoring_matrix['-'][seq_y[id_j-1]],                                         score[id_i-1][id_j-1] + scoring_matrix[seq_x[id_i-1]][seq_y[id_j-1]])
    else:
        for id_i in range(1,len_x+1):
            score[id_i][0] =  0
            
        for id_j in range(1,len_y+1):
#             print id_j
            score[0][id_j] =  0
        
        for id_i in range(1,len_x+1):
            for id_j in range(1,len_y+1):
                score[id_i][id_j] =  max(score[id_i-1][id_j]+scoring_matrix[seq_x[id_i-1]]['-'],                                         score[id_i][id_j-1]+scoring_matrix['-'][seq_y[id_j-1]],                                         score[id_i-1][id_j-1] + scoring_matrix[seq_x[id_i-1]][seq_y[id_j-1]],0)
    return score
            


# In[4]:

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


# In[14]:

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


# In[15]:

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
#                 print max_val
                max_idx = (id_i,id_j)    
    return max_idx


# In[16]:

# ## Test
# score_mat = build_scoring_matrix(["A","T","G","C","-"], 10, 4, -6)
# #for row in score_mat: print row,score_mat[row]
# dp_table_global = compute_alignment_matrix('AA', 'TAAT', score_mat, True)
# dp_table_local = compute_alignment_matrix('AA', 'TAAT', score_mat, False)
# for row in dp_table_global: print row
# print 
# for row in dp_table_local: print row
# print 
# print compute_global_alignment('AA', 'TAAT', score_mat, dp_table_global)
# print compute_local_alignment('AA', 'TAAT', score_mat, dp_table_local)

