
# coding: utf-8

# In[4]:

def make_BLOSUM62():
    mat = {}
    with open("BLOSUM62.txt") as file:
        letters = file.readline().split()
        mat = {letter:{letter:None for letter in letters} for letter in letters}
        for _ in range(len(letters)):
            alist = file.readline().split()
            letter,scores = alist[0],alist[1:]
            for i in range(len(letters)):
                mat[letter][letters[i]] = int(scores[i])
        return mat


# In[6]:

BLOSUM62 = make_BLOSUM62()

