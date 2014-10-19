
# coding: utf-8

# In[1]:

def make_PAM250():
    mat = {}
    with open("PAM250.txt") as file:
        letters = file.readline().split()
        mat = {letter:{letter:None for letter in letters} for letter in letters}
        for _ in range(len(letters)):
            alist = file.readline().split()
            letter,scores = alist[0],alist[1:]
            for i in range(len(letters)):
                mat[letter][letters[i]] = int(scores[i])
        return mat


# In[2]:

PAM250 = make_PAM250()

