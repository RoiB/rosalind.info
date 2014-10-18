
# coding: utf-8

# In[14]:

import itertools


# In[15]:

def readfile(filename):
    with open(filename) as file:
        return file.readline().split(), int(file.readline())


# In[20]:

symbols, n = readfile("/Users/richard/Downloads/rosalind_lexv.txt")
numbers = tuple(range(len(symbols)))


# In[21]:

def indices(symbols,n): return itertools.product(numbers,repeat = n)


# In[22]:

res = []
for i in range(n):
    res.extend(indices(numbers,i+1))


# In[23]:

with open("/Users/richard/Downloads/output",'w') as file:
    for tup in sorted(res):
        temp = ""
        for j in tup:
            temp+=symbols[j]
        file.write(temp+"\n")

