
# coding: utf-8

# In[9]:

def readfile(filename):
    with open(filename) as file:
        n = int(file.readline())
        chars = file.readline().rstrip() 
        GCs = map(float,file.readline().split())
        return n,chars,GCs


# In[10]:

def compute_expect(GC,chars,n_choice):
    prob_GC = GC/2
    prob_AT = (1-GC)/2
    prob_table = {'A':prob_AT,'T':prob_AT,'G':prob_GC,'C':prob_GC}
    prob = 1
    for char in chars: prob *= prob_table[char]
    return n_choice*prob


# In[15]:

n,chars,GCs = readfile("/Users/richard/Downloads/rosalind_eval.txt")
res = []
n_choice = n - len(chars)+1
for GC in GCs: res.append(compute_expect(GC,chars,n_choice))


# In[17]:

with open("/Users/richard/Downloads/output",'w') as file:
    for i in res: file.write(str(i)+" ")

