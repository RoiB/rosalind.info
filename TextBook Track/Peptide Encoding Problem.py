
# coding: utf-8

# In[58]:

import Bio as bp


# In[59]:

record = list(bp.SeqIO.parse("/Users/richard/Downloads/rosalind_2b (1).txt","fasta"))[0]


# In[60]:

seq = str(record.seq)


# In[61]:

peptide = "NIIPNECA"


# In[62]:

m = len(peptide)*3


# In[63]:

def sub_str(seq,m,peptide):
    res = []
    for i in range(len(seq)-m+1):
        if ( str(bp.Seq.Seq(seq[i:i+m]).translate()) == peptide or              str(bp.Seq.Seq(seq[i:i+m]).reverse_complement().translate())  == peptide):
            res.append(seq[i:i+m])
    return res


# In[64]:

res = sub_str(seq,m,peptide)


# In[65]:

print res


# In[66]:

with open("/Users/richard/Downloads/output","w") as file:
    for item in res: file.write(item+"\n")

