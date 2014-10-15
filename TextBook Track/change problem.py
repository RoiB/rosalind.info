
# coding: utf-8

# In[26]:

money = 40
coins = [1,5,10,20,25,50]


# In[29]:

def change(money,coins):
    table = [float("inf") for _ in range(money+1)]
    table[0] = 0
    
    def OPT(i):
        if i > 0:
            return min([table[i - j]+1 for j in coins if i-j>=0])
        else:
            return 0
    
    for i in range(1,money+1): table[i] = OPT(i)
    
    return table[money]


# In[31]:

money = 16563
coins = [1,3,5,7,8,11,20]


# In[32]:

change(money,coins)


# In[ ]:



