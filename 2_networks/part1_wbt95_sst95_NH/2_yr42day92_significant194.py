import numpy as np



##--##--##--##--     --##--##--##--##
i = 194
j = 194
random_syn = np.zeros(1000, dtype=float)   ## 重复 1000 次
for ktimes in range(1000):   ## 重复 1000 次
  PC1 = np.zeros(92*42, dtype=float)
  PC2 = np.zeros(92*42, dtype=float)
  PC1[0:i]=1.0    ## PC1 有 i 天发生事件
  PC2[0:j]=1.0    ## PC2 有 j 天发生事件
  np.random.shuffle(PC1)
  np.random.shuffle(PC2)
  random_syn[ktimes] = np.sum(np.array(PC1+PC2==2), axis=0)

matrix99 = np.nanpercentile(random_syn,99)
print("  99 th  =  ", matrix99)  ## 17
    

