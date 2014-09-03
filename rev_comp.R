require(seqinr)
data = read.fasta("rosalind_rvco (2).txt",seqonly = T)
seqs = lapply(data,s2c)
seqs

res = sapply(lapply(seqs,rev),comp,forceToLower = F)

res
sum(seqs %in% res)
