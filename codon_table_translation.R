require(seqinr)
data = readLines("rosalind_ptra.txt")
seq_char = s2c(data[1]);seq_char
ref_seq = s2c(data[2]);ref_seq

#translate according to 23 genetic code variant(codon table)
code_specific_trans = function(num) {translate(seq_char,numcode=num)}
choices = lapply(1:23,code_specific_trans)

translations = sapply(choices,paste,collapse = "")
translations
str(translations)
res = which(translations == paste(data[2],"*",sep=""))
res

#### regular expression appears to be wrong, probably because the pattern is too long
grep(data[2],translations) #translations(from 23 codon table) that match the result
