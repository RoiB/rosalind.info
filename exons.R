# read fasta file as sequence only
(dna_seq = read.fasta("~/Downloads/test.txt",seqonly = T))
require(seqinr)
#unwind dna_seq to chars
seq_char_list = sapply(dna_seq, s2c)
## give it a try for complement
seq_char_list = lapply(seq_char_list, comp)
lapply(seq_char_list,toupper)

seq_char_list
#translate seq_char_list to amino acid
amino_acids = list()
for (i in 0:2) {
  for (j in 1:3){
    amino_acids[length(amino_acids)+1] = sapply(seq_char_list, translate, frame = i)[j]
  }
  for (j in 1:3) {
    amino_acids[length(amino_acids)+1] = sapply(seq_char_list, translate, frame = i, sens = "R")[j]
  }
}
amino_acids
proteins = sapply(amino_acids,paste,collapse = "")
proteins
require(stringr)
matching = str_match(proteins,"(M\\w*)\\*")[,2]
matching
max_length = max(nchar(matching))
max_id = which(max_length == nchar(matching))[1]
res = matching[max_id]
res
write(res,"output")

