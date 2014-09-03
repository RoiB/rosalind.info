require(seqinr)
seq = read.fasta("~/Downloads/rosalind_kmer.txt",seqonly = T)[[1]]
seq_chars = s2c(seq)
##########################
table = c("A","C","G","T")
require(gtools)
perms = permutations(4,4,repeats.allowed = T)

four_mers_idx = split(perms,1:256)
four_mers= sapply(four_mers_idx, function(x) (paste(table[x],collapse = "")))
four_mers_count = rep(0,256)
names(four_mers_count) = four_mers
four_mers_count

for (i in 1:(length(seq_chars)-3)) {
  idx = paste(seq_chars[i:(i+3)],collapse=""); 
  four_mers_count[idx] = four_mers_count[idx] +1
}
write(four_mers_count,"~/Downloads/output",ncolumns = 256)

