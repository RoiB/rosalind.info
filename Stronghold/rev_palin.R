require(seqinr)
# layout original sequence and the reverse complement 
seq = s2c(read.fasta("~/Downloads/rosalind_revp (1).txt",seqonly=T)[[1]])
#rev_comp_seq = rev(comp(seq,forceToLower = F))

# build sliding window with size from 2 to length of the whole sequence
m = length(seq)
file.remove("~/Downloads/output")
for (len in 4:m) {
  for (i in 1:(m-(len-1))) {
    sub_seq = seq[i:(i+len-1)]; comp_seq = comp(sub_seq,forceToLower = F);
    #print(sub_seq);print(rev(comp_seq));print(is_rev_palin(sub_seq,comp_seq))
    if (is_rev_palin(sub_seq,comp_seq)) {
      write(c(i,len),"output",append = T);
      print(c(i,len))
    }
  }
  if (len == 12) break;
}
# function to compare if the two char arrays are reverse complement
is_rev_palin = function(seq1,seq2) {
  #print(seq1);print(seq2);print("#################")
  paste(seq1,collapse = "") == paste(rev(seq2), collapse = "")
}