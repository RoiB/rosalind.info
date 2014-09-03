dna_seq = read.fasta("~/Downloads/rosalind_splc (1).txt",seqonly = T)
exon = dna_seq[[1]]
for (i in 2:length(dna_seq)) {exon = gsub(dna_seq[[i]],"",exon)}

res = paste(translate(s2c(exon)),collapse = "")
final_res = str_match(res,"(M\\w*?)\\*")[,2]
write(final_res,"output")
