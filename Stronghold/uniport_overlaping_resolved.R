require(seqinr) # used for read fasta
ids = readLines("rosalind_mprt.txt",n = -1) #read Database access IDs
#ids = readLines("test.txt",n = -1) #read Database access IDs
main_url = "http://www.uniprot.org/uniprot/"
for (i in 1:length(ids)) {
  # make file name
  filename = paste(ids[i],".fasta",sep = "")
  #download fasta file to local
  download.file(paste(main_url, filename,sep = ""),
                destfile = filename, quiet = T)
  #read sequence from fasta file
  (seq = read.fasta(filename,seqonly = T)[[1]])
  #use regex, output indices to file "output"
  #there is no overlaping solution, need to write loop, so sad
  (match_idx = gregexpr("[N][^P][ST][^P]",seq)[[1]])
  #loop four times to make sure no, still have problem
  for (j in 1:4) {
    seq = trim_seq(seq,match_idx)
    match_idx = sort(unique(c(match_idx,gregexpr("[N][^P][ST][^P]",seq)[[1]])))
  }
  if (!is.na(match_idx[2])) {
    write(ids[i], "output", append = T); 
    write(match_idx[2:length(match_idx)],"output",append = T, sep = " ",ncolumns = 100)
    print(ids[i])
    print(match_idx[2:length(match_idx)])
  }  
}
trim_seq = function(seq,match_idx) {
  for (i in match_idx) {
    seq = paste(substr(seq,1,i-1),"-",substr(seq,i+1,nchar(seq)),sep = "")
  }
  return(seq)
}