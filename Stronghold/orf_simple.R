require(seqinr)
seq = readLines("rosalind_orfr.txt")
seq_chars = s2c(seq) #nice function transfering string to characters


# build up translations from 6 reading frame 
trans = c()
for(i in 0:2) {
  trans= c(trans,paste(translate(seq_chars,frame = i),collapse=""))
  trans= c(trans,paste(translate(seq_chars,frame = i,sens = "R"),collapse=""))
}


require(stringr) #utilize str_match function to produce matches
matches = str_match(trans,"(M\\w*?)\\*")[,2] # fisrt round matching without overlapping

#########################
# take overlapping into consideration
while(sum(grepl("M",trans))) {
  trans = sub("M","-",trans)
  matches = c(matches, str_match(trans,"(M\\w*?)\\*")[,2])
}

id = which(max(nchar(matches)) == nchar(matches))

write(matches[id],"output")
