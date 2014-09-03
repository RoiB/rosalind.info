mass_table = read.table('~/Downloads/masstable.txt')
mass_table
which(mass_table$V1 == "P")
str(mass_table$V2)
#read while string from file
protein_seq = readLines("~/Downloads/test.txt",n = -1)
require(seqinr)
protein_chars = s2c(protein_seq)
tapply(protein_chars,1:nchar(protein_seq),)

get_mass = funtion(char) {mass_table$V2[which(mass_table$V1 == char)]}

