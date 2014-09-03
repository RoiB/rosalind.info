data = readLines("~/Downloads/rosalind_lexf.txt")
require(seqinr)
#is_sorted = function(seq) {prod(sort(seq) == seq)} # function to check if char array is sorted
chars = s2c(gsub(" ","",data[1]))
k = as.numeric(data[2])
n = length(chars)

require(gtools)
perms = permutations(n,k,repeats.allowed=T) # generate all k-mer permutation
mapping = function(n) {return(chars[n])} # function to map index to char
char_groups = apply(perms,1,function(num_arr) (sapply(num_arr,mapping))) # lambda dealing with index vector
res = apply(char_groups,2,paste,collapse = "") # glue groups together
res
write(res,"output",sep = "\n")
