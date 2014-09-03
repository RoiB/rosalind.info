# protein as all abbrevs of amino acids
# codon_quant as number of condon representation of each amino acids
protein = strsplit("ARNDCQEGHILKMFPSTWYV","")[[1]]
codon_quant = as.numeric(strsplit("46222224236212464124","")[[1]]) 

#read file and store it in "amino acids" as a char array
filename = "rosalind_mrna (2).txt"
fileinfo = file.info(filename)
a = readChar(filename,fileinfo$size)
amino_acids = strsplit(a,"")[[1]]
m = length(amino_acids)
if (amino_acids[m] == "\n") {amino_acids = amino_acids[1:m-1]}
#find index of num. of representation of each chr and save it as quant
find_index = function(x) {which(x == protein)}
quant_idx = sapply(amino_acids, find_index)
num_stop = 3
#(prod(codon_quant[quant_idx]) *num_stop) %% 1000000 #first try, fail, product is too large

### Second try, split numbers to groups, multiply each group seperately
### take the sum of mods then mods 1000000 again
seq = c(3,codon_quant[quant_idx]) #array of number of possible codon representation, plus stop codon
reduce = function(seq,num) {
  levels = round(seq_along(seq)/min(num,length(seq))); #print(levels) #make levels, each has 20 numbers
  chunks = split(seq,levels);#print(chunks)
  mods = sapply(chunks,prod) %% 1000000 #multiple numbers in each chunck and take modulo
}

seq1 = reduce(seq,20)
seq2 = reduce(seq1,2)
seq3 = reduce(seq2,2);seq3
seq4 = reduce(seq3,2);seq4
seq5 = reduce(seq4,2);seq5
seq6 = reduce(seq5,2);seq6
seq7 = reduce(seq6,2);seq7

(res = (prod(seq7)) %% 1000000)


