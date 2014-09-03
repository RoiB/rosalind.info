ids = strsplit(readLines("rosalind_frmt.txt")," ")[[1]]
require(seqinr)
choosebank("genbank")

############ testing ############

info = lapply(ids,function(id){param = paste("AC=",id,sep= "");query(param)})

min_len_idx = function(info) {
  get_len = function(alist) {return(attr(alist$req[[1]],"length"))}
  len_array = sapply(info, get_len)
  which(len_array == min(len_array))
}

ids[min_len_idx(info)]
