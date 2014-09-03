require(gtools)

n = 2

all_perms = permutations(n,n)


split_cases = split(all_perms,1:factorial(n))


sign_matrix = permutations(2,n,v = c(-1,1), repeats.allowed = T)


raw_res = lapply(split_cases,function(x) (matrix(rep(x,2^n), ncol = n,byrow = T) * sign_matrix))


res = do.call(rbind,raw_res)

write(factorial(n)*2^n,"output")
write.table(res,"~/Downloads/output",row.names=F,col.names = F,append = T)

print(res)




