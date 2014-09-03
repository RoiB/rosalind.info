options(digits = 5)


p = exp(-9)-exp(-12)
p

pexp(12) - pexp(9)
lambda  = 25000 * p
lambda


res = ppois(1:25000,lambda = lambda)
plot(1:25000,res)
head(res,100)

data = data.frame(num = 1:11, prob = res[1:11])
remove("problem4_output.txt")

write.table(data,"problem4_output.txt")

for (i in 1:11) { write(data[i,],"problem4_output.txt",append = T) }

?write.table
