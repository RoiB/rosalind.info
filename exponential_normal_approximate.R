options(digits = 10)
qnorm(.94)

for (n in 1:50) {print( (40-1000/n)/(sqrt(500/n) *2) - qnorm(.94) )}

sapply(1:50,function(n) ((40-1000/n)/(sqrt(500/n) * 2) - qnorm(.94)))

for(n in 1:50) { print (40-1000/n)/(sqrt(500/n) * 2) - qnorm(.94) }


n = 35
(40-1000/n)/(sqrt(500/n) * 2) - qnorm(.94)
