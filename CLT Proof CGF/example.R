set.seed(1928)
n <- 200
M <- 10000
lambda <- 1/3
results <- matrix(0, nrow = M, ncol = n)
for(m in 1:M){
  results[m, ] <- rexp(n, lambda)
}
# symmetry of the sample mean arrives with n smaller (30ish)
hist( rowMeans(results) )
abline(v = 1/lambda, col = 'red')

values <- apply(results, MARGIN = 2, 
                function(x){quantile(x, 0.25)})
# quantile values require larger n (e.g. 200) 
hist(values)


