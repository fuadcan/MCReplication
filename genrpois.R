n <- 10; k <- 2
rpois(k,n/k) 6,3

n <- 10; k <- 3
rpois(k,n/k) 4,2,4

n <- 20; k <- 4
rpois(k,n/k) 3,2,5,5

n <- 20; k <- 5
rpois(k,n/k) 4,4,4,3,2

n <- 30; k <- 5
rpois(k,n/k) 6,6,4,4,3

n <- 30; k <- 6
rpois(k,n/k) 8,7,5,4,3,3

n <- 40; k <- 7
rpois(k,n/k) 7,7,7,6,6,4,2

n <- 40; k <- 6,6,4,4,3,3,3,2
rpois(k,n/k)
a <- 1

genpois <- function(n,k){
  out <- rpois(k,n/k)
  while(sum(out) > n | sum(out) < n*2/3){
    out <- rpois(k,n/k)
  }
  out
}

