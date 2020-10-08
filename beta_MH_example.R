betaMH <- function(n, a, b, start, proposal_sd) {
  xbeta <- list()
  x <- rep(NA, n)
  accepted <- rep(NA, n)
  
  ## initialise chain
  x[1] <- start
  accepted[1] <- 1
  
  for (i in 2:n) {
    
    y <- rnorm(1, mean = x[i-1], sd = proposal_sd)
    
    r <- (dbeta(y, a, b)*dnorm(x[i-1], mean = y, sd = proposal_sd))/
      (dbeta(x[i-1], a, b)*dnorm(y, mean = x[i-1], sd = proposal_sd))
    
    aprob <- min(1, r)
    
    if (runif(1) < aprob) {
      x[i] <- y
      accepted[i] <- 1
    } else {
      x[i] <- x[i-1]
      accepted[i] <- 0
    }
  }
  list(x=x, accepted=accepted)
}

nmcmc <- 5000
burnin <- 500

xbeta <- betaMH(nmcmc + burnin, 2, 4, start = 0.5, proposal_sd = 0.5)
mean(xbeta$accepted)

plot_idxs <- (burnin + 1):(burnin + nmcmc)
plot(plot_idxs, xbeta$x[plot_idxs], type = "l")
plot(1:(nmcmc + burnin), xbeta$x, type = "l")
xbeta[["x"]]
hist(xbeta$x[plot_idxs])
