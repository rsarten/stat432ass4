library(mixtools)
library(magrittr)

dirichletMH <- function(n, start, alpha, prop_width) {
  
  x_values <- vector("list", n)
  proposals <- vector("list", n)
  accepted <- rep(0, n)
  
  ## initialise chain
  x_values[[1]] <- start
  accepted[1] <- 1
  
  proposal_value <- function(prev_x, pos, prop_width)
    runif(1, min = prev_x[pos] - prop_width, max = prev_x[pos] + prop_width)
  
  for (i in 2:n) {
    prev_x <- x_values[[i-1]]
    x1 <- proposal_value(prev_x, 1, prop_width)
    x2 <- proposal_value(prev_x, 2, prop_width)
    x3 <- 1 - x1 - x2
    y <- c(x1, x2, x3)
    proposals[[i]] <- c(x1, x2, x3)
    
    r <- ddirichlet(x = y, alpha)/ddirichlet(x = prev_x, alpha)
    
    aprob <- min(1, r)
    if (runif(1) < aprob) {
      x_values[[i]] <- y
      accepted[i] <- 1
    } else {
      x_values[[i]] <- prev_x
      accepted[i] <- 0
    }
  }
  list(x_values = x_values, proposals = proposals, accepted = accepted)
}

burnin <- 1000
nmcmc <- 4000
lastk <- (nmcmc+1):(burnin+nmcmc)

values <- dirichletMH(burnin + nmcmc, c(0.1, 0.1, 0.8), c(1, 1, 1), 0.5)


x <- do.call(rbind, values$x_values) %>% as.data.frame()
names(x) <- paste0("x", 1:3)

pr <- do.call(rbind, values$proposals) %>% as.data.frame()
names(pr) <- paste0("x", 1:3)

plot(pr$x1, pr$x2, type = "p")

mean(values$accepted)

plot(1:burnin, x$x1[1:burnin], type = "l")
plot(1:1000, x$x1[4001:5000], type = "l")
acf(x$x1)
hist(x$x1)

plot(1:burnin, x$x2[1:burnin], type = "l")
plot(1:1000, x$x2[4001:5000], type = "l")
acf(x$x2)
hist(x$x2)

rDirichlet <- function(ndraw, alpha) {
  gamdraw <- rgamma(ndraw, shape = alpha, rate = 1)
  gamdraw/sum(gamdraw)
}

xdir.gen <- do.call(rbind, lapply(1:1000, function(i) rDirichlet(3, 1)))
hist(xdir.gen)
