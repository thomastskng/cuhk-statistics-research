# Check relationship between moneyness and expected failure rate
nsims <- 1000
noofiter <- 200
logS <- rep(0,nsims+1)
logS[1] <- log(100)

sigma <- 0.2943
mu <- 0.1289
delta <- 1/52
r <- 0.02

m <- seq(0.75,1,0.01)
tt <- matrix(NA, length(m), 3)
                                        #R1 <- read.csv("/Users/thomastskng/cuhk-statistics-research/EL_Opt_Matlab/R.csv", header = F)
R1 <- read.csv("cuhk-statistics-research\EL_Opt_Matlab\R.csv", header = F)

R1 <- as.vector(R1[,1])

for(i in 1:length(m)){
	ff <- 0
	g3 <- 0
	g31 <- 0
	# calculate BS call option price
	d11 <- (log(1/m[i]) + (r+sigma^2/2)*delta) / (sigma*sqrt(delta))
	d21 <- d11 - sigma*sqrt(delta)
	optionprice1 = pnorm(d11) - m[i]*exp(-r*delta)*pnorm(d21)

	for(k in 1:noofiter){
		fail <- 0
		for(j in 1:nsims){
			logS[j+1] =logS[j] + (mu-0.5*sigma^2)*delta + rnorm(1,0,1)*sqrt(delta)*sigma
		}

		R <- logS[2:length(logS)] - logS[1:nsims]
		ans <- exp(R) - m[i]
		ff <- sum(ans<=0) + ff

		# compute option constraint
		g3 <- mean((exp(-r*delta)*pmax(exp(R) - m[i],0)-optionprice1)*exp((-(r-mu)/(2*(sigma^2)*delta))*(-2*R*delta + (r+mu)*(delta^2) -(sigma^2)*(delta^2) ))) + g3

		#g31 <- mean(exp(-r*delta)*max(exp(R) - m[i],0) * exp(((r-mu)*R/sigma^2) - (r^2-mu^2)*delta/(2*sigma^2) + (r-mu)*delta/2) - optionprice1) + g31
		#dqdp <- exp((-delta/(2*sigma*2) ) * (mu^2 - mu*sigma^2 - r^2 + r*sigma^2))
		#g31 <- mean(exp(-r*delta)*pmax(exp(R) - m[i],0)*dqdp - optionprice1) + g31
	}

	tt[i,] <- c(m[i],ff/(noofiter*nsims), g3/(noofiter) #, g31/noofiter
				)
}
colnames(tt) <- c("m","fail(%)","E(g3)")
tt

# verify the expected value of option CONSTRAINT = 0
