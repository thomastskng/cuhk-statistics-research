function result=dY1dsigma(t,mu, sigma, delta)
result=exp(-delta*sigma^2*t^2/2)*(delta*sigma*t)*(-sin((mu-sigma^2/2)*delta*t)+t*cos((mu-sigma^2/2)*delta*t));