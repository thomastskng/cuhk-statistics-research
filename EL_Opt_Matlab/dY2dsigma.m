function result=dY2dsigma(t,mu, sigma, delta)
result=exp(-delta*sigma^2*t^2/2)*(delta*sigma*t)*(cos(t*(mu-sigma^2/2)*delta)+t*sin((mu-sigma^2/2)*t*delta));