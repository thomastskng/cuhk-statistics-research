function result=dY1dmu(t, mu, sigma, delta)
result=exp(-sigma^2*t^2*delta/2)*sin((mu-sigma^2/2)*delta*t)*(delta*t);