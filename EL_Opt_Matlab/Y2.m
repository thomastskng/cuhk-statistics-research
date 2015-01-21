function result=Y2(t, R, mu, sigma, delta)
result=sin(t*R)-exp(-delta*sigma^2*t^2)*sin(t*(mu-sigma^2/2)*delta);
