function result=dY3dmu(t, R, mu, sigma, delta, moneyness, optionprice1)

r = 0.02;
result=exp(0.5*delta*(-mu+r)+(-mu+r)*R/sigma^2-delta*(r^2-mu^2)/(2*sigma^2)).*(-delta/2-R/sigma^2 + (mu*delta/sigma^2)).*(exp(-delta*r)*max(exp(R)-moneyness,0)-optionprice1);

