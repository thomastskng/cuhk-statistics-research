function result=dY3dsigma(t, R, mu, sigma, delta, moneyness, optionprice1)

r=0.02;
result=(exp(-r*delta)*max(exp(R)-moneyness,0)-optionprice1).*exp(-(r-mu)/(2*sigma^2*delta)*(-2*R*delta+(r+mu)*delta^2-sigma^2*delta^2)).*((r-mu)*(-2*delta*R+(mu+r)*delta^2)/(sigma^3*delta));