% Calculate dl

function dl=emle_optionscore(para, nsims, t, Sn, R)

mu=para(1);
sigma=para(2);
delta=1/52;
r=0.02;

nummax=size(R,1);
dl=zeros(nummax,3);
moneyness=0.99;

% Calculate the Black-Scholes call option price
%K1=Sn*moneyness;
d11=(log(1/moneyness)+(r+0.30^2/2)*delta)/(0.30*sqrt(delta));
d21=d11-0.30*sqrt(delta);
optionprice1=(normcdf(d11)-moneyness*exp(-r*delta)*normcdf(d21))/Sn;


dl(:,1) = cos(t*R)-real(exp(delta*(i*t*(mu-sigma^2/2)-sigma^2*t^2/2)));
dl(:,2) = sin(t*R)-imag(exp(delta*(i*t*(mu-sigma^2/2)-sigma^2*t^2/2)));
dl(:,3) = (exp(-r*delta)*max(exp(R)-moneyness,0)-optionprice1).*exp(-(r-mu)/(2*sigma^2*delta)*(-2*R*delta+(r+mu)*delta^2-sigma^2*delta^2));
                                                          
end
