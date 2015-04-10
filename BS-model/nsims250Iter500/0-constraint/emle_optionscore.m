% Calculate dl

function [dl,eg3]=emle_optionscore(para, nsims, t, Sn, R, moneyness)

mu=para(1);
sigma=para(2);
delta=1/52;
r=0.02;

nummax=size(R,1);
dl=zeros(nummax,2);

% Calculate the Black-Scholes call option price
%K1=Sn*moneyness;
d11=(log(1./moneyness)+(r+0.30^2/2)*delta)/(0.30*sqrt(delta));
d21=d11-0.30*sqrt(delta);
optionprice1=normcdf(d11)-moneyness.*(exp(-r*delta)*normcdf(d21));

dl(:,1) = cos(t*R)-real(exp(delta*(i*t*(mu-sigma^2/2)-sigma^2*t^2/2)));
dl(:,2) = sin(t*R)-imag(exp(delta*(i*t*(mu-sigma^2/2)-sigma^2*t^2/2)));

repR=repmat(R,1,size(moneyness,2));
m=repmat(moneyness,size(repR,1),1);
optionprice = repmat(optionprice1, size(repR,1),1);
%dl(:,3:(2+size(moneyness,2))) = (exp(-r*delta)*bsxfun(@max,exp(repR)-m,0).*exp((2*repR*(r-mu) + delta*(mu^2  - r^2 + sigma^2 *(r-mu)))/(2*sigma^2))-optionprice);

%Z = exprnd(1,nsims,1);
%dl = [dl(:,1).*Z];
%dl = dl(:,1);
eg = mean(dl,1);
%eg3 = eg(3:(2+size(moneyness,2)));
eg3=0;
end
