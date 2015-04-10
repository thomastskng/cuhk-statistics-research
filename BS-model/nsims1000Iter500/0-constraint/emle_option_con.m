
function result=emle_option_con(lambda, dl)

temp=1+dl*lambda;
nsims=size(dl,1);

%Q1n=sum(dl(:,1)./temp);
%Q2n=sum(dl(:,2)./temp);
%Q3n=sum(dl(:,3)./temp);
Qn=sum(bsxfun(@rdivide,dl,temp));
result=Qn/nsims;
%result=[Q1n, Q2n, Q3n]/nsims;


