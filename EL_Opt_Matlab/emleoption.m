function emleoption (dummy)


nsims=1000;
logS = zeros(nsims,1);
R = zeros(nsims, 1);
noofiter=100;
paratotal=zeros(noofiter, 4);
logS(1,1)=log(100);
%
sigma = 0.30;
mu = 0.095;
delta = 1/52;
r = 0.02;
t1=[-0.50:0.1:-0.01,-0.01:0.1:0.50]*10;

options = optimset('Display','iter');%,'TolX',1e-4, 'TolFun', 1e-8);


for k=1:noofiter,
  
    % generate 1001 logS geometric brownian motion samples
  for j=1:nsims
     logS(j+1,1) = logS(j,1) + (mu-0.5*sigma^2)*delta + randn(1,1)*sqrt(delta)*sigma;
  end
  

  R=logS(2:nsims+1,1)-logS(1:(nsims),1);
  Sn=exp(logS(nsims+1));   

  %init=[mu,0.7];
  init = [0.2,0.6]
  [para, fval, exitflag] = fminunc('emle_option_Lagranian', init, options, nsims, Sn, R);
  

  paratotal(k,:)=[para, fval, exitflag];
  
  muest=para(1);
  sigmaest=para(2);
  nocon=3;
  Y1obs=zeros(size(t1,2), nsims);
  Y2obs=zeros(size(t1,2), nsims);
  Y3obs=zeros(size(t1,2), nsims);
  dY1dmuobs=zeros(size(t1,2), nsims);
  dY1dsigmaobs=zeros(size(t1,2), nsims);
  dY2dmuobs=zeros(size(t1,2), nsims);
  dY2dsigmaobs=zeros(size(t1,2), nsims);
  dY3dmuobs=zeros(size(t1,2), nsims);
  dY3dsigmaobs=zeros(size(t1,2), nsims);
  gammaobs=zeros(nocon,nocon,size(t1,2),size(t1,2));

  s11=zeros(nocon,nocon,size(t1,2));
  s12=zeros(nocon,2,size(t1,2));
  A=zeros(2,2);
  B=zeros(2,2);
  
  %Black-Scholes Formula for European call 
  moneyness = 0.99;
  d11=(log(1/moneyness)+(r+0.30^2/2)*delta)/(0.30*sqrt(delta));
  d21=d11-0.30*sqrt(delta);
  optionprice1=(normcdf(d11)-moneyness*exp(-r*delta)*normcdf(d21))/Sn;

    for j=1:size(t1,2)
     
       Y1obs(j,:)=Y1(t1(j), R, muest, sigmaest, delta);
       Y2obs(j,:)=Y2(t1(j), R, muest, sigmaest, delta);
       Y3obs(j,:)=Y3(t1(j), R, muest, sigmaest, delta, moneyness, optionprice1);
       dY1dmuobs(j,:)=dY1dmu(t1(j), muest, sigmaest, delta);
       dY1dsigmaobs(j,:)=dY1dsigma(t1(j), muest, sigmaest, delta);
       dY2dmuobs(j,:)=dY2dmu(t1(j), muest, sigmaest, delta);
       dY2dsigmaobs(j,:)=dY2dsigma(t1(j), muest, sigmaest, delta);
       dY3dmuobs(j,:)=dY3dmu(t1(j), R, muest, sigmaest, delta, moneyness, optionprice1);
       dY3dsigmaobs(j,:)=dY3dsigma(t1(j), R, muest, sigmaest, delta, moneyness, optionprice1);
      
       for n=1:nsims
           s11(:,:,j)=s11(:,:,j)+[Y1obs(j,n),Y2obs(j,n), Y3obs(j,n)]'*[Y1obs(j,n),Y2obs(j,n), Y3obs(j,n)];
           s12(:,:,j)=s12(:,:,j)+[dY1dmuobs(j,n),dY1dsigmaobs(j,n);dY2dmuobs(j,n),dY2dsigmaobs(j,n);dY3dmuobs(j,n),dY3dsigmaobs(j,n)];
       end
       %s11(:,:,j) = [Y1obs(j,:);Y2obs(j,:); Y3obs(j,:)] * [Y1obs(j,:);Y2obs(j,:); Y3obs(j,:)]';
       %s12(:,:,j) = [ones(1,nsims);zeros(1,nsims);zeros(1,nsims)]*[dY1dmuobs(j,:);dY1dsigmaobs(j,:)]' + [zeros(1,nsims);ones(1,nsims);zeros(1,nsims)]*[dY2dmuobs(j,:);dY2dsigmaobs(j,:)]' + [zeros(1,nsims);zeros(1,nsims);ones(1,nsims)]*[dY3dmuobs(j,:);dY3dsigmaobs(j,:)]'
       
       s11(:,:,j)=-s11(:,:,j)/nsims;
       s12(:,:,j)=s12(:,:,j)/nsims;   
       A=A+s12(:,:,j)'*inv(s11(:,:,j))*s12(:,:,j)/size(t1,2);
    end

    for j=1:size(t1,2)
        for l=1:size(t1,2)
            for n=1:nsims
                gammaobs(:,:,j,l)=gammaobs(:,:,j,l)+[Y1obs(j,n), Y2obs(j,n), Y3obs(j,n)]'*[Y1obs(l,n),Y2obs(l,n),Y3obs(l,n)];
            end
            %gammaobs(:,:,j,l) = [Y1obs(j,:), Y2obs(j,:), Y3obs(j,:)]*[Y1obs(l,:),Y2obs(l,:),Y3obs(l,:)]';
        end
    end

    gammaobs=gammaobs/nsims;

    for j=1:size(t1,2)
        for l=1:size(t1,2)
            B=B+s12(:,:,j)'*inv(s11(:,:,j))*gammaobs(:,:,j,l)*inv(s11(:,:,l))*s12(:,:,l);        
        end
    end  

    B=B/size(t1,2)^2;
    Z(:,:,k)=inv(A)*B*inv(A);
    toc
%
save para paratotal;
save Z Z
end
end
