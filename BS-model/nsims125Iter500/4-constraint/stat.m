function stat (moneyness,seed)

nsims=125;
logS = zeros(nsims,1);
R = zeros(nsims, 1);
noofiter=50;
paratotal=zeros(noofiter, 4);
logS(1,1)=log(100);

sigma = 0.30;
mu = 0.05;
delta = 1/52;
r = 0.02;
options = optimset('Display','iter');%;,'TolX',1e-6, 'TolFun', 1e-6);
moneyness=[moneyness, moneyness + 0.01, moneyness + 0.02, moneyness + 0.03];

    checkmat_inner = zeros(noofiter,6+ (size(moneyness,2)*3));
    rng(seed)
    for k=1:noofiter,
        % generate 1001 logS geometric brownian motion samples
        for j=1:nsims
            logS(j+1,1) = logS(j,1) + (mu-sigma^2/2)*delta + randn(1,1)*sqrt(delta)*sigma;
        end
  
        R=logS(2:nsims+1,1)-logS(1:(nsims),1);
        Sn=exp(logS(nsims+1));  
        repR=repmat(R,1,size(moneyness,2));
        test = exp(repR)-repmat(moneyness,size(R,1),1);
        fail = sum(test<=0,1) / nsims;
        
        success = false;
        while success == false 
            try
                init = [(mean(R)/delta+var(R)/2),std(R)*sqrt(1/delta)];
                %init = [mu,sigma].*unifrnd(0.95,1.15,1,2);
                csvwrite('lastFailR.csv',R);
                csvwrite('lastFailSn.csv',Sn);
                [para,fval, exitflag] = fminsearch('emle_option_Lagranian', double(init), options, nsims, Sn, R, moneyness);
                [dl,eg3_final] = emle_optionscore(para, nsims, 1, Sn, R, moneyness);
                checkmat_inner(k,:) = [moneyness,fail,eg3_final,para,init,exitflag,seed];
                csvwrite('results.csv',checkmat_inner);
                csvwrite(strcat('exitflag',num2str(exitflag),'.csv'),R);
                exitflag
                %if exitflag == 1 | exitflag  == 2
                %    success = true
                %end
                success = true;
                fprintf('exitflag\n')              
            catch err
                init = [mu,sigma] * unifrnd(0.95,1.05)
                %init = [mean(exp(R)),std(R)]
                fprintf('RESTART\n')
            end
        end
    end
    
    csvwrite('results.csv',checkmat_inner);
    %means = mean(checkmat_inner,1);
    %biasedness = mean(checkmat_inner(:,3:4) - repmat([mu,sigma],noofiter,1),1);
    %se_mu = std(checkmat_inner(:,3));
    %se_sigma = std(checkmat_inner(:,4));
    %no_converged = sum(checkmat_inner(:,5)==1 | checkmat_inner(:,5)== 2)/noofiter;
    %checkmat_outer(1,:) = [moneyness,means(1:2), biasedness,se_mu, se_sigma, no_converged];
    %csvwrite('results.csv',checkmat_outer);
    %checkmat_outer
end