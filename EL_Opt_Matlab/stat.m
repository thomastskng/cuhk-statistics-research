function stat (dummy)

nsims=1000;
logS = zeros(nsims,1);
R = zeros(nsims, 1);
noofiter=100;
paratotal=zeros(noofiter, 4);
logS(1,1)=log(100);


sigma = 0.30;
mu = 0.095;
delta = 1/52;
r = 0.02;
t1=[-0.50:0.1:-0.01,-0.01:0.1:0.50]*10;
options = optimset('Display','iter');%,'TolX',1e-4, 'TolFun', 1e-8);
mm = [0.75:0.01:1];
checkmat_outer = zeros(size(mm,2),8);

for m = 1:size(mm,2),
    checkmat_inner = zeros(noofiter,5);
    rng(1)
    for k=1:noofiter,
        % generate 1001 logS geometric brownian motion samples
        for j=1:nsims
            logS(j+1,1) = logS(j,1) + (mu-0.5*sigma^2)*delta + randn(1,1)*sqrt(delta)*sigma;
        end
  

        R=logS(2:nsims+1,1)-logS(1:(nsims),1);
        Sn=exp(logS(nsims+1));   

        test = exp(R)-mm(m);
        fail = sum(test<=0) / 1000;
  
        %init=[mu,0.7];
        %init = [0.1,0.75]     % converged, negative sigma
        %init = [0.095, 0.2]   % saddlepoint ?
        %init = [0.5,0.25]     % converged, diff mu, same sigma
        %init = [0.55,0.80]     % converged
        init = [mu,0.3]       % converged, positive sigma, same as above
        success = false;
        while success == false 
            try
                [para,fval, exitflag] = fminunc('emle_option_Lagranian', init, options, nsims, Sn, R);
                [dl,eg3_final] = emle_optionscore(para, nsims, 1, Sn, R);
                checkmat_inner(k,:) = [fail,eg3_final,para,exitflag];
                success = true;
            catch err
                for j=1:nsims
                    logS(j+1,1) = logS(j,1) + (mu-0.5*sigma^2)*delta + randn(1,1)*sqrt(delta)*sigma;
                end
 
                R=logS(2:nsims+1,1)-logS(1:(nsims),1);
                Sn=exp(logS(nsims+1));   
                
                test = exp(R)-mm(m);
                fail = sum(test<=0) / 1000;
                fprintf('RESTART\n')
            end
        end
    end
    means = mean(checkmat_inner,1);
    biasedness = mean(checkmat_inner(:,3:4) - repmat([mu,sigma],noofiter,1),1);
    se_mu = std(checkmat_inner(:,3));
    se_sigma = std(checkmat_inner(:,4));
    no_converged = sum(checkmat_inner(:,5)==1)/noofiter;
    checkmat_outer(m,:) = [mm(m),means(1:2), biasedness,se_mu, se_sigma, no_converged];
    csvwrite('results.csv',checkmat_outer);
end

checkmat_outer

end