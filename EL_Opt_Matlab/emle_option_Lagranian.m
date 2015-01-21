%Calculate the log-empirical likelihood function
function result=emle_option_Lagranian(para, nsims, Sn, R)
t1=[-0.50:0.1:-0.01,-0.01:0.1:0.50]*10;

tempsum=0;

for j = 1:size(t1,2)
    
%calculate 'g', where g is the cond char function constraint
   dl = emle_optionscore(para, nsims, t1(j), Sn, R);
   dl = dl*inv(chol(dl'*dl));

%Find the lagrange multiplier

    mlambda=inv(dl'*dl)*(ones(1, nsims)*dl)';
    
    options=optimset('fsolve');
    options=optimset(options, 'Display', 'off');
    options=optimset(options, 'diagnostics', 'off');
    options=optimset(options, 'LargeScale', 'off');
        
    %solve for lambda
    [lambda, fvals, exitflag]=fsolve('emle_option_con', mlambda, options, dl);
       
    tempsum=sum(log(1+ dl*lambda))+tempsum;
    
end
lambda
result=tempsum/size(t1,2);


