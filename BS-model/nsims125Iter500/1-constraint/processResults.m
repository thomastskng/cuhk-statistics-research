function processResults(mu,sigma)
    cd('/Users/thomastskng/cuhk-statistics-research/BS-model/nsims125Iter500/1-constraint');
    num=1
    moneyness=repmat({'moneyness'},1,num);
    fail=repmat({'m > exp(R)'},1,num);
    expected_val= repmat({'E(g)'},1,num);
    
    header = {moneyness{:}, fail{:}, expected_val{:},'bias(mu)', 'bias(sigma)', 'se(mu)', 'se(sigma)', 'E(mu_init)', 'E(sigma_init)','flags(1 and 0)','flag (1)', 'flag(0)'};
    fid = fopen('all-in-one-by-time.csv','w')
    fprintf(fid,'%s,',header{1,:});
    fprintf(fid, '\n') ;
    fclose(fid);
    file_arr={ 'nsims125Iter500constraints1.csv'};
    scol = num*3;
 for i = 1:numel(file_arr)
    mat = csvread(file_arr{i});
    moneyness = unique(mat(:,1));
    checkmat_outer = zeros(size(moneyness,1),size(header,2));
    for m =1:size(moneyness,1)
        selrows = mat(mat(:,1)==moneyness(m), :);
        %selrows = selrows(selrows(:,6) == 2 | selrows(:,6) == 1,:)
        means = mean(selrows,1)
        biasedness = mean(selrows(:,(scol+1):(scol+2)) - repmat([mu,sigma],size(selrows,1),1),1);
        se_mu = std(selrows(:,(scol+1)));
        se_sigma = std(selrows(:,(scol+2)));
        no_converged = sum(selrows(:,(scol+5))== 1)/size(selrows,1);
        no1_converged = sum(selrows(:,(scol+5))== 1)/size(selrows,1);
        no2_converged = sum(not(selrows(:,(scol+5))== 1))/size(selrows,1);
        checkmat_outer(m,:) = [means(1:scol), biasedness,se_mu, se_sigma, means((scol+3):(scol+4)),no_converged,no1_converged,no2_converged];
    end
    checkmat_outer
    fid = fopen('all-in-one-by-time.csv','a+')
    fprintf(fid,'%s',file_arr{i});
    fprintf(fid, '\n') ;
    fclose(fid);
    dlmwrite('all-in-one-by-time.csv',checkmat_outer,'delimiter',',','-append');
 end
end