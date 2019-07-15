function normality_tests(data1,data2)
    % tests if each group has a normal distribution (kstest) then tests if
    % the variance is similar across each group (vartestn)
    
    % a signficant ks would indicate the data is not normal
    % a significant bartlett would indicate the data does not have similar
    % variance
    [~,p] = kstest(data1);
    disp([9,'ks-test, g1,p=',num2str(p)]);
    [~,p] = kstest(data2);
    disp([9,'ks-test, g2,p=',num2str(p)]);
    
    group = ones(length(data1)+length(data2),1);
    group(1:length(data1)) = group(1:length(data1))+1;
    
    try
        group_data = [data1,data2];
    catch
        group_data = [data1;data2];
    end
    
    if size(group_data,1) < size(group_data,2)
        p = vartestn(group_data',group,'display','off');
    else
        p = vartestn(group_data,group,'display','off');
    end
    disp([9,9,'Bartletts for equal variances,p=',num2str(p)]);
    %[~,p] = vartest2(data1,data2);
    %disp([9,'F-test for equal variances,p=',num2str(p)]);
end

