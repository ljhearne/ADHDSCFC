function data = Supplement_stats(path,K)
%Takes a path to data txt file and returns as numbers

content = fileread(path);

stat_label = {'Deg t-test, ',...
    'Weighted deg t-test, ',...
    ['Hub STRENGTH t-test- k = ',num2str(K),', conn = 1 '],...
    ['Hub STRENGTH t-test- k = ',num2str(K),', conn = 2 '],...
    ['Hub STRENGTH t-test- k = ',num2str(K),', conn = 3 '],...
    'Connectome-wide t-test (NORMAL), ',...
    'Hub NORMAL t-test, ',...
    'Feeder NORMAL t-test, ',...
    'Periphery NORMAL t-test, '};

for j = 1:length(stat_label)
    %extract pvals from text
    substr = [stat_label{j},'pval = '];
    loc = strfind(content,substr);
    d = sscanf(content(loc+numel(substr):end), '%f', 1);
    data(j) = round(d,3);
end

substr = 'Correlation with Feeder SC-FC r = ';
loc = strfind(content,substr);
d = sscanf(content(loc+numel(substr):end), '%f', 1);
data(j+1) = round(d,2);
end

