function CI = ConfInt(x)
%Calculates confidence intervals for a single vector of data


SEM = std(x)/sqrt(length(x));               % Standard Error
ts = tinv([0.025  0.975],length(x)-1);      % T-Score
CI = mean(x) + ts*SEM;                      % Confidence Interval

end

