function [r] = calculate_corr_sc_fc(SC, FC, HBL, connectivity_type)
%% Returns the correlation coefficient between SC and FC
%
% ARGUMENTS:
%    SC  -- an N x N matrix -- structural connectivity 
%    FC  -- an N x N matrices -- functional connecitivty 
%    HBL -- an N x N x 3 array with the indices of the type of edge we should use for correlation
%        -- assumed order hub: 1, feeder: 2, periphery: 3
%    connectivity_type -- a string with the type of edge we want to use. 
%                         Note that the current conn_type index us harcoded and dependent on the order in 
%                         the dataset I got.
% OUTPUT:
%    r -- linear correlation coefficient
%
% AUTHOR:
%     Derived from Luke's FC-SC analysis at:
%     https://github.com/ljhearne/ADHDSCFC/blob/master/Scripts/StrucFunc_analysis.m
% 
%     PSL 2019, QIMR Berghofer
%
% USAGE:
%{
    
%}

% NOTE: This function returns the correlation coefficient between the 
% elements of the lower triangular parts of the SC and FC matrices. 
% Only elements for which SC is nonzero are taken into account. 

if nargin < 3
    HBL = [];
    connectivity_type = 'all';
end


   switch connectivity_type
    case {'hubs', 'hub'}
        idx = ~HBL(:,:,1);
    case {'feeders', 'feeder'}
        idx = ~HBL(:,:,2);
    case {'periphery'}
        idx = ~HBL(:,:,3);
    otherwise % Assume we use all the edges
        idx = tril(ones(size(SC)));
   end

    SC(logical(idx)) = 0;  %make lower tri = 0
    
    idx = SC(:) ~=0;       %index not zero values
    SCn = SC(idx);
    FCn = FC(idx);
    
    r = corr(SCn, FCn);

    
end % function calculate_corr_sc_fc()