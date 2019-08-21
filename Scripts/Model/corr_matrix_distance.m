function dcorr = corr_matrix_distance(A, B)
%% This function implements the correlation matrix distance
% as defined in Herdin et al 2005, IEEE
% ARGUMENTS:
%    A  -- an N x N correlation/covariance matrix
%    B  -- an N x N correlation/covariance matrix 
% OUTPUT:
%    dcorr -- a distance number between the two matrices
%
% AUTHOR:
% 
%     PSL 2019, QIMR Berghofer
%
% USAGE:
%{
    
A = rand(32);
B = rand(32);

dcorr = corr_matrix_distance(A, B);


%}

% NOTE: This function returns the correlation coefficient between the 
% elements of the lower triangular parts of the SC and FC matrices. 
% Only elements for which SC is nonzero are taken into account. 
% This distance is basically 1 - cosine similarity
% dcorr   --> 0 if A and B are different
% dcorr   --> 1 if A and B are equal


dcorr = 1 - trace(A*B) / (norm(A, 'fro')* norm(B, 'fro')); 

end % function corr_matrix_distance