%   COMPUTE ANALYTICAL FC (PCC) FROM STRUCTURAL DATA for a given value of
%   the global coupling parameter c.
%
%   This function computes the analytical functional connectivity, measured 
%   as Pearson Correlation Coefficient, starting from structural data. 
%   The analytical operation is obtained under the assumption of an 
%   underlying stochastic linear model.
%
%   Input:
%   empSC --> is the Nodes x Nodes structural connectivity matrix, where Nodes is the 
%             number of regions in the parcellation
%   global_coupling --> is value of the global coupling parameter
%   sigma_noise --> standard deviation of white noise which runs the system
%               --> can be a scalar or a vector of length n
%   
%   Output:
%   aPCC --> is a matrix Nodes x Nodes whose entries represent an estimate 
%            of the functional connectivity for the given value of the 
%            paramter global coupling

% Original from Saggio et al. 2016
% Modified by PSL 2019 QIMR to handle nonuniform sigma values

function  aPCC = calculate_analyticalFCPCC(empSC, global_coupling, sigma_noise)

% number of regions
number_of_nodes = size(empSC,1);                     

% coupling matrix for the linear system
A = - eye(number_of_nodes) + global_coupling * empSC; 

% check sigma is the same or different for every node
if ismatrix(sigma_noise)
    sigma_mat = (sigma_noise*sigma_noise')/2;
else
    disp('Sigma takes a uniform value')
    sigma_mat = (sigma_noise^2)/2;
end

aCOV= -sigma_mat * inv(A);      % analytic covariance matrix

aPCC = zeros(number_of_nodes);  % analytic PCC matrix

for i=1:number_of_nodes
    for j=1:number_of_nodes
        aPCC(i,j) = aCOV(i,j)/sqrt(aCOV(i,i)*aCOV(j,j));
    end
end

end % function calculate_analyticalFCPCC()

