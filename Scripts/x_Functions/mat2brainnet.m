function [NODE, EDGE] = mat2brainnet(matrix, MNI,network,size)
%writes a text file for the node and edge for the matrix
% you have to input the 
% (1) the matrix (weighted or otherwise depending on what you want)
% (2) a X by 3 MNI coordinate matrix (they have to be in order)
% (3) a X by 1 vector with network affilation, put 1's if none
% (4) a X by 1 vector with the size of each node

l=find(sum(matrix,2));
l=unique(l);

%% EDGE 
EDGE=matrix;
%EDGE=matrix(l,l);

fName = 'EDGE.edge';         %# A file name
fid = fopen(fName,'w');            %# Open the file

dlmwrite(fName,EDGE,'-append',...  %# Print the matrix
         'delimiter','\t',...
         'newline','pc');

%% NODE
MNI=round(MNI);
matrix(matrix~=0)=1; %binarizes matrix for NODE calculation


NODE=MNI;
NODE(:,4)=network;
NODE(:,5)=size;

%NODE=MNI(l,:);
%NODE(:,4)=network(l);
%NODE(:,5)=3;

str='-';

% Write the text file
f = fopen('NODE.node', 'w');
for n = 1:length(NODE)
    fprintf(f, '%d\t%d\t%d\t%d\t%d\t%s\n', NODE(n,1),NODE(n,2),NODE(n,3),NODE(n,4),NODE(n,5),str);
end
fclose(f);