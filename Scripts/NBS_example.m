%NBS preperation script
clearvars

%Paths, functions and toolboxes
DataPath = '/Users/luke/Documents/Projects/ADHDStrucFunc/Data/';
DocsPath = '/Users/luke/Documents/Projects/ADHDStrucFunc/Docs/';

addpath(genpath('Functions'));
addpath(genpath('Toolbox'));

Atlas = 'Schaefer214';
%Atlas = 'Shen268';
%Atlas = 'Brainnetome_246';

%---------------------------------%
%---------------------------------%
% Load data
load([DataPath,Atlas,'/',Atlas,'_Info/',Atlas,'_coordinates.mat']);
load([DataPath,Atlas,'/',Atlas,'_Info/',Atlas,'parcellation_Yeo8Index.mat']);
load([DataPath,Atlas,'/','SC/CTRLSC118.mat']);
load([DataPath,Atlas,'/','SC/ADHDSC78.mat']);

%Matrices
MAT = cat(3,ADHDSC,CTRLSC);

%Design
load('/Users/luke/Documents/Projects/ADHDStrucFunc/Data/Schaefer214/NBS/Design_ADHD_CTRL_FD.txt');
DESIGN = Design_ADHD_CTRL_FD(:,1:2); % no need for FD

% Cordinates are COG

% save all
save NBS/NBS_MAT.mat MAT
save NBS/NBS_COG.mat COG
save NBS/NBS_Design.mat DESIGN