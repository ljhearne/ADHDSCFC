function surficeEdgeVis(MAT,COG,net,Nsize,thresh,out)
%function surficeVisbundleVis(MAT,COG,net,Nsize,thresh,out)
%surficeVis - generates brain vis for StrokeNet project using surfice and
%R. Mofidy this script to change all renderings in the project.
%MAT =node by node matrix
%COG =coordinates (node by 3)
%net =network affilations (or ones)
%size=size of nodes (or ones)
%thresh=lower and upper threshold for edges.
%out = output path and prefix

% generate NODE and EDGE files
mat2brainnet(MAT,COG,net,Nsize,[out,'.EDGE'],[out,'.NODE']);

edgefile = [out,'.EDGE'];
nodefile = [out,'.NODE'];
mesh = '/Applications/Surfice/sample/mni152_2009.mz3';
exe = '/Applications/Surfice/surfice.app/Contents/MacOS/surfice';

%list of commands
BACKCOLOR = 'BACKCOLOR(255, 255, 255);';
MESHLOAD = ['meshload(''',mesh,''');'];
NODELOAD = ['nodeload(''',nodefile,''');'];
NODESIZE = 'NODESIZE(1.0, true);';
EDGESIZE = 'EDGESIZE(0, true);';
EDGETHRESH = ['EDGETHRESH(',num2str(thresh(1)),',',num2str(thresh(2)),');'];
SHADERXRAY = 'SHADERXRAY(0.50, 0);';
NODETHRESH = 'NODETHRESH(1, 5);';
ORIENT = 'ORIENTCUBEVISIBLE(false);';
%options
CB = 'colorbarvisible(true);';

%plot colorbar for post processing

% plot three angles
VIEW{1} = 'VIEWAXIAL(true);';
VIEW{2} = 'VIEWSAGITTAL(true);';
VIEW{3} = 'VIEWSAGITTAL(true);AZIMUTH(180);';

clear outfile
outfile{1} = [out,'_surficeEdgeAxial.bmp'];
outfile{2} = [out,'_surficeEdgeRight.bmp'];
outfile{3} = [out,'_surficeEdgeLeft.bmp'];
CB = 'colorbarvisible(false);';

for i = 1:3
    cmd = [exe,' -S "begin RESETDEFAULTS;',...
        BACKCOLOR,MESHLOAD,NODELOAD,NODESIZE,EDGESIZE,SHADERXRAY,NODETHRESH,...
        ORIENT,CB,VIEW{i},...
        'SAVEBMP(''',outfile{i},''');',...
        'quit;',...
        'end."'];
    system(cmd);
end
