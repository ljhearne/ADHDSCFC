function setCmap(cols)
%sets colormap for matrices.

x = [0 255/2 255];
map = interp1(x/255,cols,linspace(0,1,255));
colormap(map)

end

