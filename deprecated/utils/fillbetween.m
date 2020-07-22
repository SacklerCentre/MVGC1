function h = fillbetween(x,yhi,ylo,color,edge,transparency)

% fill between two curves

if nargin<4; color        = 'b'; end % default color is blue
if nargin<5; edge         = 'k'; end % default edge color is black
if nargin<6; transparency = 0.5; end % default is to have a transparency of .5

assert(length(yhi)==length(ylo) && length(ylo)==length(x),'must have the same number of points in each vector');

filled = [yhi; flipud(ylo)];
x = [x; flipud(x)];
h = fill(x,filled,color);
set(h,'EdgeColor',edge,'FaceAlpha',transparency);
