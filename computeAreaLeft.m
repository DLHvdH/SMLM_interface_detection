function [area] = computeAreaLeft(p1,p2,varargin)
% 
% computeAreaLeft
% 
% version:  1.1
% authors:  Dingeman van der Haven and Pim van der Hoorn
%
% Description:
% Computes the area corresponding to the part of a given rectangle to the
% left of a straight line through the points p1 and p2. If this line is
% horizontal then the area below it is computed. 
%
% REMARK: inputs of this function are not validated.
%
% Input:
% p1,p2     Two 1x2 arrays (row vectors) corresponding to the coordinates
%           of the points through which the line goes.
% extBound  (Optional)  The external bound of the region given as
%           [[xmin, ymin]; [xmax, ymax]]
%           Default     The unit square [[0,0];[1,1]]
%
% Output
% area      The area to the left of the line
%
%{
DEPENDENCIES:
 - lineIntersections
%}

%% Parse arguments

defaultSquare = [0 0; 1 1];

if isempty(varargin)
    extBounds = defaultSquare;
else
    extBounds = varargin{1};
end

% Compute the intersections of the line with the extended boundaries of the
% domain
[xt,xb,yl,yr] = lineIntersections(p1,p2,extBounds);

% Limit the intersection points to bounds of domain
xt = min(max(xt,extBounds(1,1)),extBounds(2,1));
xb = min(max(xb,extBounds(1,1)),extBounds(2,1));
yl = min(max(yl,extBounds(1,2)),extBounds(2,2));
yr = min(max(yr,extBounds(1,2)),extBounds(2,2));

% Domain bounds
yt = extBounds(2,2);
yb = extBounds(1,2);
xl = extBounds(1,1);
xr = extBounds(2,1);

% The line is given by the function f(z) = a z + b
% where a = (y2 - y1)/(x2 - x1) and b = y1 - x1(y2 - y1)/(x2 - x1).

% Note that if yr = yl, then we consider the area below the line to be
% 'left' of the line.

if yl >= yr
    area = 0.5*abs(yl-yr)*abs(xt-xb) + abs(xt-xl)*abs(yl-yb) + ...
        abs(xt - xr)*abs(yr - yb);
else
    area = 0.5*abs(yl-yr)*abs(xt-xb) + abs(xb-xl)*abs(yr-yb) + ...
        abs(xt - xl)*abs(yr - yt);
end

end

