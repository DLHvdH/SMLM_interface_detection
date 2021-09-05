function [xt,xb,yl,yr] = lineIntersections(p1,p2,varargin)
%
% lineIntersections
%
% version:  1.1
% authors:  Dingeman van der Haven and Pim van der Hoorn
%
% Description: 
% Computes the x-coordinates where the line through the points p1 and p2
% intersects the top and bottom of the given region. In addition, it also
% computes the y-coordinates of the intersection of this line with the left 
% and right boundary of the region. 
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
% Output:
% xt  the x coordinate of the intersection of the line with the top
% xb  the x coordinate of the intersection of the line with the bottom
% yl  the y coordinate of the intersection of the line with the left
% yr  the y coordinate of the intersection of the line with the right
%
%
%         \ xt
%       _ _\_ __________
%           \|          |
%            \ yl       |
%            |\         |
%            | \        |
%            |  \       |
%            |___\______|
%               xb\           
%                  \    |
%                   \
%                    \  |
%                     \
%                      \| yr
%                       \

%% Parse arguments

defaultSquare = [0 0; 1 1];

if isempty(varargin)
    extBounds = defaultSquare;
else
    extBounds = varargin{1};
end

%% Compute the line intersections

% First we create the points (x1, y1) and (x2, y2), according to p1 and p2.
    x1 = p1(1);
    y1 = p1(2);

    x2 = p2(1);
    y2 = p2(2);
    
% The line is given by the function f(z) = a z + b
% where a = (y2 - y1)/(x2 - x1) and b = y1 - x1(y2 - y1)/(x2 - x1).
% So the inverse is given by: f^(-1)(y) = (y - y1)*(x2 - x1)/(y2 - y1) + x1
    
% x coordinate where f(x) crosses the top horizontal line 
% y = extBounds(2,2) is given by (f^(-1)(extBounds(2,2)))
    
xt = (extBounds(2,2)-y1)*((x2-x1)/(y2-y1))+x1;
    
% y coordinate of the intersection with the left boundary is
% f(extBounds(1,1))
    
yl =  (extBounds(1,1)-x1)*(y2 - y1)/(x2 - x1) + y1;

% x coordinate where f(x) crosses the bottom horizontal line 
% y = extBounds(1,2) is given by (f^(-1)(extBounds(1,2)))
    
xb = (extBounds(1,2)-y1)*((x2-x1)/(y2-y1))+x1;
    
% y coordinate of the intersection with the right boundary is
% f(extBounds(2,1))
    
yr = (extBounds(2,1)-x1)*(y2 - y1)/(x2 - x1) + y1;

end

