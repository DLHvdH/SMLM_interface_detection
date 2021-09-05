function [P1,P2] = separatePointsByLine(coords,p1,p2,varargin)
%
% separatePointsByLine
%
% version:  1.1
% authors:  Dingeman van der Haven and Pim van der Hoorn
%
% Description:
% Divides the points in coords into two sets based on whether the point is
% to the left or right of the line going through p1 and p2.
%
% REMARK: inputs of this function are not validated.
%
% Input:
% coords    The coordinates of the points as a Mx2 matrix
% p1,p2     Two 1x2 arrays (row vectors) corresponding to the coordinates
%           of the points through which the line goes.
% extBound  (Optional)  The external bound of the region given as
%           [[xmin, ymin]; [xmax, ymax]]
%           Default     The unit square [[0,0];[1,1]]
%
% Output:
% P1        The point that are to the left of the line
% P2        The point that are to the right of the line
%
%{
DEPENDENCIES:
 - lineIntersections
%}

%% Parse arguments

defaultSquare = [0 1; 0 1];

if isempty(varargin)
    extBounds = defaultSquare;
else
    extBounds = varargin{1};
end


% Compute the x and y-coordinates at which the line through p1 and p2 
% intersects the boundaries of the region.

[xt,xb,yl,yr] = lineIntersections(p1,p2,extBounds);

% We now compute the logical vector that is 1 if a point is to the left 
% of the line and 0 otherwise.
%
% For this we have to make some case distinctions.

if (xt == xb)

    % In this case we have a straight line down. So we need to check
    % if the x coordinate is less or equal to the xt (= xb).

    L1 = (coords(:,1) <= xt);

elseif (yl == yr)

    % Here we have a horizontal line. In this case we define left of 
    % the line as meaning below the line.

    L1 = (coords(:,2) <= yl);

else

    % We are now in the situation where the line is given by a function
    % f(z) = az + b. To separate points to the left and right of this 
    % line we must first compute the coefficients a and b. For this we use 
    % that f(xt) = ymax and f(xb) = ymin.
    
    x1 = xt;
    y1 = extBounds(2,2);

    x2 = xb;
    y2 = extBounds(1,2);

    % The line that intersects these points is given by the function 
    % f(z) = a z + b where a = (y2 - y1)/(x2 - x1) and 
    % b = y1 - x1(y2 - y1)/(x2 - x1).

    a = (y2 - y1)/(x2-x1);
    b = y1 - x1*a;

    % We now compute the value f(x) for each x-coordinate in coords.
    FP = a.*coords(:,1)+b;

    % In order to separate the points we still need to make one more
    % case distinction.
    % 
    % If xt <= xb then a point p = (x,y) in coords is to the left of 
    % the line when y <= f(x). When xt > xb the point is to the left 
    % when y >= f(x). Points on the line are considered to be left of 
    % the line.

    if (xt < xb)
        L1 = (FP - coords(:,2))>=0;
    else
        L1 = (FP - coords(:,2))<=0;
    end  
end

% The logical vector L1 can now be used to select the points P1 left of the
% line and P2 right of the line.

P1 = coords(L1,:);
P2 = coords(~L1,:);

end

