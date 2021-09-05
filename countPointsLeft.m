function [numPoints] = countPointsLeft(coords,p1,p2)
%
% countPointsLeft
%
% version:  1.1
% authors:  Dingeman van der Haven and Pim van der Hoorn
%
% Description:
% Counts the number of points in coords that are to the left of the straight line
% going through the points p1 and p2.
%
% REMARK: inputs of this function are not validated.
%
% Input:
% coords    The coordinates of the points as a Mx2 matrix
% p1,p2     Two points, represented as a 1x2 matrix (row vector) that 
%           determine the line for which the points to the left are 
%           counted.
%
% Output:
% numPoints The number of points in P that lie to the left of the line
%           through p1 and p2.


%% Setup the x and y coordinates of the two points
% We construct the straight line through the points (x1,y1) and (x2,y2).
% These are selected such that x1 <= x2. 

    if (p1(1)<=p2(1))
        x1 = p1(1);
        y1 = p1(2);

        x2 = p2(1);
        y2 = p2(2);
    else
        x1 = p2(1);
        y1 = p2(2);

        x2 = p1(1);
        y2 = p1(2);
    end
    
% The line is given by the function f(z) = a z + b
% where a = (y2 - y1)/(x2 - x1) and b = y1 - x1(y2 - y1)/(x2 - x1)

    
%% Count all points to the left of the line

% First we compute f(x) for all points in coords

    FP = (coords(:,1)-x1).*((y2-y1)/(x2-x1))+y1;

% Compute the logical vector that is 1 if a point is to the left of the
% line and zero otherwise.
%
% If y1 < y2 then a point p = (x,y) in coords is to the left of the line  
% when y >= f(x). When y1 => y2 the point is to the left when y <= f(x).

% Note that if y1 = y2, then we consider points below the line f(x) = b to 
% be to the 'left' of the line.

    if y1 < y2
        L1 = (coords(:,2)-FP)>=0;
    else
        L1 = (coords(:,2)-FP)<=0;
    end

% Then the number of points to the left is the total number of non-zero
% elements in L1.
    
    numPoints = sum(L1);

end

