function [points] = buildBandWidthSteps(extBounds,bandWidth,numSteps)
% 
% buildStepMatrices 
% 
% version:  1.1
% authors:  Dingeman van der Haven and Pim van der Hoorn
%
% Description:
% Support function for mleBoundaryEstimation.
%
% Creates a collection of points (as a M x 2 matrix) that are equally 
% spaced on the intersection of the boundaries of the provided region and
% bandwidth.
%
% REMARK: inputs of this function are not validated.
%
% Input:
% 
% extBounds     The external bound of the region given as
%               [[xmin, ymin]; [xmax, ymax]]
% bandWidth     The bandwidth for the points used to compute the 
%               [[minHorizontal,minVertical];[maxHorizontal,maxVertical]]
% numSteps      The number of steps between the points on the boundaries 
%               of the bandwidth
%
% Output:
% points        The coordinates of the numSteps + 1 boundaries equally 
%               spaced intervals.
%
%               bandwidth with two intersecting boundaries
%
% 					   p1
% 					+-|-|-|----------
% 					T	  |			|
% 				 p1 T	  |			|
% 					T_____|			|
% 					|	   			|
% 					|  	 _____		|
% 					|	|	  |		|
% 					|	|	  |		|
% 					----|-|-|-|------
% 						   p2
%
%                bandwidth with one intersecting boundary

% First, check with what boundaries of the region the bandWidth intersects 
% and compute the span or length of the intersecting region.
leftIntersect = false;
rightIntersect = false;
topIntersect = false;
bottomIntersect = false;

% Set the initinal span or length of the intersections to zero.
spanV = 0;
spanH = 0;

% Left boundary of the bandwidth intersects the the region?
if abs(bandWidth(1,1) - extBounds(1,1)) < (10*eps)
    leftIntersect = true;
    spanV = abs(bandWidth(2,2)-bandWidth(1,2));
end

% Right boundary of the bandwidth intersects the the region?
if abs(bandWidth(2,1) - extBounds(2,1))  < (10*eps)
    rightIntersect = true;
    spanV = abs(bandWidth(2,2)-bandWidth(1,2));
end
    
% Top boundary of the bandwidth intersects the the region?
if abs(bandWidth(2,2) - extBounds(2,2))  < (10*eps)
    topIntersect = true;
    spanH = abs(bandWidth(2,1)-bandWidth(1,1));
end

% Bottom boundary of the bandwidth intersects the the region?
if abs(bandWidth(1,2) - extBounds(1,2))  < (10*eps)
    bottomIntersect = true;
    spanH = abs(bandWidth(2,1)-bandWidth(1,1));
end

% Compute the step width. The exact location of the bandWidth does not 
% matter for this calculation. If the bandWidth is in one of the corners
% the result still holds because bandWidth is always rectangular.
stepWidth = (spanH+spanV)/numSteps;

% Check if the bandWidth intersects with multiple boundaries. If so, we
% need to compute the points differently.
if (leftIntersect + rightIntersect + topIntersect + bottomIntersect) == 1 
    % BandWidth region is NOT in one of the corners
    
    if leftIntersect
        % Iterate vertically on the left boundary
        points = [ bandWidth(1,1).*ones(numSteps+1,1), ...
            (bandWidth(1,2):stepWidth:bandWidth(2,2))' ];
    elseif rightIntersect
        % Iterate vertically on the right boundary
        points = [ bandWidth(2,1).*ones(numSteps+1,1), ...
            (bandWidth(1,2):stepWidth:bandWidth(2,2))' ];
    elseif topIntersect
        % Iterate hortizontally on the top boundary
        points = [ (bandWidth(1,1):stepWidth:bandWidth(2,1))' ,...
            bandWidth(2,2).*ones(numSteps+1,1) ];
    elseif bottomIntersect
        % Iterate hortizontally on the bottom boundary
        points = [ (bandWidth(1,1):stepWidth:bandWidth(2,1))' ,...
            bandWidth(1,2).*ones(numSteps+1,1) ];
    end
    
elseif (leftIntersect + rightIntersect + ...
        topIntersect + bottomIntersect) == 2 
    % BandWidth region is in one of the corners
    
    % Compute how many of the points fall on each intersection
    numH = floor(spanH/stepWidth-10*eps)+1;
    numV = floor(spanV/stepWidth-10*eps)+1;
 
    % Compute the increments in the x direction
    % Start increments at the point furthest away from the corner
    % Set the x increments for the points on the vertical intersction to
    % be zero.
    hSpaces = [ ( spanH:-stepWidth:0 )' ; zeros(numV,1) ];

    % Compute the increments in the y direction
    % Continue the increments at the corner
    % Take into account that the last x increment ended with an incomplete
    % stepWidth. xIncr(numH) gives the remainder.
    vSpaces = [ zeros(numH,1) ;...
        ( (stepWidth-hSpaces(numH)):stepWidth:spanV )' ];

    % Detect in which corner the bandWidth is. Then, 
    if leftIntersect
        % bandWidth is in left corner
        if topIntersect
        	% bandWidth is in top-left corner
            points = [bandWidth(1,1), bandWidth(2,2)] + ...
                [hSpaces, -vSpaces];
        else
            % bandWidth is in bottom-left corner
            points = bandWidth(1,:) + [hSpaces, vSpaces];
        end
    else
        % bandWidth is in right corner
        if topIntersect
            % bandWidth is in top-right corner
            points = bandWidth(2,:) + [-hSpaces, -vSpaces];
        else
            % bandWidth is in bottom-right corner
            points = [bandWidth(2,1), bandWidth(1,2)] + ...
                [-hSpaces, vSpaces];
        end
    end

elseif (leftIntersect + rightIntersect + ...
        topIntersect + bottomIntersect) >= 3
    error("Bandwidth region using step function intersects more than "+...
        "two boundaries. Please redefine the bandwidth region and try"+...
        "again.");
end


end