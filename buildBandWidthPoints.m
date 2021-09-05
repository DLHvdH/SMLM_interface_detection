function [bottomPoints,topPoints] = buildBandWidthPoints(coords,...
    bottomBandWidth,topBandWidth)
% 
% buildBandWidthMatrices 
% 
% version:  1.1
% authors:  Dingeman van der Haven and Pim van der Hoorn
%
% Description:
% Support function for mleBoundaryEstimation.
%
% Selects from the points in coords those that fall into the provided 
% bottom and top bandwidth.
%
% REMARK: inputs of this function are not validated.
%
% Input:
% coords            The coordinates of the points as a Mx2 matrix
% topBandWidth      The bandwidth for the top as
%                   [[minHorizontal,minVertical];[maxHorizontal,maxVertical]]
% bottomBandWidth   The bandwidth for the bottom as 
%                   [[minHorizontal,minVertical];[maxHorizontal,maxVertical]]
%
% Output:
% bottomPoints      All point that fall into the bottom bandwidth
% topPoints         All point that fall into the top bandwidth
%
% 					top bandwidth		
% 						-----------------
% 						| x	  |			|
% 		    topPoints	|x  x |			|
% 						|____x|			|
% 						|	   			|
% 						|  	 _______	|
% 						|	| x	x   |	| bottomPoints
% 						|	|  x  x	|	|
% 						-----------------
% 						  bottom bandwidth


%% Select points for the top bandwidth

% Horizontal bounds
topHMin = topBandWidth(1,1);
topHMax = topBandWidth(2,1);

% Vertical bounds
topVMin = topBandWidth(1,2);
topVMax = topBandWidth(2,2);
    
% Find elements of P in the horizontal bandwidth. These are all points 
% p = (x,y) in P such that 
% topHMin <= x <= topHMax and topVMin <= y <= topVMax.
% First we construct the logical vector that corresponds to this criteria
TCond = (coords(:,1)>=topHMin & coords(:,1)<=topHMax &...
    coords(:,2)>=topVMin & coords(:,2)<=topVMax);

% Apply the condition TCond to the coordinates to obtain the top points.
topPoints = coords(TCond,:);

%% Select points in the vertical bandwidths

% Horizontal bounds
bottomHMin = bottomBandWidth(1,1);
bottomHMax = bottomBandWidth(2,1);

% Vertical bounds
bottomVMin = bottomBandWidth(1,2);
bottomVMax = bottomBandWidth(2,2);
    
% Find elements of P in the bottom bandwidth. These are all points 
% p = (x,y) in P such that 
% bottomHMin <= x <= bottomHMax and bottomVMin <= y <= bottomVMax.
% First we construct the logical vector that corresponds to this criteria
BCond = (coords(:,1)>=bottomHMin & coords(:,1)<=bottomHMax &...
    coords(:,2)>=bottomVMin & coords(:,2)<=bottomVMax);

% Apply the condition BCond to the coordinates to obtain the bottom points.
bottomPoints = coords(BCond,:);

    
end

