function [paramMax] = mleBoundaryEstimation(coords,varargin)
%
% mleBoundaryEstimation
%
% version:  1.1
% authors:  Dingeman van der Haven and Pim van der Hoorn
%
% Description:
% Estimates the boundary between two regions of a Poisson Point Process.
% Please see corresponding README for explanation of usage.
%
% Input
% coords            The coordinates of the points as a Mx2 matrix
% TopBandWidth      (Optional)  The bandwidth for the top points used to 
%                               compute the ML as 
%                               [[minHorizontal,minVertical]; [maxHorizontal,maxVertical]]
%                   Default     The top half of the region
% BottomBandWidth   (Optional)  The bandwidth for the bottom points used to
%                               compute the ML as 
%                               [[minHorizontal,minVertical];[maxHorizontal,maxVertical]]
%                   Default     The bottom half of the region
% ExtBound          (Optional)  The external bound of the region given as
%                               [[xmin, ymin]; [xmax, ymax]]
%                   Default     The unit square [[0,0];[1,1]]
% IterationMethod   (Parameter) Determines the iteration method for finding
%                               the maximum likelihood. Two version are
%                               available: steps or points. Please see
%                               corresponding README for information on the
%                               two methods.
%                   Default     steps
%
% Output:
% paramMax  An 2x2 matrix [p1; p2] corresponding to the two points
%           that maximized the likelihood function.
%
%{
DEPENDENCIES:
 - buildBandWidthMatrices
 - buildStepMatrices
 - computeAreaLeft
 - countPointsLeft
 - lineIntersections
 - mleBoundaryEstimation
%}

%% Parse arguments

checkPoints = @(x) isnumeric(x) && size(x,2)==2;
checkNonNegative = @(x) isnumeric(x) && x >= 0;
check2DMatrix = @(x) isnumeric(x) && (size(x,1) == 2 && size(x,2) == 2);

defaultIteration = 'steps';
expectedIterations = {'steps','points'};
checkMethod = @(x) mustBeMember(x,expectedIterations);

defaultIterationSteps = 50;
defaultTopBand = [[0,0.5];[1,1]];
defaultBottomBand = [[0,0];[1,0.5]];
defaultSquare = [0 0; 0 0];

p = inputParser;

addRequired(p,'coords',checkPoints);

addOptional(p,'TopBandwidth',defaultTopBand,check2DMatrix);
addOptional(p,'BottomBandwidth',defaultBottomBand,check2DMatrix);
addOptional(p,'ExtBounds',defaultSquare,check2DMatrix);

addParameter(p,'IterationMethod',defaultIteration,checkMethod);
addParameter(p,'IterationSteps',defaultIterationSteps,checkNonNegative)

parse(p,coords,varargin{:});

%% Set boundaries
% If both the xmin and xmax are zero we assume that no boundaries are
% provided. In this case we take the boundaries to be the smallest
% rectangle containing the points, extended by an offset of 0.05.

extBounds = p.Results.ExtBounds;
if (extBounds(1,1) == 0 && extBounds(2,1) == 0)
    offset = 0.05;
    extBounds = [[min(coords(:,1))-offset,min(coords(:,2))-offset];...
        [max(coords(:,1))+offset,max(coords(:,2))+offset]];
end

%% Determine bandwidths

bottomBandWidth = p.Results.BottomBandwidth;
topBandWidth = p.Results.TopBandwidth;

% Make sure the bandwidths are contained within the specified region

if (max(bottomBandWidth(1,1),bottomBandWidth(2,1)) > extBounds(2,1)...
    || min(bottomBandWidth(1,1),bottomBandWidth(2,1)) < extBounds(1,1)...
    || max(bottomBandWidth(1,2),bottomBandWidth(2,2)) > extBounds(2,2)...
    || min(bottomBandWidth(1,2),bottomBandWidth(2,2)) < extBounds(1,2))
    
    warningStr = ['The provided bottom bandwidth is not contained'...
        'within the specified region. The bottom bandwidths will be'...
        'recomputed to fit inside the region'];
    warning(warningStr);
    
    % Correct the bottom bandwidth
    bottomBandWidth(1,1) = min(max(bottomBandWidth(1,1),extBounds(1,1)),...
        extBounds(2,1));
    bottomBandWidth(2,1) = max(min(bottomBandWidth(1,2),extBounds(2,1)),...
        extBounds(1,1));
    bottomBandWidth(1,2) = min(max(bottomBandWidth(2,1),extBounds(1,2)),...
        extBounds(2,2));
    bottomBandWidth(2,2) = max(min(bottomBandWidth(2,2),extBounds(2,2)),...
        extBounds(1,2));
end

if (max(topBandWidth(1,1),topBandWidth(2,1)) > extBounds(2,1)...
    || min(topBandWidth(1,1),topBandWidth(2,1)) < extBounds(1,1)...
    || max(topBandWidth(1,2),topBandWidth(2,2)) > extBounds(2,2)...
    || min(topBandWidth(1,2),topBandWidth(2,2)) < extBounds(1,2))

    warningStr = ['The provided top bandwidth is not contained'...
        'within the specified region. The top bandwidths will be'...
        'recomputed to fit inside the region'];
    warning(warningStr);

    % Correct the top bandwidth
    topBandWidth(1,1) = min(max(topBandWidth(1,1),extBounds(1,1)),...
        extBounds(2,1));
    topBandWidth(2,1) = max(min(topBandWidth(1,2),extBounds(2,1)),...
        extBounds(1,1));
    topBandWidth(1,2) = min(max(topBandWidth(2,1),extBounds(1,2)),...
        extBounds(2,2));
    topBandWidth(2,2) = max(min(topBandWidth(2,2),extBounds(2,2)),...
        extBounds(1,2));
end


%% Build the matrix of points in the top and bottom through which the 
%% procedure has to iterate

if(strcmp(p.Results.IterationMethod,'points')==1)
    
    [bottomPoints,topPoints] = buildBandWidthPoints(coords,...
        bottomBandWidth,topBandWidth);
    
else
    
    numSteps = p.Results.IterationSteps;
    
    topPoints = buildBandWidthSteps(extBounds,topBandWidth,numSteps);
    bottomPoints = buildBandWidthSteps(extBounds,bottomBandWidth,numSteps);
    
end

sizeBottom = size(bottomPoints,1);
sizeTop = size(topPoints,1);
    
%% Compute area of the region and total number of points
area = abs(extBounds(1,1)-extBounds(2,1))*...
    abs(extBounds(1,2)-extBounds(2,2));
nPtotal = size(coords,1);

%% Initialize variables

maxValue = -Inf;
paramMax = zeros(2,2);

%% Run over all pairs in the bandwidths and compute the loglikelihood

% The procedure considers lines that go through point pairs (p1,p2) where 
% p1 is in the top bandwidth and p2 is in the bottom bandwidth.

for i=1:sizeTop
    for j=1:sizeBottom     
        p1 = topPoints(i,:);
        p2 = bottomPoints(j,:);
        
        % Compute the loglikelihood function based on the line through the
        % points p1 and p2
        
        nP1 = countPointsLeft(coords,p1,p2);        % #points to the left
        nP2 = nPtotal-nP1;                          % #points to the right
        areaLeft = computeAreaLeft(p1,p2,extBounds);    % area left
        areaRight = area-areaLeft;                      % area right

        % The loglikelihood
        value = nP1*log(nP1/areaLeft) + nP2*log(nP2/areaRight); 

        % Update the current maximum value and the corresponding points
        % used for the line.

        if value >= maxValue
            maxValue = value;
            paramMax = [p1; p2];
        end

    end
end

end
