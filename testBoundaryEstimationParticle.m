%
% testBoundaryEstimationParticle
%
% version:  1.1
% authors:  Dingeman van der Haven and Pim van der Hoorn
%
% Description:
% This script is a minimal example of how to use the 
% mleBoundaryEstimationParticle function to estimate the boundary between
% two homogeneous Poisson processes in the presence of one or more
% particles. 
%
% This script first generates a Poisson process on a unit square. The
% process consist of two homogeneous processes with different intensities
% mu1 and mu2, that are seperated by a line going through two points a and
% b.
% 							a
% 						-----------------
% 						|	\			|
% 						|	 \			|
% 						|	  \	   mu2	|
% 						|	   \		|
% 						|  mu1	\		|
% 						|		 \		|
% 						|		  \		|
% 						-----------------
% 									b
%
% Then we use the mleBoundary estimation to estimate the boundary (the line
% through a and b) from the generated data.
%
% The results a presented by plotting the points of the generated Poisson
% process, the true boundary in green and the estimated boundary in red.

%% Setup the parameters for the Poisson point process

domain = [[0,0];[1,1]]; % The region in which the Poisson Point Process 
                        % should be generated given as 
                        % [[xmin, ymin];[xmax, ymax]]

a = [0.4,1.0];            % These are the two points that determine the 
b = [0.6,0.0];            % boundary between the two processes.

M = 2000;               % Expected total number of generated points.

delta = 0.3;	% delta determines the fractional difference between the 
                % densities mu1 and mu2 in the left area A1 and the right 
                % area A2, respectively: mu2 = delta x mu1. 

centers = [     % Center coordinates of all the particles in 2xM format
    0.44 0.85
    0.45 0.5
    0.57 0.20
    ];


radii = [       % Radii of all the particles
    0.04
    0.10
    0.06
    ];



%% Simulate the process
P = generatePoisson2D(a,b,M,delta,domain);

%% Set bandwidths for estimating the boundary.

% The bandwidths are set to be a rectangle of width 2*LV and height 2*LH
% that have the point a or b at the center of one of their four boundaries

% 						
%
%                            LH
%    top bandwidth          |---|---|
%                               a
% 					   ___	-----------------
% 				   LV  _|_	|  	  	|		|
% 					   _|_	|_______|		|
% 							|       		|
% 							|	   			|
% 							|  	 _______	|
% 							|	|	    |	|
% 							|	|    	|	|   bottom bandwidth
% 							-----------------
%                                   b
% 							  

LV = 0.1;   %Half of the vertical bandwidth
LH = 0.1;   %Half of the horizontal bandwidth

topHmin = max(a(1)-LH,domain(1,1));
topHMax = min(a(1)+LH,domain(2,1));
topVmin = max(a(2)-LV,domain(1,2));
topVmax = min(a(2)+LV,domain(2,2));

bottomHmin = max(b(1)-LH,domain(1,1));
bottomHMax = min(b(1)+LH,domain(2,1));
bottomVmin = max(b(2)-LV,domain(1,2));
bottomVmax = min(b(2)+LV,domain(2,2));

topBand = [[topHmin,topVmin];[topHMax,topVmax]];
bottomBand = [[bottomHmin,bottomVmin];[bottomHMax,bottomVmax]];

%% Estimate the boundary using the non-Particle MLE

PnoPart = removeParticles(P,centers,radii);

[paramMax] = mleBoundaryEstimation(PnoPart,...
    topBand,bottomBand,domain,'IterationMethod','steps');

p1 = paramMax(1,:);
p2 = paramMax(2,:);

% Compute the x coordinates of the intersections of this line with the top
% and bottom of the region
[xt,xb,yl,yr] = lineIntersections(p1,p2,domain);

if (xt < domain(1,1))
    x1 = domain(1,1);
    y1 = yl;
elseif (xt > domain(2,1))
    x1 = domain(2,1);
    y1 = yr;
else
    x1 = xt;
    y1 = domain(2,2);
end

if (xb < domain(1,1))
    x2 = domain(1,1);
    y2 = yl;
elseif (xb > domain(2,1))
    x2 = domain(2,1);
    y2 = yr;
else
    x2 = xb;
    y2 = domain(1,2);
end


%% Estimate the boundary using the particle MLE

[paramMaxPart,Pnew] = mleBoundaryEstimationParticle(P,centers,radii,...
    topBand,bottomBand,domain,'IterationMethod','steps');

p1part = paramMaxPart(1,:);
p2part = paramMaxPart(2,:);

% Compute the x coordinates of the intersections of this line with the top
% and bottom of the region
[xt,xb,yl,yr] = lineIntersections(p1part,p2part,domain);

if (xt < domain(1,1))
    x1part = domain(1,1);
    y1part = yl;
elseif (xt > domain(2,1))
    x1part = domain(2,1);
    y1part = yr;
else
    x1part = xt;
    y1part = domain(2,2);
end

if (xb < domain(1,1))
    x2part = domain(1,1);
    y2part = yl;
elseif (xb > domain(2,1))
    x2part = domain(2,1);
    y2part = yr;
else
    x2part = xb;
    y2part = domain(1,2);
end


%% Initialize figure (without points within particles removed)
if false
    figure;

    % Populate figure with points in P the true line and the line corresponding 
    % to the maxmimum value of the likelihood function.
    hold on;

    h1 = scatter(P(:,1),P(:,2),'.');
    h2 = line([a(1) b(1)], [a(2) b(2)], 'Color','green','LineWidth',2);
    h3 = line([x1 x2], [y1 y2], 'Color','red','LineWidth',1,'LineStyle','--');
    h4 = line([x1part x2part], [y1part y2part], 'Color','red','LineWidth',1);
    h5 = viscircles(centers, radii,'Color','black','LineWidth',1);

    legend([h1 h2 h3 h4 h5],["Data","True boundary",...
        "MLE estimate","Particle MLE estimate",...
        "Particle boundary"]);
    box on; hold off;
end

%% Initialize figure (with points within particles removed)
figure;

% Populate figure with points in P the true line and the line corresponding 
% to the maxmimum value of the likelihood function.
hold on;

h1 = scatter(Pnew(:,1),Pnew(:,2),'.');
h2 = line([a(1) b(1)], [a(2) b(2)], 'Color','green','LineWidth',2);
h3 = line([x1 x2], [y1 y2], 'Color','red','LineWidth',1,'LineStyle','--');
h4 = line([x1part x2part], [y1part y2part], 'Color','red','LineWidth',1);
h5 = viscircles(centers, radii,'Color','black','LineWidth',1);

legend([h1 h2 h3 h4 h5],["Data","True boundary",...
    "MLE estimate","Particle MLE estimate",...
    "Particle boundary"]);
box on; hold off;
