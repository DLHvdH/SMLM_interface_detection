function [coords] = removeParticles(coords,centers,radii)
% Removes points based on a particle with a given position and radius. All
% points that are within the disc defined by the particle are removed.
%
% version:  1.1
% authors:  Dingeman van der Haven and Pim van der Hoorn
%
%   Input
%   coords      The coordinates of the points as a 2xM matrix
%   centers     The coordinates of the centers of the particles
%               as a 2xM matrix
%   radii       The radii of the discs corresponding to the particles
%               as a vector
%
%   Output
%   newCoords   The new set of coordinates with the points inside the
%               particle disc removed.
%

% Determine the total number of particles
nParticles = numel(radii);

% Initialise logical flags that indicates whether a point is outside all
% particles that have been given.
passTotal = true(size(coords,1),1);

% Loop over all particles
for i = 1:nParticles
    % Subtract particle center from all coordinates
    % Square thedifference
    % Sum along the second dimension (i.e. dx^2 + dy^2)
    % Take the square root to get the radial distances
    pass = sqrt(sum( (coords - centers(i,:)).^2 ,2)) > radii(i);
    % Update flag indicating if point is outside all particles
    passTotal = passTotal & pass;
end

% Apply the selection of the coordinates
coords = coords(passTotal,:);
        
end

