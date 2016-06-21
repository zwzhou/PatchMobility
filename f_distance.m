function [ distance ] = f_distance( lx, ly, nx, ny )
% Calculate patch distance
%   size: nx^2 * ny^2

patchAmt = nx * ny ;
delx = lx / nx;  dely = ly / ny; % patch size

% initialization
patchnumber = 1:patchAmt;
patch = f_findcoor(patchnumber, nx, delx, dely);
xdis = bsxfun(@minus, patch.x2.', patch.x2);% patch x distance
ydis = bsxfun(@minus, patch.y2.', patch.y2);% patch y distance
distance = sqrt( xdis.^2 + ydis.^2 );% patch distance

end

