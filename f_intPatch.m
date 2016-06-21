function [ int_patch ] = f_intPatch( lx, ly, nx, ny, idx_p )
% integral over patch 
%    idx_p: panel index. idx_p.x and idx_p.y size: Nmodal*1
%   size: Nmodal*patchAmt

patchAmt = nx * ny;
delx = lx/nx;
dely = ly/ny;
patchnumber = 1:patchAmt;
% patch location
ploc = f_findcoor(patchnumber,nx,delx,dely);% size: 1*patchAmt

% tic
ppilx2 = idx_p.x*pi/lx *ploc.x2;% size: Nmodal*patchAmt
ppilx1 = idx_p.x*pi/lx *ploc.x1;

ppily2 = idx_p.y*pi/ly *ploc.y2;
ppily1 = idx_p.y*pi/ly *ploc.y1;

int_p1 = cos(ppilx2) - cos(ppilx1);% Nmodal*patchAmt
int_p2 = cos(ppily2) - cos(ppily1);

int_times = int_p1 .* int_p2;
idx_times = idx_p.x .* idx_p.y;
tmp = lx*ly/(pi^2);

int_patch = tmp * bsxfun(@rdivide, int_times, idx_times);% Nmodal*patchamt

% toc

end

