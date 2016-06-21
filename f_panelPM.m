function [ yp ] = f_panelPM( lx,ly,h, nx,ny, rho_p,Dc, freq, Nmodal,np,nq )
% return patch mobility martix of a panel: patchamt*patchamt*freqNUM
%   lx, ly & h: length, width and thickness of the plate
%   nx & ny: patch number along length and width
%   rho_p: density of the plate
%   Dc: complex bending stiffness 
%   freq: frequency vector. Example: freq=1:0.1:100
%      see: test_panelPM.m


% disp('Calculating Panel Patch Mobility')
% panel size
% lx = 1.5;
% ly = 0.96;
p_area = lx*ly;
% h = 0.002;

% panel property
% rho = 2700 ; % kg/m^3; density of the panel
% E = 6.9e10 ; % Pa;  Young's modulus
% eta = 0.01 ; % damping loss factor
% mu = 0.334 ;  % Poisson ratio
% Ec = E *( 1 + 1i*eta ); % complex Young's modulus
% Dc = Ec *h^3 /( 12*(1-mu^2) ); 

M_pq = rho_p*h*lx*ly/4 ;

% panel meshgrid
% nx = 19;
% ny = 13;
patchAmt = nx*ny;
delx = lx/nx;
dely = ly/ny;
delarea = lx*ly/patchAmt;
% np = 50;
% nq = 50;
% Nmodal = 500;

% fprintf('step %3d/2  integral patch i and j \n',1)

% patch location
patchnumber = 1:patchAmt;
ploc = f_findcoor(patchnumber,nx,delx,dely);% size: 1*patchAmt

% modal frequency
% [omega_c_squa,idx_p]  size: Nmodal*1
[omegac_pq_squa,idx_p] = f_plateOmegaSqua_simspt(lx,ly,Dc,rho_p*h,Nmodal,np,nq);

% caculated frequency
% freq = 1:1:600;

% tic
pi_lx = pi/lx;
pi_ly = pi/ly;
intij1 = cos(idx_p.x*ploc.x2*pi_lx) - cos(idx_p.x*ploc.x1*pi_lx);% Nmodal*patchAmt
intij2 = cos(idx_p.y*ploc.y2*pi_ly) - cos(idx_p.y*ploc.y1*pi_ly);
int_patch = -p_area/(pi^2)*bsxfun(@rdivide,(intij1.*intij2),idx_p.x.*idx_p.y);% Nmodal*patchamt
clear intij1 intij2
% toc

freqNum = length(freq);
% if freqNum>1
    disp('Calculating Panel Patch Mobility')
	omega = ones(1,1,freqNum);
	omega(1,1,:) = 2*pi*freq;
	omega_squa = omega.^2;% size: 1*1*freqNUM
	% fprintf('step %3d/3  Apqr\n',1)
	% paper function (31)
	% integral patch i and j

	
% 	fprintf('step %3d/2  Assembling \n',2)
	
	% yp1 = bsxfun(@minus,omegac_pq_squa,omega_squa);% size: Nmodal*1*freqNUM
	% yp2 = bsxfun(@rdivide,intij,yp1);% size: Nmodal*patchamt*freqNUM
	% % equal to the two lines above.  size: Nmodal*patchamt*freqNUM
	% % yp2 = bsxfun(@rdivide,intij,bsxfun(@minus,omegac_pq_squa,omega_squa));
	% yp3 = reshape(yp2,Nmodal,patchamt*freqNUM);% size: Nmodal*(patchamt*freqNUM); two dimensional maxtrix
	
	% equal to the three lines above.
	% size: Nmodal*(patchamt*freqNUM); two dimensional maxtrix
	yp3 = reshape(...
        bsxfun(@rdivide,int_patch, bsxfun(@minus,omegac_pq_squa,omega_squa) ),...
        Nmodal,patchAmt*freqNum);
	yp = bsxfun( @times, 1i/delarea^2/M_pq*omega ,...
        reshape(int_patch.'*yp3,patchAmt,patchAmt,freqNum) );
	% clear yp3
	
	fprintf('complete\n')
% else
%     omega = 2*pi*freq;
% end

end

