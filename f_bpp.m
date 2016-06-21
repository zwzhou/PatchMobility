function [ P ] = f_bpp( lx,ly,lz, Xs,Ys,Zs,S_o, Lx,Lz,nx,nz, Ax,Az, freq,c,rho_a,Nmodal,np,nq,nr )
% return Blocked Patch Pressure matrix, size: nPatch*freqNUM
%   nPatch = nx*nz;       freqNUM = length(freq);
% see test_bpp.m
%
%   lx,ly,lz: Cavity size
%   Xs,Ys,Zs: Source location
%   S_o: Source amplitude
%   Lx,Lz: panel size
%   nx,ny: panel meshgrid
%   Ax,Az: panel location. see f_findcoorL.m
%   freq: frequency vector.  size: 1*freqNUM
%   c: sound speed
%   Nmodal: total frequency order
%   np,nq,nr: frequency order in x,y and z axis.
%       !!  Nmodal must less than np*nq*nr  !!

fprintf('Calculating Blocked Patch Pressure\n')

% lx = 11.5;% Cavity size
% ly = 8.69;
% lz = 4.03; 

% Ax = 5.25; % panel location
% Az = 1.27; 
% Lx = 1.5;  % panel size
% Lz = 0.96;
% nx = 19;  % panel meshgrid
% nz = 13;
nPatch = nx*nz;
delx = Lx/nx; % patch size
delz = Lz/nz; 

% c = 343.6;  % m/s; sound speed

% Xs = 2;% Source location
% Ys = 4;
% Zs = 1; 
% S_o = 2; % Source amplitude

% freq = 50:2:1000 ; % frequency range
% freqNUM = length( freq );


%----------------------------------------------------------------------------------------
% fprintf('step %3d/4  cavity frequency\n',1)
% 
% % np = 350;
% % nq = 250;
% % nr = 150;
% % Nmodal = 100000 ;
% % [omegasqua,idx_c]   size: Nmodal*1
% [k_pqr_squa,idx_c] = f_cavityksqua(lx,ly,lz,Nmodal,np,nq,nr);
% 
% fprintf('step %3d/4  Apqr\n',2)
% 
% % source   size: Nmodal*1
% source = S_o* cos(idx_c.x*pi*Xs/lx).*cos(idx_c.y*pi*Ys/ly).*cos(idx_c.z*pi*Zs/lz);
% 
% lxcomp(1:Nmodal,1) = lx;
% lxcomp(idx_c.x~=0) = 0.5*lx;
% lycomp(1:Nmodal,1) = ly;
% lycomp(idx_c.y~=0) = 0.5*ly;
% lzcomp(1:Nmodal,1) = lz;
% lzcomp(idx_c.z~=0) = 0.5*lz;
% Npqr = lxcomp.*lycomp.*lzcomp;% Npqr   size: Nmodal*1
% 
% T_r = 10;
% eta_r = 2.2./(freq*T_r);
% cc = c* (1 + 1i*eta_r);
% kc_squa = (2*pi*freq./cc).^2;% kc_squa    size: 1*freqNUM
% 
% % A_pqr   size: Nmodal*freqNUM
% % A_pqr = bsxfun( @rdivide,intphis./Npqr,(bsxfun(@minus,kc_squa,k_pqr_squa)) );
% omega = 2*pi*freq;
% A_pqr = bsxfun(@times,1i*omega*rho_a,source./Npqr)./bsxfun(@minus,kc_squa,k_pqr_squa) ;
%--------------------------------------------------------------------------------------------

fprintf('    step %3d/3  Apqr\n',1)
[A_pqr,idx_c] = f_Apqr( lx,ly,lz, Xs,Ys,Zs,S_o, c,freq,rho_a,Nmodal,np,nq,nr );


fprintf('    step %3d/3  psi(x,y0,z)\n',2)

patchnumber = 1:nPatch;
[ coor ] = f_findcoorL(patchnumber,nx,delx,delz,Ax,Az);% size: 1*nPatch

psi1 = sin(idx_c.x*coor.x2*pi/lx) - sin(idx_c.x*coor.x1*pi/lx) ;% size: Nmodal*nPatch
psi1(idx_c.x==0,:) = 1;
psi1(idx_c.x~=0,:) = bsxfun(@rdivide,psi1(idx_c.x~=0,:),idx_c.x(idx_c.x~=0)) *lx/pi/delx;
psi2 = sin(idx_c.z*coor.y2*pi/lz) - sin(idx_c.z*coor.y1*pi/lz) ;% size: Nmodal*nPatch
psi2(idx_c.z==0,:) = 1;
psi2(idx_c.z~=0,:) = bsxfun(@rdivide,psi2(idx_c.z~=0,:),idx_c.z(idx_c.z~=0)) *lz/pi/delz;

fprintf('    step %3d/3  integral\n',3)

P = (psi1.*psi2).' * A_pqr;% size: nPatch*freqNUM

fprintf('complete\n')

end

