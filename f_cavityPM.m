function [ Ysam,Ydif, Zsam,Zdif ] = f_cavityPM( lx,ly,lz, nx,ny, rho,cc, freq, Nmodal,np,nq,nr )
% return patch mobility of a small thickness air cavity: patchamt*patchamt*freqamt
% see: test_cavityPM.m
%  Ysam,Zsam: patch mobility and impedance in the same face
%  Ydif,Zdif: patch mobility and impedance in different faces
% Acoustic patch mobility of a cavity -> FIG.4 & FIG.5
disp('Calculating Cavity Acoustic Patch Mobility')

% cavity size
% lx = 1.5; 
% ly = 0.96; 
% lz = 0.01;

% rho = 1.293 ; % kg/m^3
% c=343.6; %m/s, speed of sound in air 20
% eta = 0.01;
% cc = c* sqrt(1+1i*eta);

% patch number of x and y
% nx = 19 ; 
% ny = 13 ;
patchamt = nx * ny ;%  amount of patch
delx = lx / nx ;  % patch dimension
dely = ly / ny ; 

% size: 1*patchamt
% patchnumber = 1:patchamt;
patch = f_ploc(patchamt,nx,delx,dely);


% % modal frequency order at x, y and z axis
% np = 30;
% nq = 30;
% nr = 5; % np nq & nr must >= 2 !
% % totle used order
% Nmodal = 1200 ;

% kpqr_squa and idx_c   size: Nmodal*1
[kpqr_squa,idx_c] = f_cavityksqua(lx,ly,lz,Nmodal,np,nq,nr);

% kc_squa
% freq = 100:1:600 ;
freqNUM = length( freq );
% kc_squa = ones(1,1,freqNUM);
omega = ones(1,1,freqNUM);
omega(1,1,:) = 2*pi*freq;
kc_squa = (omega./cc).^2;% size: 1*1*freqNUM

fprintf('    step %3d/4  prepare data martices \n',1)

% integral on patch i and j, size: Nmodal*patchamt 
intSij1 = sin(idx_c.x*patch.x2*pi/lx) - sin(idx_c.x*patch.x1*pi/lx);% Nmodal*patchamt
intSij1(idx_c.x==0,:) = delx;
intSij1(idx_c.x~=0,:) = bsxfun(@rdivide,lx*intSij1(idx_c.x~=0,:),pi*idx_c.x(idx_c.x~=0));
intSij2 = sin(idx_c.y*patch.y2*pi/ly) - sin(idx_c.y*patch.y1*pi/ly);% Nmodal*patchamt
intSij2(idx_c.y==0,:) = dely;
intSij2(idx_c.y~=0,:) = bsxfun(@rdivide,ly*intSij2(idx_c.y~=0,:),pi*idx_c.y(idx_c.y~=0));
intSij3 = ones(Nmodal,1);
intSij3(mod(idx_c.z,2)==1) = -1;
intSij_s = intSij1.*intSij2;
intSij_d = bsxfun(@times,intSij_s,intSij3);
clear intSij1 intSij2 intSij3


% Npqr    size: Nmodal*1
Npqr = lx*ly*lz*ones(Nmodal,1);
Npqr(idx_c.x~=0) = Npqr(idx_c.x~=0)/2;
Npqr(idx_c.y~=0) = Npqr(idx_c.y~=0)/2;
Npqr(idx_c.z~=0) = Npqr(idx_c.z~=0)/2;
% Np(idx_c.x==0) = lx;
% Np(idx_c.x~=0) = 0.5*lx;
% Nq(idx_c.y==0) = ly;
% Nq(idx_c.y~=0) = 0.5*ly;
% Nr(idx_c.z==0) = lz;
% Nr(idx_c.z~=0) = 0.5*lz;
% Npqr = (Np.*Nq.*Nr).';
% clear Np Nq Nr


% Assembly
% Z1 = bsxfun(@minus,kc_squa,kpqr_squa);% size: Nmodal*1*freqNUM
% Z2 = bsxfun(@rdivide,intSij_s,Npqr);% size: Nmodal*patchamt 
% Z3 = bsxfun(@rdivide,Z2,Z1);% Nmodal*patchamt*freqNUM
Z3 = bsxfun(@rdivide,bsxfun(@rdivide,intSij_s,Npqr),bsxfun(@minus,kc_squa,kpqr_squa));


fprintf('    step %3d/4  Zsam Assembling \n',2)

% Zsam and Zdif  size: patchamt*patchamt*freqNUM
Zsam = reshape(intSij_s.'*reshape(Z3,Nmodal,patchamt*freqNUM), patchamt,patchamt,freqNUM);
Zsam = bsxfun(@times,-1i*rho*omega,permute(Zsam,[2,1,3]));

fprintf('    step %3d/4  Zdif Assembling \n',3)

Zdif = reshape(intSij_d.'*reshape(Z3,Nmodal,patchamt*freqNUM), patchamt,patchamt,freqNUM);
Zdif = bsxfun(@times,-1i*rho*omega,permute(Zdif,[2,1,3]));
clear Z3

fprintf('    step %3d/4  Z Inversing \n',4)

% Ysam and Ydif  size: patchamt*patchamt*freqNUM
Ysam = ones(patchamt,patchamt,freqNUM);
Ydif = Ysam;
% tic
for ii = 1:freqNUM
%     disp(ii)
    Ysam(:,:,ii) = inv(squeeze(Zsam(:,:,ii)));
    Ydif(:,:,ii) = inv(squeeze(Zdif(:,:,ii)));
end
% toc

disp('Complete')

end

