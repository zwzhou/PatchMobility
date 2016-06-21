% clc;
clear;
% panel patch mobility. FIG.3
% disp('Calculating Panel Patch Mobility')
% panel size
lx = 1.5;
ly = 0.96;
p_area = lx*ly;
h = 0.002;

% panel property
rho_p = 2700 ; % kg/m^3; density of the panel
E = 70e9 ; % Pa;  Young's modulus
eta = 0.01 ; % damping loss factor
mu = 0.3 ;  % Poisson ratio
Ec = E*(1 + 1i*eta); % complex Young's modulus
Dc = Ec*h^3 /( 12*(1-mu^2) ); 

M_pq = rho_p*h*lx*ly/4 ;

% panel meshgrid
nx = 19;
ny = 13;
patchamt = nx*ny;
np = 100;
nq = 100;
Nmodal = 200;

freq = 100:1:500;

yp = f_panelPM( lx,ly,h, nx,ny, rho_p,Dc, freq, Nmodal,np,nq );

%---------------------------------------------------------------------------------------------
% delx = lx/nx;
% dely = ly/ny;
% delarea = delx*dely;
% 
% fprintf('step %3d/2  integral patch i and j \n',1)
% 
% % patch location
% patchnumber = 1:patchamt;
% ploc = f_findcoor(patchnumber,nx,delx,dely);% size: 1*patchamt
% 
% % modal frequency
% % [omega_c_squa,idx_p]  size: Nmodal*1
% [omegac_pq_squa,idx_p] = f_plateOmegaSqua_simspt(lx,ly,Dc,rho_p*h,Nmodal,np,nq);
% 
% % caculated frequency
% 
% freqNUM = length(freq);
% omega = ones(1,1,freqNUM);
% omega(1,1,:) = 2*pi*freq;
% % omega_squa = ones(1,1,freqNUM);
% omega_squa = omega.^2;% size: 1*1*freqNUM
% 
% % fprintf('step %3d/3  Apqr\n',1)
% 
% % integral patch i and j
% intij1 = cos(idx_p.x*ploc.x2*pi/lx) - cos(idx_p.x*ploc.x1*pi/lx);% Nmodal*patchamt
% intij2 = cos(idx_p.y*ploc.y2*pi/ly) - cos(idx_p.y*ploc.y1*pi/ly);
% intij = -p_area/pi^2*bsxfun(@rdivide,(intij1.*intij2),idx_p.x.*idx_p.y);% Nmodal*patchamt
% clear intij1 intij2
% 
% fprintf('step %3d/2  Assembling \n',2)
% 
% % yp1 = bsxfun(@minus,omegac_pq_squa,omega_squa);% size: Nmodal*1*freqNUM
% % yp2 = bsxfun(@rdivide,intij,yp1);% size: Nmodal*patchamt*freqNUM
% % % equal to the two lines above.  size: Nmodal*patchamt*freqNUM
% % % yp2 = bsxfun(@rdivide,intij,bsxfun(@minus,omegac_pq_squa,omega_squa));
% % yp3 = reshape(yp2,Nmodal,patchamt*freqNUM);% size: Nmodal*(patchamt*freqNUM); two dimensional maxtrix
% 
% % equal to the three lines above.
% % size: Nmodal*(patchamt*freqNUM); two dimensional maxtrix
% yp3 = reshape(bsxfun(@rdivide,intij,bsxfun(@minus,omegac_pq_squa,omega_squa)),Nmodal,patchamt*freqNUM);
% yp = bsxfun(@times,1i/delarea^2/M_pq*omega,reshape(intij.'*yp3,patchamt,patchamt,freqNUM));
% clear yp3
% 
% disp('complete')
%-----------------------------------------------------------------------------------------------------------


figure(1)
% tr_patchmobility = zeros(1,freqamt);
tr_patchmobility = squeeze(yp(42,72,:));
plot(freq,20*log10(abs(tr_patchmobility)),'b')
hold on
% ip_patchmobility = zeros(1,freqamt);
ip_patchmobility = squeeze(yp(42,42,:));
plot(freq,20*log10(abs(ip_patchmobility)),'m')
legend('transfer patch mobility','input patch mobility')
xlabel('Frequency (HZ)');
ylabel('Magnitude (dB)');





