% Double panel structure, Transmission loss, Patch-mobility method.
% Paper: "Prediction of transmission loss of double panels with a
% patch-moblility method, Jean-Louis Guyader"
clear

% cavity parameters
c = 343.6; % sound speed
rho_a = 1.204 ; % kg/m^3; density of air on 20 degrees Celsius
eta_a = 0.01;
cc = c* (1+1i*eta_a);

% big cavity. The paper Fig.8
Lx = 0.8; %11.5; % m
Ly = 1; %8.69;
Lz = 0.6; %4.03;

% Source in big cavity
S0 = 2; % Amplitude
Xs = 0.08; %2;% m; Location
Ys = 0.05; %4;
Zs = 0.1; %1;

% panel size and small cavity( Cavity C, The paper Fig.2 ) depth lz
lx = Lx; %1.5; % m
ly = Lz; %0.96;
h1 = 0.003; % m; thickness of the panel
h2 = 0.004; %0.0015;% m
lz = 0.1; %0.03;% m; depth of the small cavity
Ax = 0; %6-lx/2;% small cavity location
Ay = 0; %1.75-ly/2;

% panel parameters
rho_p = 2700 ; % kg/m^3
E = 7e10 ; % Pa 
eta_p = 0.01 ; % damping loss factor
mu = 0.3 ;  % Poisson ratio
Ec = E *( 1 + 1i*eta_p ); % complex Young's modulus
Dc1 = Ec*h1^3 /( 12*(1-mu^2) ); 
Dc2 = Ec*h2^3 /( 12*(1-mu^2) ); 

% panel patch
nXp = 20; %12;% meshgrid
nYp = 15; %20;
delXp = lx/nXp;% m
delYp = ly/nYp;
patchAmtp = nXp*nYp;
patchp = f_ploc(patchAmtp,nXp,delXp,delYp);

% small air cavity patch
nXa = nXp;%20; %12;
nYa = nYp;%15; %20;
delXa = lx/nXa;
delYa = ly/nYa;
patchAmta = nXa*nYa;
patcha = f_ploc(patchAmta,nXa,delXa,delYa);

% study frequency
freq = 1:1:500;
freqNum = length(freq);

disp('YP1')
YP1 = f_panelPM(lx,ly,h1,nXp,nYp,rho_p,Dc1,freq,800,60,60);
disp('YP2')
YP2 = f_panelPM(lx,ly,h2,nXp,nYp,rho_p,Dc2,freq,800,60,60);
[YAsam,YAdif,~,~] = f_cavityPM(lx,ly,lz,nXa,nYa,rho_a,cc,freq,2000,60,60,15);

bpp = f_bpp(Lx,Ly,Lz,Xs,Ys,Zs,S0,lx,ly,nXp,nYp,Ax,Ay,freq,c,rho_a,3000,300,280,230);
% V_B = f_m3vp2(YP,bpp);

% Zrad = f_radZ(lx,ly,nXp,nYp,rho_a,c,freq);


nxMax = max(nXp,nXa);
nyMax = max(nYp,nYa);
nAmt = nxMax*nyMax;


% if patchamtp~=namt
%     YP_e = zeros(namt,namt,freqNUM);
%     for ii = 1:freqNUM
%         YP_e(:,:,ii) = imresize(squeeze(YP(:,:,ii)),[namt namt],'lanczos3');
%     end
% else
%     YP_e = YP;
% end
% clear YP



A = zeros(2*nAmt);
V_B = zeros(2*nAmt,1);
mx1 = 1:nAmt;
mx2 = (nAmt+1):(2*nAmt);
F_BC = zeros(nAmt,freqNum);
F_DC = zeros(nAmt,freqNum);
V_D = zeros(nAmt,freqNum);
Irad = zeros(1,freqNum);
fprintf('Calculating coupling forces and Radiated power \nloop: %5d\n',1)
for ii = 1:freqNum
    if ~mod(ii,50)
        fprintf('loop: %5d\n',ii)
    end
    V_B(mx1) = squeeze(YP1(:,:,ii)) * squeeze(bpp(:,ii));% *delXp*delYp;
    A(mx1,mx1) = squeeze(YAsam(:,:,ii))+squeeze(YP1(:,:,ii));
    A(mx1,mx2) = squeeze(YAdif(:,:,ii));
    A(mx2,mx1) = squeeze(YAsam(:,:,ii));
    A(mx2,mx2) = squeeze(YAdif(:,:,ii))+squeeze(YP2(:,:,ii));
    F = A\V_B;
    F_BC(:,ii) = F(mx1);
    F_DC(:,ii) = F(mx2);
    
%     Radiated power
    V_D_t = squeeze(YP1(:,:,ii)) * F(mx2);
    V_D(:,ii) = V_D_t;
    Zrad = f_radZ(lx,ly,nXp,nYp,rho_a,c,freq(ii));
%     Irad(ii) = 0.5* ((V_D_t)'*squeeze(Zrad(:,:,ii))*V_D_t);
    Irad(ii) = 0.5* ((V_D_t)'*Zrad*V_D_t);
end
disp('Complete!')


quadPr = f_quadPr(Lz,Ly,Lx,Zs,Ys,Xs,S0,c,freq,rho_a,20000,350,250,150);


% transmission loss
% TL = Irad./(lx*ly*quadPr/4/rho_a/c);
TL = (lx*ly*quadPr/4/rho_a/c)./Irad;

figure(1)
fac = 1;%2e-5;
plot(freq,10*log10(abs(TL)/fac));
title('TL')
ylabel('db')
xlabel('Frequency (Hz)')

figure(2)
plot(freq,10*log10(Irad));
title('Radiated Power')
ylabel('db')
xlabel('Frequency (Hz)')

figure(3)
plot(freq,10*log10(quadPr));
title('quadPr')
ylabel('db')
xlabel('Frequency (Hz)')