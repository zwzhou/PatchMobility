function [ bpp ] = f_bpp_alpha( lx,ly,lz, Xs,Ys,Zs,S_0, Lx,Lz,nx,nz, Ax,Az, freq,c_air )
% retune bpp: patchamt*freqamt matrix, blocked patch pressure of the patch number versus frequency
%   Detailed explanation goes here
%clear;clc;
% VII. SOURCE ROOM MODELING --> B.Blocked patch pressure --> FIG.10
% c_air = 343.6; % sound speed
%S_0 = 2;
%lx=11.5; ly=8.69; lz=4.03;% cavity dimension
%Xs = 2; Ys = 4; Zs = 1;% source location
%Lx = 1.5;   Lz = 0.96;
%nx = 19;    nz = 13;
%Ax = 5.245; Az = 1.27;
delx = Lx/nx; delz = Lz/nz; % patch size
patchamt = nx*nz;
% patch(patchamt) = struct('x1',0,'x2',0,'z1',0,'z2',0,'xc',0,'yc',0,'area',0);
% for loop= 1:patchamt
%     [x1,x2,z1,z2] = findcoor(loop,nx,delx,delz);
%     xc = (x1+x2) / 2;
%     zc = (z1+z2) / 2;
%     p_area = abs(x2-x1)*abs(y2-y1);
%     patch(loop) = struct('x1',x1,'x2',x2,'z1',z1,'z2',z2,'xc',xc,'zc',zc,'area',p_area);
% end
pi_lx = pi/lx;  pi_ly = pi/ly;  pi_lz = pi/lz;

np = 70;   nq = 55;   nr = 45;

%freq = 10:100 ;
freqamt = length( freq );
T_r = 10;

% (1)  p,r > 0     q>0
[p1,q1,r1] = meshgrid(1:np,1:nq,1:nr);
p1pi_lx = p1*pi_lx;
q1pi_ly = q1*pi_ly;
r1pi_lz = r1*pi_lz;
N_pqr_n0 = lx*ly*lz /8;
k_pqr_sq_1 = (p1pi_lx).^2 + (q1pi_ly).^2 +(r1pi_lz).^2;
A_pqr_comp_1 = S_0.*cos(Xs*p1pi_lx).*cos(Ys*q1pi_ly).*cos(Zs*r1pi_lz)./ N_pqr_n0./ (p1pi_lx.*r1pi_lz) ;
% ipo_1 = 1./ (p1pi_lx.*r1pi_lz);

% (2)  p=0, r>0     q>0
[q2,r2] = meshgrid(1:nq,1:nr);
q2pi_ly = q2*pi_ly;
r2pi_lz = r2*pi_lz;
N_pqr_10 = lx*ly*lz /4;
k_pqr_sq_2 =  (q2pi_ly).^2 +(r2pi_lz).^2;
A_pqr_comp_2 = S_0.*cos(Ys*q2pi_ly).*cos(Zs*r2pi_lz)./ N_pqr_10./ r2pi_lz ;
Pxcomp_2 = delx ;

% (3)  p>0, r=0     q>0
[p3,q3] = meshgrid(1:np,1:nq);
p3pi_lx = p3*pi_lx;
q3pi_ly = q3*pi_ly;
k_pqr_sq_3 = (p3pi_lx).^2 + (q3pi_ly).^2 ;
A_pqr_comp_3 = S_0.*cos(Xs*p3pi_lx).*cos(Ys*q3pi_ly)./ N_pqr_10./ p3pi_lx ;
Pzcomp_3 = delz;

% (4)  p,r = 0     q>0
q4 = 1:nq ;
q4pi_ly = q4*pi_ly;
k_pqr_sq_4 = q4pi_ly.^2 ;
N_pqr_20 = lx*ly*lz /2;
A_pqr_comp_4 = S_0.*cos(Ys*q4pi_ly) ./N_pqr_20 ;
Pxcomp_4 = delx;
Pzcomp_4 = delz;

% (5)  p,r > 0     q=0
[p5,r5] = meshgrid(1:np,1:nr);
p5pi_lx = p5*pi_lx;
r5pi_lz = r5*pi_lz;
k_pqr_sq_5 = (p5pi_lx).^2 +(r5pi_lz).^2;
A_pqr_comp_5 = S_0.*cos(Xs*p5pi_lx).*cos(Zs*r5pi_lz)./ N_pqr_10./ (p5pi_lx.*r5pi_lz) ;

% (6)  p=0, r>0     q=0
r6 = 1:nr ;
r6pi_lz = r6*pi_lz;
k_pqr_sq_6 = r6pi_lz.^2 ;
A_pqr_comp_6 = S_0.*cos(Zs*r6pi_lz)./ N_pqr_20./ r6pi_lz ;
Pxcomp_6 = delx;

% (7)  p>0, r=0     q=0
p7 = 1:np ;
p7pi_lx = p7*pi_lx;
k_pqr_sq_7 = p7pi_lx.^2;
A_pqr_comp_7 = S_0.*cos(Xs*p7pi_lx)./ N_pqr_20./ p7pi_lx ;
Pzcomp_7 = delz;

% (8)  p,q,r = 0
N_pqr_a0 = lx*ly*lz;
A_pqr_comp_8 = S_0/N_pqr_a0;
Pxcomp_8 = delx;
Pzcomp_8 = delz;

P = ones(patchamt,freqamt);
for loop = 1:freqamt
    fprintf('%i\n',freq(loop))
    eta_r = 2.2/(freq(loop)*T_r);
    cc = c_air* sqrt(1 + 1i*eta_r);
    omega = 2* pi* freq( loop );
    kc_sq = ( omega / cc )^2;
    % P_sum = 0;
    
    % (1)  p,r > 0     q>0
    A_pqr_1 = A_pqr_comp_1 ./ (kc_sq - k_pqr_sq_1) ;
    % (2)  p=0, r>0     q>0
    A_pqr_2 = A_pqr_comp_2 ./ (kc_sq - k_pqr_sq_2) ;
    % (3)  p>0, r=0     q>0
    A_pqr_3 = A_pqr_comp_3 ./ (kc_sq - k_pqr_sq_3) ;
    % (4)  p,r = 0     q>0
    A_pqr_4 = A_pqr_comp_4 ./ (kc_sq - k_pqr_sq_4) ;
    % (5)  p,r > 0     q=0
    A_pqr_5 = A_pqr_comp_5 ./ (kc_sq - k_pqr_sq_5) ;
    % (6)  p=0, r>0     q=0
    A_pqr_6 = A_pqr_comp_6 ./ (kc_sq - k_pqr_sq_6) ;
    % (7)  p>0, r=0     q=0
    A_pqr_7 = A_pqr_comp_7 ./ (kc_sq - k_pqr_sq_7) ;
    % (8)  p,q,r = 0
    A_pqr_8 = A_pqr_comp_8 / kc_sq  ;
    P_8 = A_pqr_8 * Pxcomp_8 * Pzcomp_8;

    for i = 1:patchamt
        %disp([loop i])
        [ix1,ix2,iz1,iz2] = findcoor(i,nx,delx,delz);
        ix1 = ix1+Ax;    ix2 = ix2+Ax;
        iz1 = iz1+Az;    iz2 = iz2+Az;
        
        % (1)  p,r > 0     q>0
        Pxcomp_1 = sin(p1pi_lx*ix2) - sin(p1pi_lx*ix1);
        Pzcomp_1 = sin(r1pi_lz*iz2) - sin(r1pi_lz*iz1);
        P_1 = A_pqr_1.* Pxcomp_1.* Pzcomp_1;
        
        % (2)  p=0, r>0     q>0
        Pzcomp_2 = sin(r2pi_lz*iz2) - sin(r2pi_lz*iz1);
        P_2 = A_pqr_2.* Pxcomp_2.* Pzcomp_2;
        
        % (3)  p>0, r=0     q>0
        Pxcomp_3 = sin(p3pi_lx*ix2) - sin(p3pi_lx*ix1);
        P_3 = A_pqr_3.* Pxcomp_3.* Pzcomp_3;
        
        % (4)  p,r = 0     q>0
        P_4 = A_pqr_4.* Pxcomp_4.* Pzcomp_4;
        
        % (5)  p,r > 0     q=0
        Pxcomp_5 = sin(p5pi_lx*ix2) - sin(p5pi_lx*ix1);
        Pzcomp_5 = sin(r5pi_lz*iz2) - sin(r5pi_lz*iz1);
        P_5 = A_pqr_5.* Pxcomp_5.* Pzcomp_5;
        
        % (6)  p=0, r>0     q=0
        Pzcomp_6 = sin(r6pi_lz*iz2) - sin(r6pi_lz*iz1);
        P_6 = A_pqr_6.* Pxcomp_6.* Pzcomp_6;
        
        % (7)  p>0, r=0     q=0
        Pxcomp_7 = sin(p7pi_lx*ix2) - sin(p7pi_lx*ix1);
        P_7 = A_pqr_7.* Pxcomp_7.* Pzcomp_7;
        
        P_sum = sum(P_1(:))+sum(P_2(:))+sum(P_3(:))+sum(P_4(:))+sum(P_5(:))+sum(P_6(:))+sum(P_7(:))+P_8 ;
        P(i,loop) = P_sum /delx/delz;
    end
    
end
bpp = P;

end

