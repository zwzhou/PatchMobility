clear;clc;
% VII. SOURCE ROOM MODELING --> B.Blocked patch pressure --> FIG.10
c_air = 343.6; % sound speed
S_0 = 2;
lx=11.5; ly=8.69; lz=4.03;% cavity dimension
Xs = 2; Ys = 4; Zs = 1;% source location
Lx = 1.5;   Lz = 0.96;
nx = 19;    nz = 13;
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

np = 18;   nq = 10;   nr = 8;

%N_pqr_n0 = lx*ly*lz /8;
freq = 50:500 ;
freqamt = length( freq );
T_r = 10;
for loop = 1:freqamt
    disp(freq(loop))
    eta_r = 2.2/(freq(loop)*T_r);
    cc = c_air* sqrt(1 + 1i*eta_r);
    omega = 2* pi* freq( loop );
    kc_sq = ( omega / cc )^2;
    P_sum = 0;

    for i = 1:patchamt
        %disp([loop i])
        [ix1,ix2,iz1,iz2] = findcoor(i,nx,delx,delz);
        ix1 = ix1+5.245;    ix2 = ix2+5.245;
        iz1 = iz1+1.27;     iz2 = iz2+1.27;
        P_sum = 0;
        for r = 0:nr
            if r==0
                Nzcomp = lz;
                Pzcomp = delz;
            else
                Nzcomp = lz/2;
                Pzcomp = (sin(r*iz2*pi_lz) - sin(r*iz1*pi_lz))* lz/(r*pi) ;
            end
            for p = 0:np
                if p==0
                    Nxcomp = lx;
                    Pxcomp = delx;
                else
                    Nxcomp = lx/2;
                    Pxcomp = (sin(p*ix2*pi_lx) - sin(p*ix1*pi_lx))* lx/(p*pi);
                end
                for q = 0:nq
                    if q==0
                        Nycomp = ly;
                    else
                        Nycomp = ly/2;
                    end
                    k_pqr_squared =  ( (p*pi_lx)^2 + (q*pi_ly)^2 + (r*pi_lz)^2 ); 
                    N_pqr = Nzcomp* Nycomp* Nxcomp;
                    A_pqr = S_0*cos(Xs*p*pi_lx)*cos(Ys*q*pi_ly)*cos(Zs*r*pi_lz) / ((kc_sq-k_pqr_squared)*N_pqr);
                    preasure = A_pqr* Pxcomp* Pzcomp;
                    P_sum = P_sum + preasure ;
                end
            end
        end
        P(i,loop) = P_sum;
    end
    
end
figure(1);
P1 = ones(1,freqamt);
P1(:) = P(1,:);
plot(freq,20*log10(abs(P1)),'b');
hold on;
P124 = ones(1,freqamt);
P124(:) = P(124,:);
plot(freq,20*log10(abs(P124)),'r');
legend('patch1','patch124')
title( 'Blocked Patch Pressure' );
xlabel('Frequency (HZ)');
ylabel('Pressure Level (dB)');
set(gca,'XScale','log')