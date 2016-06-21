clear;clc;
% VII. SOURCE ROOM MODELING --> B.Blocked patch pressure --> FIG.10
c = 340; % sound speed
S_0 = 2;
lx=11.5; ly=8.69; lz=4.03;% cavity dimension
Xs = 2; Ys = 4; Zs = 1;% source location
x124=6; x1=5.245; z124=1.75; z1=1.27; % patch 124 and 1 location
delta_x = 0.08; delta_z = 0.074; % patch size
p124_x0=x124-delta_x/2; p124_x1=x124+delta_x/2; p1_x0=x1-delta_x/2; p1_x1=x1+delta_x/2; % integral from x0 to x1
p124_z0=z124-delta_z/2; p124_z1=z124+delta_z/2; p1_z0=z1-delta_z/2; p1_z1=z1+delta_z/2; % integral from z0 to z1
np = 144; nq = 117; nr = 55;
N_pqr = lx*ly*lz /8;
freq = 51:1:100 ;
freqNUM = length( freq );
P_r_124 = ones( 1,freqNUM );
P_r_1 = ones( 1,freqNUM );
for i = 1:freqNUM
    disp(i)
    T_r = 10;
    eta_r = 2.2/(freq(i)*T_r);
    cc = c* sqrt(1 + 1i*eta_r);
    omega = 2* pi* freq( i );
    kc = omega / cc;
    p124_sum = 0;
    p1_sum = 0;
    for r = 1:nr
        for q = 1:nq
            for p = 1:np
                k_pqr_squared =  ( (pi*p/lx)^2 + (pi*q/ly)^2 + (pi*r/lz)^2 ); % * ( c/(2*pi) )^2 ;
                A_pqr = cos(Xs*p*pi/lx)*cos(Ys*q*pi/ly)*cos(Zs*r*pi/lz)*S_0 / ((kc^2-k_pqr_squared)*N_pqr);
                psi_tmp = lx*lz/(p*r*pi^2);
                i_Psi_pqr124 = psi_tmp * (sin(p*pi*p124_x1/lx)-sin(p*pi*p124_x0/lx))*(sin(p*pi*p124_z1/lz)-sin(r*pi*p124_z0/lz));
                i_Psi_pqr1 = psi_tmp * (sin(p*pi*p1_x1/lx)-sin(p*pi*p1_x0/lx))*(sin(p*pi*p1_z1/lz)-sin(r*pi*p1_z0/lz));
                p124_sum = p124_sum + A_pqr * i_Psi_pqr124 ;
                p1_sum = p1_sum + A_pqr * i_Psi_pqr1 ;
            end
        end
    end
    P_r_124(i) = p124_sum ;%* vs / 8 ;
    P_r_1(i) = p1_sum ;
end
figure(1);
plot(freq,20*log10(abs(P_r_124)/2e-5),'b');
hold on;
plot(freq,20*log10(abs(P_r_1)/2e-5),'r');
legend('patch124','patch1')
title( 'Blocked Patch Pressure' );
xlabel('Frequency (HZ)');
ylabel('Pressure Level (dB)');