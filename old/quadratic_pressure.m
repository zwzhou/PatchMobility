clear;clc;
% VII. SOURCE ROOM MODELING --> A.Quadratic room pressure --> FIG.9
c = 343.6;  % m/s
S_0 = 2; % Source amplitude
%
%
lx = 11.5;ly = 8.69;lz = 4.03; % Cavity size
Xs = 2;Ys = 4;Zs = 1; % Source location
pi_lx=pi/lx;pi_ly=pi/ly;pi_lz=pi/lz;
%
np = 70;nq = 55;nr = 45;
%
%
freq = 50:1:100 ; 
freqNUM = length( freq );

%vs = (lx*ly*lz * S_0 / N_pqr )^2 ;
P_r_squared = ones( 1,freqNUM );
for i = 1:freqNUM
    fprintf('%4.2f\n',freq(i))
    T_r = 10;
    eta_r = 2.2/(freq(i)*T_r);
    cc = c* sqrt(1 + 1i*eta_r);
    omega = 2* pi* freq( i );
    kc_sq = (omega / cc)^2;
    sig_sum = 0;
    for r = 0:nr
        if r==0
            epsilon_r = 1;
            Nzcomp = lz;
        else
            epsilon_r = 2;
            Nzcomp = lz/2;
        end
        A_pqr_zcomp = cos(Zs*r*pi_lz);
        for q = 0:nq
            if q==0
                epsilon_q = 1;
                Nycomp = ly;
            else
                epsilon_q = 2;
                Nycomp = ly/2;
            end
            A_pqr_ycomp = cos(Ys*q*pi_ly);
            for p = 0:np
                if p==0
                    epsilon_p = 1;
                    Nxcomp = lx;
                else
                    epsilon_p = 2;
                    Nxcomp = lx/2;
                end
                N_pqr = Nxcomp*Nycomp*Nzcomp;
                k_pqr_sq =  ( (p*pi_lx)^2 + (q*pi_ly)^2 + (r*pi_lz)^2 ); 
                A_pqr = S_0* cos(Xs*p*pi_lx)*A_pqr_ycomp*A_pqr_zcomp / ((kc_sq-k_pqr_sq)*N_pqr);
                Pr2_tmp = abs(A_pqr^2) * epsilon_p*epsilon_q*epsilon_r /8 ;
                % tmp = (abs( cos(Xs*p*pi/lx)*cos(Ys*q*pi/ly)*cos(Zs*r*pi/lz)/(kc^2-k_pqr_squared) ))^2 * epsilon_p*epsilon_q*epsilon_r ;
                sig_sum = sig_sum + Pr2_tmp ;
            end
        end
    end
    P_r_squared(i) = sqrt(sig_sum) ;%* vs / 8 ;
end
figure(2)
plot(freq,20*log10(abs(P_r_squared)))
title( 'Quadratic Room Pressure' );
xlabel('Frequency (HZ)');
ylabel('Pressure Level (dB)');
set(gca,'XScale','log')
