clear
% VII. SOURCE ROOM MODELING --> A.Quadratic room pressure --> FIG.9
% disp('Calculating Quadratic Room Pressure')
c = 343.6;  % m/s
rho_a = 1.293 ; % kg/m^3; air density
S_0 = 2; % Source amplitude
%
%
lx = 11.5;ly = 8.69;lz = 4.03; % Cavity size
Xs = 2;Ys = 4;Zs = 1; % Source location
% pi_lx=pi/lx;pi_ly=pi/ly;pi_lz=pi/lz;
%

freq = 50:0.5:1000 ; 
freqNUM = length( freq );

% fprintf('step %3d/3  Modal Frequency\n',1)

np = 350;
nq = 250;
nr = 150;
Nmodal = 100000 ;

% Quadratic Room Pressure
Pr_squa = f_quadPr(lx,ly,lz,Xs,Ys,Zs,S_0,c,freq,rho_a,Nmodal,np,nq,nr );

% ------------------------------------------------------------------------------------------------
% % [omegasqua,idx_c] = f_cavityfreqsqua(lx,ly,lz,c,Nmodal,np,nq,nr);
% [kpqr_squa,idx_c] = f_cavityksqua(lx,ly,lz,Nmodal,np,nq,nr);
% 
% fprintf('step %3d/3  Apqr\n',2)
% 
% % intphis = repmat(cos(idx_c.x*pi*Xs/lx).*cos(idx_c.y*pi*Ys/ly).*cos(idx_c.z*pi*Zs/lz),1,freqNUM);
% intphis = S_0 * cos(idx_c.x*pi*Xs/lx).*cos(idx_c.y*pi*Ys/ly).*cos(idx_c.z*pi*Zs/lz);
% 
% lxcomp = lx*ones(Nmodal,1);
% lxcomp(idx_c.x~=0) = 0.5*lx;
% lycomp = ly*ones(Nmodal,1);
% lycomp(idx_c.y~=0) = 0.5*ly;
% lzcomp = lz*ones(Nmodal,1);
% lzcomp(idx_c.z~=0) = 0.5*lz;
% % Npqr = repmat(lxcomp.*lycomp.*lzcomp,1,freqNUM);
% Npqr = lxcomp.*lycomp.*lzcomp;
% 
% % k_pqr_squa = repmat(omegasqua/c^2,1,freqNUM);
% % k_pqr_squa = omegasqua/c^2;
% 
% T_r = 10;
% eta_r = 2.2./(freq*T_r);
% cc = c* sqrt(1 + 1i*eta_r);
% % kc_squa = repmat( (2*pi*freq./cc).^2, Nmodal,1);
% kc_squa = (2*pi*freq./cc).^2;
% 
% % A_pqr size: Nmodal*freqNUM
% % A_pqr = bsxfun( @rdivide,intphis./Npqr,(bsxfun(@minus,kc_squa,kpqr_squa)) );
% omega = 2*pi*freq;
% A_pqr = bsxfun(@times,1i*omega*rho_a,intphis./Npqr)./bsxfun(@minus,kc_squa,kpqr_squa) ;
% 
% fprintf('step %3d/3  Pr^2\n',3)
% 
% epsnp = ones(1,Nmodal);
% epsnq = ones(1,Nmodal);
% epsnr = ones(1,Nmodal);
% 
% % epsnp(idx_c.x~=0) = 2;
% % epsnq(idx_c.y~=0) = 2;
% % epsnr(idx_c.z~=0) = 2;
% 
% epsnp(idx_c.x==0) = 2;
% epsnq(idx_c.y==0) = 2;
% epsnr(idx_c.z==0) = 2;
% 
% 
% % size: 1*freqNUM
% Pr_squa = (epsnp.*epsnq.*epsnr) * (abs(A_pqr).^2) ./8;
% 
% disp('Complete')
% ---------------------------------------------------------------------------------------------



% %vs = (lx*ly*lz * S_0 / N_pqr )^2 ;
% P_r_squared = ones( 1,freqNUM );
% for i = 1:freqNUM
%     fprintf('%4.2f\n',freq(i))
%     T_r = 10;
%     eta_r = 2.2/(freq(i)*T_r);
%     cc = c* sqrt(1 + 1i*eta_r);
%     omega = 2* pi* freq( i );
%     kc_squa = (omega / cc)^2;
%     sig_sum = 0;
%     
%     for r = 0:nr
%         if r==0
%             epsilon_r = 1;
%             Nzcomp = lz;
%         else
%             epsilon_r = 2;
%             Nzcomp = lz/2;
%         end
%         A_pqr_zcomp = cos(Zs*r*pi_lz);
%         for q = 0:nq
%             if q==0
%                 epsilon_q = 1;
%                 Nycomp = ly;
%             else
%                 epsilon_q = 2;
%                 Nycomp = ly/2;
%             end
%             A_pqr_ycomp = cos(Ys*q*pi_ly);
%             for p = 0:np
%                 if p==0
%                     epsilon_p = 1;
%                     Nxcomp = lx;
%                 else
%                     epsilon_p = 2;
%                     Nxcomp = lx/2;
%                 end
%                 N_pqr = Nxcomp*Nycomp*Nzcomp;
%                 k_pqr_sq =  ( (p*pi_lx)^2 + (q*pi_ly)^2 + (r*pi_lz)^2 ); 
%                 A_pqr = S_0* cos(Xs*p*pi_lx)*A_pqr_ycomp*A_pqr_zcomp / ((kc_sq-k_pqr_sq)*N_pqr);
%                 Pr2_tmp = abs(A_pqr^2) * epsilon_p*epsilon_q*epsilon_r /8 ;
%                 % tmp = (abs( cos(Xs*p*pi/lx)*cos(Ys*q*pi/ly)*cos(Zs*r*pi/lz)/(kc^2-k_pqr_squared) ))^2 * epsilon_p*epsilon_q*epsilon_r ;
%                 sig_sum = sig_sum + Pr2_tmp ;
%             end
%         end
%     end
% 
%     P_r_squared(i) = sqrt(sig_sum) ;%* vs / 8 ;
% end
figure(4)
plot(freq,20*log10((Pr_squa)/2e-5))
title( 'Quadratic Room Pressure' );
xlabel('Frequency (HZ)');
ylabel('Pressure Level (dB)');
set(gca,'XScale','log')
xlim([freq(1) freq(freqNUM)]);
