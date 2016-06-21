clear;clc;
% IV. PANEL PATCH MOBILITY  --> FIG.3
tic;
%profile on;
Lx = 1.50 ;  % m
Ly = 0.96 ;  % m
sz = [ Lx, Ly ];
h = 0.002 ;  % m
rho = 2700 ; % kg/m^3
E = 6.9e10 ; % Pa
eta = 0.01 ; % damping loss factor
mu = 0.334 ;  % Poisson ratio
Ec = E *( 1 + 1i*eta ); % complex Young's modulus
Dc = Ec *h^3 /( 12*(1-mu^2) ); 

np = 19;
nq = 13;
pq = [np, nq];
%
%
patch_42 = [ 0.237  0.316  0.148  0.222   ]; % m  [ x1 x2  y1 y2  z ] 
patch_72 = [ 1.106  1.185  0.222  0.296   ];

delS = 0.079*0.074; % m^2;
%
%
M_pq = Lx *Ly *rho *h /4;
% Kc_pq = ( p.^4*Ly/Lx^3 + q.^4*Lx/Ly^3 + 2*p.^2.*q.^2/(Lx*Ly) )* Dc*pi^4 /4;
%
%-----------------------------------------------------------------------------------------
int_42 = int_i_S( pq, sz, patch_42 );
int_72 = int_i_S( pq, sz, patch_72 );
[ p, q ] = meshgrid( 1:np, 1:nq );
omegac_pq_square = Dc* (pi^4) / (rho*h) *( (p/Lx).^2 + (q/Ly).^2 ).^2 ;%Kc_pq / M_pq ;
mid_a = int_42 .* int_72 / M_pq ;
mid_b = int_42 .* int_42 / M_pq ;
%------------------------------------------------------------------------------------------
%
%
%
freq = 10:0.5:600;
freqnum = length( freq );
yp_t = ones( 1, freqnum );
yp_i = ones( 1, freqnum );
 for ii = 1:freqnum
	omega = 2 *pi *freq( ii );

%-----------------------------------------------------------------------------------------
    mid = 1i*omega/delS^2 ;
    mid_t = mid_a ./ (omegac_pq_square - omega^2) ;
    yp_t(ii) = mid * sum( mid_t(:) );
    mid_i = mid_b ./ (omegac_pq_square - omega^2) ;
    yp_i(ii) = mid * sum( mid_i(:) );
%-----------------------------------------------------------------------------------------
% 
%
% ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
%     y_sum = 0;
%     for q = 1:nq
%         for p = 1:np
%             omegac_pq_square = Dc* (pi^4) / (rho*h) *( (p/Lx)^2 + (q/Ly)^2 )^2;
%             % Kc_pq = Dc*pi^4 *( Ly*p^4/Lx^3 +Lx*q^4/Ly^3 + 2*p^2*q^2/(Lx*Ly) ) /4 ;
%             y_sum = y_sum + ...
%                     Lx*Ly/(p*q*pi^2)*( cos(p*pi*patch_42(2)/Lx)-cos(p*pi*patch_42(1)/Lx) )*( cos(q*pi*patch_42(4)/Ly)-cos(q*pi*patch_42(3)/Ly) )*...
%                     Lx*Ly/(p*q*pi^2)*( cos(p*pi*patch_72(2)/Lx)-cos(p*pi*patch_72(1)/Lx) )*( cos(q*pi*patch_72(4)/Ly)-cos(q*pi*patch_72(3)/Ly) )/...
%                     ( M_pq * ( omegac_pq_square - omega^2 ) ) ;
%                     %( M_pq *( Kc_pq/M_pq - omega^2) ) ; 
%         end
%     end
%     yp_t( ii ) = 1i * omega * y_sum / (delS*delS) ;   
% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% 
 end
 
toc
 
%profile viewer
figure(2)
tr_db_t = 20* log10( abs(yp_t) );   % transfer
tr_db_i = 20* log10( abs(yp_i) );   % input
plot(freq,tr_db_t,'b')
hold on
plot(freq,tr_db_i,'m')
legend('transfer patch mobility','input patch mobility')
xlabel('Frequency (HZ)');
ylabel('Magnitude (dB)');
hold off
