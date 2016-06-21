clear;clc;
Lx = 1.5 ;   % m
Ly = 0.96 ;  % m
rho = 2373 ; % kg/m^3
h = 0.002 ;  % m
E = 6e10 ;   % Pa
eta = 0.01 ; % damping loss factor
Ec = E *( 1 + 1i*eta ); % complex Young's modulus
mu = 0.22 ;  % Poisson ratio
Dc = Ec *h^3 /( 12*(1-mu^2) ); 
np = 19;
nq = 13;
delS = 0.079*0.074; % m^2;


A42x = 0.24;%0.237    ; % m  ;
A42y = 0.15;%0.148    ;
B42x = 0.32;%0.316    ;
B42y = 0.22;%0.222    ;
A72x = 1.11;%1.106    ;
A72y = 0.22;%0.222    ;
B72x = 1.18;%1.185    ;
B72y = 0.30;%0.296    ;

tic;
M_pq = Lx *Ly *rho *h /4;

%-----------------------------------------------------------------------------------------
int_42 = int_i_S( np, nq, Lx, Ly, A42x, A42y, B42x, B42y );
int_72 = int_i_S( np, nq, Lx, Ly, A72x, A72y, B72x, B72y );
[ p, q ] = meshgrid( 1:np, 1:nq );
Kc_pq = ( p.^4*Ly/Lx^3 + q.^4*Lx/Ly^3 + 2*p.^2.*q.^2/(Lx*Ly) )* Dc*pi^4 /4;
omegac_pq_square = Dc* (pi^4) / (rho*h) *( (p/Lx).^2 + (q/Ly).^2 ).^2 ;%Kc_pq / M_pq ;
mid_a = int_42 .* int_72 / M_pq ;
mid_b = int_42 .* int_42 / M_pq ;
%------------------------------------------------------------------------------------------



freq = 1:1:600;
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


%     y_sum = 0;
%     for q = 1:nq
%         for p = 1:np
%             omegac_pq_square = Dc* (pi^4) / (rho*h) *( (p/Lx)^2 + (q/Ly)^2 )^2;
%             % Kc_pq = Dc*pi^4 *( Ly*p^4/Lx^3 +Lx*q^4/Ly^3 + 2*p^2*q^2/(Lx*Ly) ) /4 ;
%             y_sum = y_sum + ...
%                     Lx*Ly/(p*q*pi^2)*( cos(p*pi*B42x/Lx)-cos(p*pi*A42x/Lx) )*( cos(q*pi*B42y/Ly)-cos(q*pi*A42y/Ly) )*...
%                     Lx*Ly/(p*q*pi^2)*( cos(p*pi*B72x/Lx)-cos(p*pi*A72x/Lx) )*( cos(q*pi*B72y/Ly)-cos(q*pi*A72y/Ly) )/...
%                     ( M_pq * ( omegac_pq_square - omega^2 ) ) ;
%                     %( M_pq *( Kc_pq/M_pq - omega^2) ) ; 
%         end
%     end
%     yp_t( ii ) = 1i * omega * y_sum / (delS*delS) ;   


 end
 
toc
 
tr_db_t = 20* log10( abs(yp_t) );
tr_db_i = 20* log10( abs(yp_i) );
plot(freq,tr_db_t,'g')
hold on
plot(freq,tr_db_i,'b')
