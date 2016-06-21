clear;clc;
Lx = 1500;  % mm
Ly = 960;   % mm
rho = 2.7e-3;   % g/mm^3
h = 2;      % mm
E = 6e4;    % MPa
eta = 0.5 ; % damping loss factor
Ec = E *( 1 + 1i*eta ); % complex Young's modulus
mu = 0.22;  % Poisson ratio
Dc = Ec *h^3 /( 12*(1-mu^2) );
np = 19;
nq = 13;
delS = 79*74;


M_pq = Lx *Ly *rho *h /4;


int_42 = int_i_S( np, nq, Lx, Ly,  237, 148,  316, 222 );
int_72 = int_i_S( np, nq, Lx, Ly, 1106, 222, 1185, 296 );
[ p, q ] = meshgrid( 1:np, 1:nq );
Kc_pq = ( p.^4*Ly/Lx^3 + q.^4*Lx/Ly^3 + 2*p.^2.*q.^2/(Lx*Ly) )* Dc*pi^4 /4;
omegac_pq = sqrt( Kc_pq / M_pq );


freq = 1:10:601;
freqnum = length( freq );
yp_t = ones( 1, freqnum );
yp_i = ones( 1, freqnum );
for ii = 1:freqnum
    omega = 2 *pi *freq( ii );
    mid = -1i*omega/delS^2 ;
    mid_t = int_42 .* int_72 ./ (omegac_pq - omega^2) / M_pq ;
    yp_t(ii) = mid * sum( mid_t(:) );
    mid_i = int_42 .* int_42 ./ (omegac_pq - omega^2) * M_pq;
    yp_i(ii) = mid * sum( mid_i(:) );
end

plot(freq,yp_i)


% for q = 1:nq
%     for p = 1:np
%         ok
%         
%         
%         
%     end
% end