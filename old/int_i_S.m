function intiS = int_i_S( od, sz, patch )
% od: order of two or three dimision i.e. np, nq & nr;   sz: patch size
% patch: patch location
%
% ( Ax Ay ) : left-bottom coordinate values
% ( Bx By ) : right-top coordinate values
%
tmp = length( patch );
switch tmp
    case 4
% 
% if  patch = [ Ax Bx Ay By ]
%
% solve the equation below:
%-----------------------------------------------------
%               _Bx      P*pi      _By      q*pi
%     int_i_S = |   sin(----- x)dx |   sin(----- y)dy
%               `Ax       Lx       `Ay       Lx
%-----------------------------------------------------
%
%
        np = od(1);
        nq = od(2);
        Lx = sz(1);
        Ly = sz(2);
        [ p, q ] = meshgrid( 1:np, 1:nq );
        int_dx = ( - cos( p.*pi*patch(2)/Lx ) + cos( p.*pi*patch(1)/Lx ) ) ./ p ;
        int_dy = ( - cos( q.*pi*patch(4)/Ly ) + cos( q.*pi*patch(3)/Ly ) ) ./ q ;
        intiS = int_dx .* int_dy *Lx*Ly / pi^2 ;
    
    case 5
% 
% if  patch = [ Ax Bx Ay By z ]
%
% solve the equation below:
%-----------------------------------------------------------------
%               _Bx      P*pi      _By      q*pi           r*pi
%     int_i_S = |   cos(----- x)dx |   cos(----- y)dy*cos(----- z)
%               `Ax       Lx       `Ay       Lx              Lz
%-----------------------------------------------------------------
%
%
        np = od(1);
        nq = od(2);
        nr = od(3);
        Lx = sz(1);
        Ly = sz(2);
        Lz = sz(3);
        [ p, q, r ] = meshgrid( 1:np, 1:nq, 1:nr );
        int_dx = ( sin( p.*pi*patch(2)/Lx ) - sin( p.*pi*patch(1)/Lx ) ) ./ p ;
        int_dy = ( sin( q.*pi*patch(4)/Ly ) - sin( q.*pi*patch(3)/Ly ) ) ./ q ;
        intiS = int_dx .* int_dy *Lx*Ly / pi^2 .* cos( r.*pi*patch(5)/Lz ) ;
%         p(:,:,1)
end