function [ Zrad,dis ] = f_radZ( lx,ly, nx,ny, rho_a,c,freq )
% radiation impedance matrix Zrad, Zrad=Yrad^-1; patchamt*patchamt*freqNum
% paprer function (39)
%   see: test_radZ.m
% Zrad size: patchamt*pathamt*freqNum



% rho = 1.293 ; % kg/m^3
% c = 340  ;  % m/s

% lx = 1.5; ly = 0.96; % panel size
% nx = 19;  ny = 13; % x y patch count
patchAmt = nx * ny ;
% delx = lx / nx;  dely = ly / ny; % patch size
S_j = lx*ly / patchAmt; % patch aera
ai = 1.0*abs( sqrt(S_j/pi) ); % equivalent patch radius


% initialization
dis = f_distance( lx, ly, nx, ny );
% patchnumber = 1:patchAmt;
% patch = f_findcoor(patchnumber,nx,delx,dely);
% xdis = bsxfun(@minus,patch.x2.',patch.x2);% patch x distance
% ydis = bsxfun(@minus,patch.y2.',patch.y2);% patch y distance
% dis = sqrt( xdis.^2 + ydis.^2 );% patch distance

% for i = 1:xn
%     for j = 1:yn
%         ptn = j+(i-1)*yn;   % patch number
%         lbx = (i-1)*delta_x;
%         lby = j*delta_y;    % lbx lby: left-bottom;
%         rtx = i*delta_x;
%         rty = (j-1)*delta_y; %  rtx rty: right-top;
%         cx = (lbx + rtx) / 2;
%         cy = (lby + rty) / 2; %  cx cy: center;
%         patch( ptn ) = struct( 'lbx',lbx,'lby',lby,'rtx',rtx,'rty',rty,'cx',cx,'cy',cy );
%     end
% end


% patch( patchamt ) = struct( 'x1',0,'x2',0,'y1',0,'y2',0,'xc',0,'yc',0 ); % initialization
% for i=1:patchamt
%     [x1,x2,y1,y2] = findcoor(i,nx,delx,dely);
%     xc = (x1+x2) / 2;
%     yc = (y1+y2) / 2;
%     patch(i) = struct( 'x1',x1,'x2',x2,'y1',y1,'y2',y2,'xc',xc,'yc',yc );
% end


% dis = zeros( patchamt ); % distance between patch i and patch j
% for i = 1 : patchamt-1
%     for j = i+1 : patchamt
%         dis(i,j) = abs( sqrt( (patch(i).xc-patch(j).xc)^2 + (patch(i).yc-patch(j).yc)^2 ) );
%         dis(j,i) = dis(i,j);
%     end
% end

% freq = 1:1:1000;

omega = 2*pi*freq;
k = omega / c;


freqNum = length( freq );

if freqNum>1
    disp('Calculating Radiation Impedance Matrix Z')
    
    Z = zeros( patchAmt ); % initialization
    Zrad = ones( patchAmt, patchAmt, freqNum );


    fprintf(' ---> current: %6.1f / ( start: %4.1f   step: %3.1f   end: %5.1f )\n',...
        freq(1), freq(1), freq(2)-freq(1), freq(freqNum) )% didplay
    
    for ii = 1:freqNum
        
        if ~mod(ii,100)
            fprintf(' ---> current: %6.1f / ( start: %4.1f   step: %3.1f   end: %5.1f )\n',...
                freq(ii), freq(1), freq(2)-freq(1), freq(freqNum) )% didplay
        end
    
%        omega = 2*pi*freq(l);
%        k = omega / c;
%        Z(dis==0) = rho*c*(1-exp(-1i*k*ai));
%        Z(dis~=0) = omega*rho*1i/2/pi*S_j* exp(-1i*k*dis(dis~=0)) ./ dis(dis~=0);
    
        Z(dis==0) = rho_a*c*(1-exp(-1i*k(ii)*ai));
        Z(dis~=0) = omega(ii)*rho_a*1i/2/pi*S_j* exp(-1i*k(ii)*dis(dis~=0)) ./ dis(dis~=0);
    
        Zrad(:,:,ii) = Z;
   
    %     YP(:,:,ii) = inv(Z);
    end
    disp('Complete!')
    
else
    Zrad = zeros( patchAmt ); % initialization
    Zrad(dis==0) = rho_a*c*(1-exp(-1i*k*ai));
    Zrad(dis~=0) = omega*rho_a*1i/2/pi*S_j* exp(-1i*k*dis(dis~=0)) ./ dis(dis~=0);

end


end

