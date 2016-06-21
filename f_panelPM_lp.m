function [ Ypanel ] = f_panelPM_lp( lx,ly,h, nx,ny, rho,Dc, freq, Nmodal,np,nq )
% return patch mobility martix of a panel: patchamt*patchamt*freqamt
%   lx, ly & h: length, width and thickness of the plate
%   nx & ny: patch number along length and width
%   freq: frequency vector. Example: freq=1:0.1:100
% 
disp('Calculating Panel Patch Mobility (loop)')
%lx = 1.5;
%ly = 0.96;
%p_area = lx*ly;
%h = 0.002;
% nx = 19;
% ny = 13;
patchAmt = nx*ny;
delx = lx/nx;
dely = ly/ny;
%c_air = 343.6;% m/s
% rho = 2700 ; % kg/m^3
% E = 6.9e10 ; % Pa
% eta = 0.01 ; % damping loss factor
% mu = 0.334 ;  % Poisson ratio
% Ec = E *( 1 + 1i*eta ); % complex Young's modulus
% Dc = Ec *h^3 /( 12*(1-mu^2) ); 


% patch location
% patch(patchAmt) = struct('x1',0,'x2',0,'y1',0,'y2',0,'xc',0,'yc',0,'area',0);
% for i= 1:patchAmt
%     [x1,x2,y1,y2] = findcoor(i,nx,delx,dely);
%     xc = (x1+x2) / 2;
%     yc = (y1+y2) / 2;
%     p_area = abs(x2-x1)*abs(y2-y1);
%     patch(i) = struct('x1',x1,'x2',x2,'y1',y1,'y2',y2,'xc',xc,'yc',yc,'area',p_area);
% end
patchnumber = 1:patchAmt;
patch = f_findcoor(patchnumber,nx,delx,dely);% size: 1*patchAmt

%freq = 10:0.5:600 ;
freqamt = length( freq );
% Nmodal = 500;
% np = 20;
% nq = 18;


% modal frequency
% [omega_c_squa,idx_p]  size: Nmodal*1
[omega_pq_squa,idx_p] = f_plateOmegaSqua_simspt(lx,ly,Dc,rho*h,Nmodal,np,nq);
%-----------------------------------
% accelerate computational speed
% [idx_p.x,idx_p.y] = meshgrid( 1:np, 1:nq );
two_int_comp = ( lx*ly./(idx_p.x.*idx_p.y*pi^2) ).^2 ;%size: Nmodal*1
p_pidivlx = idx_p.x*pi/lx;
q_pidivly = idx_p.y*pi/ly;
% omega_pq_squa = pi^4*Dc/(rho*h)* ( (idx_p.x/lx).^2 + (idx_p.y/ly).^2 ).^2;
idivs = 1i/(delx*dely)^2 ;
%-----------------------------------
M_pq = rho*h*lx*ly/4 ;
YP0 = zeros(patchAmt);
YP = zeros(patchAmt,patchAmt,freqamt);

for l = 1:freqamt
    fprintf('%4.2f\n',freq(l))
    omega = 2*pi*freq(l);
    YP_deno = M_pq*omega_pq_squa - M_pq*omega^2 ;
    ipo = idivs*omega ;
    %YPsum = 0 ;
    for i=1:patchAmt
        intsi = (cos(p_pidivlx*patch.x2(i))-cos(p_pidivlx*patch.x1(i))) .*...
            (cos(q_pidivly*patch.y2(i))-cos(q_pidivly*patch.y1(i)));
        YP_mid = two_int_comp.*intsi./YP_deno ;
        for j=i:patchAmt
            intsj = (cos(p_pidivlx*patch.x2(j))-cos(p_pidivlx*patch.x1(j))) .*...
                (cos(q_pidivly*patch.y2(j))-cos(q_pidivly*patch.y1(j)));
            YP_comp = YP_mid .* intsj ;
%             YP0(i,j) = 1i*omega/(patch(i).area*patch(j).area)* sum(YP_comp(:)) ; 
            YP0(i,j) = ipo* sum(YP_comp(:)) ; 
            if (i~=j)
                YP0(j,i) = YP0(i,j) ;
            end
        end
    end
    YP(:,:,l) = YP0 ;
end


Ypanel = YP;


end

