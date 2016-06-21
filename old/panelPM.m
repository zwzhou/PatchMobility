clear;clc;
% panel patch mobility matrix
tic
lx = 1.5;
ly = 0.96;
p_area = lx*ly;
h = 0.002;
nx = 19;
ny = 13;
patchamt = nx*ny;
delx = lx/nx;
dely = ly/ny;
c_air = 343.6;% m/s
rho = 2700 ; % kg/m^3
E = 6.9e10 ; % Pa
eta = 0.01 ; % damping loss factor
mu = 0.334 ;  % Poisson ratio
Ec = E *( 1 + 1i*eta ); % complex Young's modulus
Dc = Ec *h^3 /( 12*(1-mu^2) ); 


patch(patchamt) = struct('x1',0,'x2',0,'y1',0,'y2',0,'xc',0,'yc',0,'area',0);
for i= 1:patchamt
    [x1,x2,y1,y2] = findcoor(i,nx,delx,dely);
    xc = (x1+x2) / 2;
    yc = (y1+y2) / 2;
    p_area = abs(x2-x1)*abs(y2-y1);
    patch(i) = struct('x1',x1,'x2',x2,'y1',y1,'y2',y2,'xc',xc,'yc',yc,'area',p_area);
end


freq = 10:1:100 ;
freqamt = length( freq );
np = 20;
nq = 18;

%-----------------------------------
% accelerate computational speed
[p,q] = meshgrid( 1:np, 1:nq );
two_int_comp = ( lx*ly./(p.*q*pi^2) ).^2 ;
p_pidivlx = p*pi/lx;
q_pidivly = q*pi/ly;
omega_pq_sq = pi^4*Dc/(rho*h)* ( (p/lx).^2 + (q/ly).^2 ).^2;
idivs = 1i/(delx*dely)^2 ;
%-----------------------------------
M_pq = rho*h*lx*ly/4 ;
YP0 = zeros(patchamt);
YP = zeros(patchamt,patchamt,freqamt);

for l = 1:freqamt
    fprintf('%4.2f\n',freq(l))
    omega = 2*pi*freq(l);
    YP_deno = M_pq*omega_pq_sq - M_pq*omega^2 ;
    ipo = idivs*omega ;
    %YPsum = 0 ;
    for i=1:patchamt
        intsi = (cos(p_pidivlx*patch(i).x2)-cos(p_pidivlx*patch(i).x1)) .* (cos(q_pidivly*patch(i).y2)-cos(q_pidivly*patch(i).y1));
        YP_mid = two_int_comp.*intsi./YP_deno ;
        for j=i:patchamt
            intsj = (cos(p_pidivlx*patch(j).x2)-cos(p_pidivlx*patch(j).x1)) .* (cos(q_pidivly*patch(j).y2)-cos(q_pidivly*patch(j).y1));
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

toc

figure(1)
tr_patchmobility = zeros(1,freqamt);
tr_patchmobility(:) = YP(42,72,:);
plot(freq,20*log10(abs(tr_patchmobility)),'b')
hold on
ip_patchmobility = zeros(1,freqamt);
ip_patchmobility(:) = YP(42,42,:);
plot(freq,20*log10(abs(ip_patchmobility)),'m')
legend('transfer patch mobility','input patch mobility')
xlabel('Frequency (HZ)');
ylabel('Magnitude (dB)');
