close all;clc;clear;
% 未完成
Lx = 1.5; Ly = 0.96; Lz = 0.01;
n_x = 19 ; n_y = 13 ;% patch number of x and y
patch_amt = n_x * n_y ;%  amount of patch
del_x = Lx / n_x ;  del_y = Ly / n_y ; % patch dimension
rho = 1.293 ; % kg/m^3
c_air=343.6; %m/s, speed of sound in air 20
eta = 0.01;
cc = c_air* sqrt(1+1i*eta);
%
% if matlabpool('size')<=0 %判断并行计算环境是否已然启动
% matlabpool('open','local',7); %若尚未启动，则启动并行环境
% else
% disp('Already initialized'); %说明并行环境已经启动。
% end
% matlabpool local 4
tic
% l_patch( patch_amt ) = struct( 'x1',0,'y1',0,'x2',0,'y2',0,'xc',0,'yc',0,'z',0 ); % initialization
% r_patch( patch_amt ) = struct( 'x1',0,'y1',0,'x2',0,'y2',0,'xc',0,'yc',0,'z',0 );
% for i = 1:n_x
%     for j = 1:n_y
%         ptn = j+(i-1)*n_y;   % patch number
%         lbx = (j-1)*del_x;
%         lby = (i-1)*del_y;    % lbx lby: left-bottom;
%         rtx = j*del_x;
%         rty = i*del_y; %  rtx rty: right-top;
%         cx = (lbx + rtx) / 2;
%         cy = (lby + rty) / 2; %  cx cy: center;
%         l_patch( ptn ) = struct( 'x1',lbx,'y1',lby,'x2',rtx,'y2',rty,'xc',cx,'yc',cy,'z',0 );
%         r_patch( ptn ) = struct( 'x1',lbx,'y1',lby,'x2',rtx,'y2',rty,'xc',cx,'yc',cy,'z',Lz );
%     end
% end

np = 30;    nq = 30;    nr = 12; % np nq & nr must >= 2 !

freq = 100:101 ;
freqamt = length( freq );

patch_z0 = 0;
pqtch_z1 = Lz;
Z = ones( patch_amt, patch_amt, freqamt );
Y = ones( patch_amt, patch_amt, freqamt );
[p,q,r] = meshgrid(0:np,0:nq,0:nr);
pi_lx = pi/Lx; pi_ly = pi/Ly; pi_lz = pi/Lz;
k_pqr_sq = (p*pi_lx).^2 + (q*pi_ly).^2 + (r*pi_lz).^2 ;
N_pqr = ones(np+1,nq+1,nr+1)*Lx*Ly*Lz/8 ;
for m=1:np+1
    for n=1:nq+1
        for l=1:nr+1
            if m==1
                N_pqr(m,n,l) = N_pqr(m,n,l)*2;
            end
            if n==1
                N_pqr(m,n,l) = N_pqr(m,n,l)*2;
            end
            if l==1
                N_pqr(m,n,l) = N_pqr(m,n,l)*2;
            end
        end
    end
end

intsi = cos(r*pi_lz*patch_z0) ; % initialization
for m=1:np+1
    for n=1:nq+1
        for l=1:nr+1
            if m==1
                intsi(m,n,l) = intsi(m,n,l)* del_x ;
            end
            if n==1
                intsi(m,n,l) = intsi(m,n,l)* del_y ;
            end
        end
    end
end
for n=2:nq+1
    for r=1:nr+1
        intsi(1,n,r) = intsi(1,n,r) *Ly/(n-1)/pi;
    end
end
for m=2:np+1
    for r=1:nr+1
        intsi(m,1,r) = intsi(m,1,r) *Lx/(m-1)/pi;
    end
end
intsj = intsi ;
intsj0=intsi(2:np+1,2:nq+1,:);
intsi0=intsi(2:np+1,2:nq+1,:);
[p1,q1,r1] = meshgrid(1:np,1:nq,0:nr);
ipo = (Lx*Ly./(p1.*q1*pi^2)).^2 ;
p_pidivlx = p1*pi/Lx;
q_pidivly = q1*pi/Ly;
r_pidivlz = r1*pi/Lz;

for l = 1:freqamt

    Z0 = ones( patch_amt );
    %disp(l)
    omega = 2*pi*freq(l);
    kc_sq = (omega / cc)^2;
    for i = 1:patch_amt
        disp([freq(l),i])
        [ix1,ix2,iy1,iy2] = findcoor(i,n_x,del_x,del_y);
        int_xcomp1 = ( sin(p_pidivlx*ix2)-sin(p_pidivlx*ix1) ) ;
        int_ycomp1 = ( sin(q_pidivly*iy2)-sin(q_pidivly*iy1) ) ;
        intsi1 = intsi0.* int_xcomp1.* int_ycomp1;
        intsi(2:np+1,2:nq+1,:) = ipo.* intsi1;

        for j = i:patch_amt
            [jx1,jx2,jy1,jy2] = findcoor(j,n_x,del_x,del_y);
          
            % p,q,r > 0
            int_xcomp2 = ( sin(p_pidivlx*jx2)-sin(p_pidivlx*jx1) ) ;
            int_ycomp2 = ( sin(q_pidivly*jy2)-sin(q_pidivly*jy1) ) ;
            intsj1 = intsj0.* int_xcomp2.* int_ycomp2;
            intsj(2:np+1,2:nq+1,:) = intsj1;
            
            impmat = intsi.*intsj./(N_pqr.*(kc_sq - k_pqr_sq));
            
            
            impsum = sum(impmat(:));
            Z0(i,j) = -1i*rho*omega*impsum;
            if i~=j
                Z0(j,i) = Z0(i,j);
            end
        end
    end
%     for i = 1:patch_amt
%         for j = i:patch_amt
%             if i~=j
%                 Z0(j,i) = Z0(i,j);
%             end
%         end
%     end
    Z(:,:,l) = Z0;
    Y(:,:,l) = inv(Z0);
    
end

toc
%[row, column] = findlocation(42,n_x,n_y);
% matlabpool close

fac = 2e-5 ;
figure(1)
impedance = zeros(1,freqamt);
impedance(:) = Z(42,72,:);
plot(freq,20*log10(abs(impedance)/fac),'b')
hold on
impedance2 = zeros(1,freqamt);
impedance2(:) = Z(42,42,:);
plot(freq,20*log10(abs(impedance2)/fac),'m')
legend('transfer impedance','input impedance')
xlabel('Frequency (HZ)');
ylabel('Magnitude (dB)');
figure(2)
patchmobility = zeros(1,freqamt);
patchmobility(:) = Y(42,72,:);
plot(freq,20*log10(abs(patchmobility)),'b')
hold on
patchmobility2 = zeros(1,freqamt);
patchmobility2(:) = Y(42,42,:);
plot(freq,20*log10(abs(patchmobility2)),'m')
legend('transfer patch mobility','input patch mobility')
xlabel('Frequency (HZ)');
ylabel('Magnitude (dB)');

% yr = zeros(patch_amt,patch_amt,freqamt);
% for l=1:freqamt
%     Z0(:,:) = abs(Z(:,:,l));
%     yr(:,:,l) = inv(Z0);
% end
% figure(3)
% patchmobility3 = zeros(1,freqamt);
% patchmobility3(:) = yr(42,72,:);
% plot(freq,20*log10(abs(patchmobility3)),'b')
% hold on
% patchmobility4 = zeros(1,freqamt);
% patchmobility4(:) = yr(42,42,:);
% plot(freq,20*log10(abs(patchmobility4)),'m')
% legend('transfer patch mobility','input patch mobility')
% xlabel('Frequency (HZ)');
% ylabel('Magnitude (dB)');