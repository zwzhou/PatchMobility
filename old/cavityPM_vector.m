close all;clc;clear;
% Acoustic patch mobility of a cavity -> FIG.4 & FIG.5
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


freq = 100:400 ;
freqamt = length( freq );

Z = ones( patch_amt, patch_amt, freqamt );
Y = ones( patch_amt, patch_amt, freqamt );
pi_lx = pi/Lx;
pi_ly = pi/Ly;
pi_lz = pi/Lz;
N_pqr_n0 = Lx*Ly*Lz/8;
N_pqr_10 = Lx*Ly*Lz/4;
N_pqr_20 = Lx*Ly*Lz/2;
N_pqr_a0 = Lx*Ly*Lz;
np = 30;    nq = 30;    nr = 12; % np nq & nr must >= 2 !
lefz = 0; rigz = Lz;
for l = 1:freqamt

    Z0 = ones( patch_amt );
    %disp(l)
    omega = 2*pi*freq(l);
    kc_sq = (omega / cc)^2;
    for i = 1:patch_amt
        [ix1,ix2,iy1,iy2] = findcoor(i,n_x,del_x,del_y);
        disp([freq(l),i])
        for j = i:patch_amt
            [jx1,jx2,jy1,jy2] = findcoor(j,n_x,del_x,del_y);
            %j
            %impsum = 0;
            
            % p,q,r > 0
            [p,q,r] = meshgrid(1:np,1:nq,1:nr);

            int_xcomp1 = ( sin(p*pi_lx*ix2)-sin(p*pi_lx*ix1) ) ./ (p*pi_lx) ;
            int_xcomp2 = ( sin(p*pi_lx*jx2)-sin(p*pi_lx*jx1) ) ./ (p*pi_lx) ;
%             int_xcomp1_c = ( sin(p*pi_lx*l_patch(i).x2)-sin(p*pi_lx*l_patch(i).x1) ) ./ (p*pi_lx) ;
%             int_xcomp2_c = ( sin(p*pi_lx*l_patch(j).x2)-sin(p*pi_lx*l_patch(j).x1) ) ./ (p*pi_lx) ;
%             debug1=int_xcomp1_c-int_xcomp1;
%             sumdebug1=sum(debug1(:));
%             debug2=int_xcomp2_c-int_xcomp2;
%             sumdebug2=sum(debug2(:));
            
            int_ycomp1 = ( sin(q*pi_ly*iy2)-sin(q*pi_ly*iy1) ) ./ (q*pi_ly) ;
            int_ycomp2 = ( sin(q*pi_ly*jy2)-sin(q*pi_ly*jy1) ) ./ (q*pi_ly) ;
            intsi = cos(r*pi_lz*lefz).* int_xcomp1.* int_ycomp1;
            intsj = cos(r*pi_lz*lefz).* int_xcomp2.* int_ycomp2;
            k_pqr_sq = (p*pi_lx).^2 + (q*pi_ly).^2 + (r*pi_lz).^2 ;
            impmat1 = intsi.*intsj./(N_pqr_n0*(kc_sq - k_pqr_sq));
            
            % r = 0 & p,q > 0
            [p,q] = meshgrid(1:np,1:nq);
            int_xcomp1 = ( sin(p*pi_lx*ix2)-sin(p*pi_lx*ix1) ) ./ (p*pi_lx) ;
            int_xcomp2 = ( sin(p*pi_lx*jx2)-sin(p*pi_lx*jx1) ) ./ (p*pi_lx) ;
            int_ycomp1 = ( sin(q*pi_ly*iy2)-sin(q*pi_ly*iy1) ) ./ (q*pi_ly) ;
            int_ycomp2 = ( sin(q*pi_ly*jy2)-sin(q*pi_ly*jy1) ) ./ (q*pi_ly) ;
            intsi = int_xcomp1.* int_ycomp1;
            intsj = int_xcomp2.* int_ycomp2;
            k_pqr_sq = (p*pi_lx).^2 + (q*pi_ly).^2 ;
            impmat2 = intsi.*intsj./(N_pqr_10*(kc_sq - k_pqr_sq));
            
            % q = 0 & p,r > 0 
            [p,r] = meshgrid(1:np,1:nr);
            int_xcomp1 = ( sin(p*pi_lx*ix2)-sin(p*pi_lx*ix1) ) ./ (p*pi_lx) ;
            int_xcomp2 = ( sin(p*pi_lx*jx2)-sin(p*pi_lx*jx1) ) ./ (p*pi_lx) ;
            int_ycomp1 = del_y ;
            int_ycomp2 = del_y ;
            intsi = cos(r*pi_lz*lefz).* int_xcomp1* int_ycomp1;
            intsj = cos(r*pi_lz*lefz).* int_xcomp2* int_ycomp2;
            k_pqr_sq = (p*pi_lx).^2 + (r*pi_lz).^2 ;
            impmat3 = intsi.*intsj./(N_pqr_10*(kc_sq - k_pqr_sq));
            
            % p = 0 & q,r > 0
            [q,r] = meshgrid(1:nq,1:nr);
            int_xcomp1 = del_x ;
            int_xcomp2 = del_x ;
            int_ycomp1 = ( sin(q*pi_ly*iy2)-sin(q*pi_ly*iy1) ) ./ (q*pi_ly) ;
            int_ycomp2 = ( sin(q*pi_ly*jy2)-sin(q*pi_ly*jy1) ) ./ (q*pi_ly) ;
            intsi = cos(r*pi_lz*lefz).* int_ycomp1* int_xcomp1;
            intsj = cos(r*pi_lz*lefz).* int_ycomp2* int_xcomp2;
            k_pqr_sq = (q*pi_ly).^2 + (r*pi_lz).^2 ;
            impmat4 = intsi.*intsj./(N_pqr_10*(kc_sq - k_pqr_sq));
            
            % r,q = 0 & p > 0
            p = 1:np;
            int_xcomp1 = ( sin(p*pi_lx*ix2)-sin(p*pi_lx*ix1) ) ./ (p*pi_lx) ;
            int_xcomp2 = ( sin(p*pi_lx*jx2)-sin(p*pi_lx*jx1) ) ./ (p*pi_lx) ;
            int_ycomp1 = del_y ;
            int_ycomp2 = del_y ;
            intsi = int_xcomp1* int_ycomp1;
            intsj = int_xcomp2* int_ycomp2;
            k_pqr_sq = (p*pi_lx).^2 ;
            impmat5 = intsi.*intsj./(N_pqr_20*(kc_sq - k_pqr_sq));
            
            % p,q = 0 & r > 0 
            r = 1:nr;
            int_xcomp1 = del_x ;
            int_xcomp2 = del_x ;
            int_ycomp1 = del_y ;
            int_ycomp2 = del_y ;
            intsi = cos(r*pi_lz*lefz)* int_xcomp1* int_ycomp1;
            intsj = cos(r*pi_lz*lefz)* int_xcomp2* int_ycomp2;
            k_pqr_sq = (r*pi_lz).^2 ;
            impmat6 = intsi.*intsj./(N_pqr_20*(kc_sq - k_pqr_sq));
            
            % p,r = 0 & q > 0
            q = 1:nq;
            int_xcomp1 = del_x ;
            int_xcomp2 = del_x ;
            int_ycomp1 = ( sin(q*pi_ly*iy2)-sin(q*pi_ly*iy1) ) ./ (q*pi_ly) ;
            int_ycomp2 = ( sin(q*pi_ly*jy2)-sin(q*pi_ly*jy1) ) ./ (q*pi_ly) ;
            intsi =  int_ycomp1* int_xcomp1;
            intsj =  int_ycomp2* int_xcomp2;
            k_pqr_sq = (q*pi_ly).^2 ;
            impmat7 = intsi.*intsj./(N_pqr_20*(kc_sq - k_pqr_sq));
            
            % p,q,r = 0
            impmat8 = (del_x*del_y)^2 / (N_pqr_a0*kc_sq);
            
            impsum = sum(impmat1(:))+sum(impmat2(:))+sum(impmat3(:))+sum(impmat4(:))+sum(impmat5(:))+sum(impmat6(:))+sum(impmat7(:))+impmat8;
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
figure(3)
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