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

np = 20;    nq = 20;    nr = 10; % np nq & nr must >= 2 !

freq = 100:0.5:120 ;
freqamt = length( freq );

Z = ones( patch_amt, patch_amt, freqamt );
Y = ones( patch_amt, patch_amt, freqamt );
lefz = 0; rigz = Lz;

% ------------------------
pi_lx = pi/Lx;pi_ly = pi/Ly;pi_lz = pi/Lz;
%---------------------------
% (0)  p,q,r > 0
[p,q,r] = meshgrid(1:np,1:nq,1:nr);
ppi_lx = p*pi_lx;
qpi_ly = q*pi_ly;
rpi_lz = r*pi_lz;
k_pqr_sq_0 = (ppi_lx).^2 + (qpi_ly).^2 + (rpi_lz).^2 ;
N_pqr_n0 = Lx*Ly*Lz/8;
ipo_0 = ( cos(rpi_lz*lefz)./ (ppi_lx.* qpi_ly) ).^2 ;

% (1)  r = 0 & p,q > 0
N_pqr_10 = Lx*Ly*Lz/4;
[p1,q1] = meshgrid(1:np,1:nq);
p1pi_lx = p1*pi_lx;
q1pi_ly = q1*pi_ly;
k_pqr_sq_1 = (p1pi_lx).^2 + (q1pi_ly).^2 ;
ipo_1 = ( 1./ (p1pi_lx.* q1pi_ly) ).^2 ;

% (2)  q = 0 & p,r > 0
[p2,r2] = meshgrid(1:np,1:nr);
p2pi_lx = p2*pi_lx;
r2pi_lz = r2*pi_lz;
k_pqr_sq_2 = (p2pi_lx).^2 + (r2pi_lz).^2 ;
ipo_2 = ( cos(r2pi_lz*lefz)./ (p2pi_lx) ).^2 ;

% (3)  p = 0 & q,r > 0
[q3,r3] = meshgrid(1:nq,1:nr);
q3pi_ly = q3*pi_ly;
r3pi_lz = r3*pi_lz;
k_pqr_sq_3 = (q3pi_ly).^2 + (r3pi_lz).^2 ;
ipo_3 = ( cos(r3pi_lz*lefz)./ (q3pi_ly) ).^2 ;

% (4.1)  r,q = 0 & p > 0
N_pqr_20 = Lx*Ly*Lz/2;
p4 = 1:np;
p4pi_lx = p4*pi_lx;
k_pqr_sq_41 = (p4pi_lx).^2 ;
ipo_41 = 1./ (p4*pi_lx).^2  ;

% (4.2)  p,q = 0 & r > 0
r4 = 1:nr;
r4pi_lz = r4*pi_lz;
k_pqr_sq_42 = (r4pi_lz).^2 ;
ipo_42 = ( cos(r4pi_lz*lefz) ).^2 ;

% (4.3)  p,r = 0 & q > 0
q4 = 1:nq;
q4pi_ly = q4*pi_ly;
k_pqr_sq_43 = (q4pi_ly).^2 ;
ipo_43 = ( 1./ (q4pi_ly) ).^2 ;

% (final)  p,q,r = 0
N_pqr_a0 = Lx*Ly*Lz;

for l = 1:freqamt

    Z0 = ones( patch_amt );
    %disp(l)
    omega = 2*pi*freq(l);
    kc_sq = (omega / cc)^2;
    
    % (0)  p,q,r > 0
    immid0 = ipo_0./ (N_pqr_n0*(kc_sq - k_pqr_sq_0));
    % (1)  r = 0 & p,q > 0
    immid1 = ipo_1./ (N_pqr_10*(kc_sq - k_pqr_sq_1));
    % (2)  q = 0 & p,r > 0
    immid2 = ipo_2./ (N_pqr_10*(kc_sq - k_pqr_sq_2));
    % (3)  p = 0 & q,r > 0
    immid3 = ipo_3./ (N_pqr_10*(kc_sq - k_pqr_sq_3));
    % (4.1)  r,q = 0 & p > 0
    immid41 = ipo_41./ (N_pqr_20*(kc_sq - k_pqr_sq_41));
    % (4.2)  p,q = 0 & r > 0
    immid42 = ipo_42./ (N_pqr_20*(kc_sq - k_pqr_sq_42));
    % (4.3)  p,r = 0 & q > 0
    immid43 = ipo_43./ (N_pqr_20*(kc_sq - k_pqr_sq_43));
    % (final)  p,q,r = 0
    impmat8 = (del_x*del_y)^2 / (N_pqr_a0*kc_sq);
            
    for i = 1:patch_amt
        disp([freq(l),i])
        [ix1,ix2,iy1,iy2] = findcoor(i,n_x,del_x,del_y);
        
        % (0)  p,q,r > 0
        int_xcomp1_0 = ( sin(ppi_lx*ix2)-sin(ppi_lx*ix1) ) ;
        int_ycomp1_0 = ( sin(qpi_ly*iy2)-sin(qpi_ly*iy1) ) ;
        intsi_0 = int_xcomp1_0.* int_ycomp1_0 ;
        
        % (1)  r = 0 & p,q > 0
        int_xcomp1_1 = ( sin(p1pi_lx*ix2)-sin(p1pi_lx*ix1) )  ;
        int_ycomp1_1 = ( sin(q1pi_ly*iy2)-sin(q1pi_ly*iy1) )  ;
        intsi_1 =  int_xcomp1_1.* int_ycomp1_1;
        
        % (2)  q = 0 & p,r > 0
        int_xcomp1_2 = ( sin(p2pi_lx*ix2)-sin(p2pi_lx*ix1) )  ;
        int_ycomp1_2 = del_y ;
        int_ycomp2_2 = del_y ;
        intsi_2 =  int_xcomp1_2* int_ycomp1_2;
        
        % (3)  p = 0 & q,r > 0
        int_xcomp1_3 = del_x ;
        int_xcomp2_3 = del_x ;
        int_ycomp1_3 = ( sin(q3pi_ly*iy2)-sin(q3pi_ly*iy1) )  ;
        intsi_3 =  int_ycomp1_3* int_xcomp1_3;
        
        % (4.1)  r,q = 0 & p > 0
        int_xcomp1_41 = ( sin(p4pi_lx*ix2)-sin(p4pi_lx*ix1) )  ;
        int_ycomp1_41 = del_y ;
        int_ycomp2_41 = del_y ;
        intsi_41 = int_xcomp1_41* int_ycomp1_41;
        
        % (4.2)  p,q = 0 & r > 0
        int_xcomp1_42 = del_x ;
        int_xcomp2_42 = del_x ;
        int_ycomp1_42 = del_y ;
        int_ycomp2_42 = del_y ;
        intsi_42 =  int_xcomp1_42* int_ycomp1_42;
        intsj_42 =  int_xcomp2_42* int_ycomp2_42;
        impmat6 = immid42.* intsi_42.*intsj_42;
        sumimpmat6 = sum(impmat6(:));
        
        % (4.3)  p,r = 0 & q > 0
        int_xcomp1_43 = del_x ;
        int_xcomp2_43 = del_x ;
        int_ycomp1_43 = ( sin(q4pi_ly*iy2)-sin(q4pi_ly*iy1) )  ;
        intsi_43 =  int_ycomp1_43* int_xcomp1_43;
        
        % (final)  p,q,r = 0
        
        for j = i:patch_amt
            [jx1,jx2,jy1,jy2] = findcoor(j,n_x,del_x,del_y);

            % (0)  p,q,r > 0
            int_xcomp2_0 = ( sin(ppi_lx*jx2)-sin(ppi_lx*jx1) )  ;
            int_ycomp2_0 = ( sin(qpi_ly*jy2)-sin(qpi_ly*jy1) )  ;
            intsj_0 =  int_xcomp2_0.* int_ycomp2_0;
            impmat0 = immid0.* intsi_0.*intsj_0 ;
            
            % (1)  r = 0 & p,q > 0
            int_xcomp2_1 = ( sin(p1pi_lx*jx2)-sin(p1pi_lx*jx1) )  ;
            int_ycomp2_1 = ( sin(q1pi_ly*jy2)-sin(q1pi_ly*jy1) )  ;
            intsj_1 =  int_xcomp2_1.* int_ycomp2_1;
            impmat1 = immid1.* intsi_1.*intsj_1 ;
            
            % (2)  q = 0 & p,r > 0 
            int_xcomp2_2 = ( sin(p2pi_lx*jx2)-sin(p2pi_lx*jx1) )  ;
            intsj_2 =  int_xcomp2_2* int_ycomp2_2;
            impmat2 = immid2.* intsi_2.*intsj_2 ;
            
            % (3)  p = 0 & q,r > 0
            int_ycomp2_3 = ( sin(q3pi_ly*jy2)-sin(q3pi_ly*jy1) )  ;
            intsj_3 =  int_ycomp2_3* int_xcomp2_3;
            impmat3 = immid3.* intsi_3.*intsj_3;
            
            % (4.1)  r,q = 0 & p > 0
            int_xcomp2_41 = ( sin(p4pi_lx*jx2)-sin(p4pi_lx*jx1) )  ;
            intsj_41 = int_xcomp2_41* int_ycomp2_41;
            impmat5 = immid41.* intsi_41.*intsj_41 ;
            
            % (4.2)  p,q = 0 & r > 0 
            
            % (4.3)  p,r = 0 & q > 0
            int_ycomp2_43 = ( sin(q4pi_ly*jy2)-sin(q4pi_ly*jy1) )  ;
            intsj_43 =  int_ycomp2_43* int_xcomp2_43;
            impmat7 = immid43.* intsi_43.*intsj_43;
            
            % (final)  p,q,r = 0
            
            impsum = sum(impmat0(:))+sum(impmat1(:))+sum(impmat2(:))+sum(impmat3(:))+sum(impmat5(:))+sumimpmat6+sum(impmat7(:))+impmat8;
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