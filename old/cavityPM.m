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
%matlabpool local 4
tic
l_patch( patch_amt ) = struct( 'x1',0,'y1',0,'x2',0,'y2',0,'xc',0,'yc',0,'z',0 ); % initialization
r_patch( patch_amt ) = struct( 'x1',0,'y1',0,'x2',0,'y2',0,'xc',0,'yc',0,'z',0 );
for i = 1:n_y
    for j = 1:n_x
        ptn = j+(i-1)*n_y;   % patch number
        lbx = (i-1)*del_x;
        lby = j*del_y;    % lbx lby: left-bottom;
        rtx = i*del_x;
        rty = (j-1)*del_y; %  rtx rty: right-top;
        cx = (lbx + rtx) / 2;
        cy = (lby + rty) / 2; %  cx cy: center;
        l_patch( ptn ) = struct( 'x1',lbx,'y1',lby,'x2',rtx,'y2',rty,'xc',cx,'yc',cy,'z',0 );
        r_patch( ptn ) = struct( 'x1',lbx,'y1',lby,'x2',rtx,'y2',rty,'xc',cx,'yc',cy,'z',Lz );
    end
end
freq = 150:152 ;
freqamt = length( freq );
Z0 = ones( patch_amt );
Z = ones( patch_amt, patch_amt, freqamt );
Y = ones( patch_amt, patch_amt, freqamt );
pi_lx = pi/Lx;
pi_ly = pi/Ly;
pi_lz = pi/Lz;
for l = 1:freqamt
    l
    omega = 2*pi*freq(l);
    kc_sq = (omega / cc)^2;
    for i = 1:patch_amt
        
        i
        for j = i:patch_amt
            %j
            impsum = 0;
            
            %[p,q,r] = meshgrid(1:30,1:30,1:20);
            for p = 0:30
%                 kpcomp = (p*pi/Lx)^2;
                kpcomp = (p*pi_lx)^2;
                if (p==0)
                    int_xcomp1 = del_x;
                    int_xcomp2 = del_x;
                    N_xcomp = Lx;
                else
%                     int_xcomp1 = Lx/(p*pi) * ( sin(p*pi*l_patch(i).x2/Lx)-sin(p*pi*l_patch(i).x1/Lx) );
%                     int_xcomp2 = Lx/(p*pi) * ( sin(p*pi*l_patch(j).x2/Lx)-sin(p*pi*l_patch(j).x1/Lx) );
                    int_xcomp1 = ( sin(p*pi_lx*l_patch(i).x2)-sin(p*pi_lx*l_patch(i).x1) ) / (p*pi_lx) ;
                    int_xcomp2 = ( sin(p*pi_lx*l_patch(j).x2)-sin(p*pi_lx*l_patch(j).x1) ) / (p*pi_lx) ;

                    N_xcomp = Lx/2;
                end
                for q = 0:30
%                     kqcomp = (q*pi/Ly)^2;
                    kqcomp = (q*pi_ly)^2;
                    if (q==0)
                        int_ycomp1 = del_y;
                        int_ycomp2 = del_y;
                        N_ycomp = Ly;
                    else
%                         int_ycomp1 = Ly/(q*pi) * ( sin(q*pi*l_patch(i).y2/Ly)-sin(q*pi*l_patch(i).y1/Ly) );
%                         int_ycomp2 = Ly/(q*pi) * ( sin(q*pi*l_patch(j).y2/Ly)-sin(q*pi*l_patch(j).y1/Ly) );
                        int_ycomp1 =  ( sin(q*pi_ly*l_patch(i).y2)-sin(q*pi_ly*l_patch(i).y1) ) / (q*pi_ly) ;
                        int_ycomp2 =  ( sin(q*pi_ly*l_patch(j).y2)-sin(q*pi_ly*l_patch(j).y1) ) / (q*pi_ly) ;
                        N_ycomp = Ly/2;
                    end
                    for r = 0:20
                        if (r==0)
                            N_zcomp = Lz;
                        else
                            N_zcomp = Lz/2;
                        end
                        N_pqr = N_xcomp * N_ycomp * N_zcomp ;
                        k_pqr_sq =  kpcomp +kqcomp +(r*pi_lz)^2 ;
                        
%                         intsi = cos(r*pi*l_patch(i).z/Lz)* int_xcomp1* int_ycomp1;
%                         intsj = cos(r*pi*l_patch(j).z/Lz)* int_xcomp2* int_ycomp2;
                        intsi = cos(r*pi_lz*l_patch(i).z)* int_xcomp1* int_ycomp1;
                        intsj = cos(r*pi_lz*l_patch(j).z)* int_xcomp2* int_ycomp2;
                        
                        impsum = impsum + intsi*intsj/(N_pqr*(kc_sq-k_pqr_sq));
                    end
                end
            end
            
            Z0(i,j) = -1i*rho*omega*impsum;
            if i~=j
                Z0(j,i) = Z0(i,j);
            end
        end
    end
    Z(:,:,l) = Z0;
    Y(:,:,l) = inv(Z0);
    
end

toc
%[row, column] = findlocation(42,n_x,n_y);
%matlabpool close

figure(1)
impedance = zeros(1,freqamt);
impedance(:) = Z(42,72,:);
plot(freq,20*log10(abs(impedance)),'b')
hold on
impedance2 = zeros(1,freqamt);
impedance2(:) = Z(42,42,:);
plot(freq,20*log10(abs(impedance2)),'m')
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