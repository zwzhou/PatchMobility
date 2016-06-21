close all;clc;clear;
% Acoustic patch mobility of a cavity -> FIG.4 & FIG.5
Lx = 1.5; Ly = 0.96; Lz = 0.01;
n_x = 19 ; n_y = 13 ;% patch number of x and y
patch_amt = n_x * n_y ;%  amount of patch
del_x = Lx / n_x ;  del_y = Ly / n_y ; % patch dimension
rho = 1.293 ; % kg/m^3
c_air=343.6; %m/s, speed of sound in air 20
eta = 0.01;
cc = c_air*sqrt(1+1i*eta);

%Initialize Matlab Parallel Computing Enviornment by Xaero | Macro2.cn
%CoreNum=8; %设定机器CPU核心数量
% if matlabpool('size')<=0 %判断并行计算环境是否已然启动
% matlabpool('open','local',CoreNum); %若尚未启动，则启动并行环境
% else
% disp('Already initialized'); %说明并行环境已经启动。
% end
%
%matlabpool local 4

%matlabpool(7)

l_patch( patch_amt ) = struct( 'x1',0,'y1',0,'x2',0,'y2',0,'xc',0,'yc',0,'z',0 ); % initialization
r_patch( patch_amt ) = struct( 'x1',0,'y1',0,'x2',0,'y2',0,'xc',0,'yc',0,'z',0 );
lpatch_x1=zeros(1,patch_amt);lpatch_x2=zeros(1,patch_amt);
lpatch_y1=zeros(1,patch_amt);lpatch_y2=zeros(1,patch_amt);
lpatch_z=zeros(1,patch_amt);
for i = 1:n_y
    for j = 1:n_x
        ptn = j+(i-1)*n_y;   % patch number
        lbx = (j-1)*del_x;
        lby = (i-1)*del_y;    % lbx lby: left-bottom;
        rtx = j*del_x;
        rty = i*del_y; %  rtx rty: right-top;
        cx = (lbx + rtx) / 2;
        cy = (lby + rty) / 2; %  cx cy: center;
        l_patch( ptn ) = struct( 'x1',lbx,'y1',lby,'x2',rtx,'y2',rty,'xc',cx,'yc',cy,'z',0 );
        r_patch( ptn ) = struct( 'x1',lbx,'y1',lby,'x2',rtx,'y2',rty,'xc',cx,'yc',cy,'z',Lz );
        lpatch_x1( ptn ) = lbx;
        lpatch_x2( ptn ) = rtx;
        lpatch_y1( ptn ) = lby;
        lpatch_y2( ptn ) = rty;
        lpatch_z( ptn ) = 0 ;
    end
end
freq = 100:600 ;
freqamt = length( freq );
%N_pqr = Lx*Ly*Lz /8;
impedance = zeros(1,freqamt);
i = 42;
lpatch_x2_42 = lpatch_x2(i);
lpatch_x1_42 = lpatch_x1(i);
lpatch_y2_42 = lpatch_y2(i);
lpatch_y1_42 = lpatch_y1(i);
lpatch_z_42 = lpatch_z(i);
% temp = 20:10:120;
% tempamt = length(temp);
% for ii = 1:tempamt
% nn = temp(ii);
nn = 50;
tic
for l = 1:freqamt
    l
    omega = 2*pi*freq(l);
    kc = omega / cc;
    impsum = 0;
            for p = 0:nn
                if (p==0)
                    xcomp = del_x;
                    N_xcomp = Lx;
                else
                    xcomp = Lx/(p*pi) * ( sin(p*pi*lpatch_x2_42/Lx)-sin(p*pi*lpatch_x1_42/Lx) );
                    N_xcomp = Lx/2;
                end
                for q = 0:nn
                    if (q==0)
                        ycomp = del_y;
                        N_ycomp = Ly;
                    else
                        ycomp = Ly/(q*pi) * ( sin(q*pi*lpatch_y2_42/Ly)-sin(q*pi*lpatch_y1_42/Ly) );
                        N_ycomp = Ly/2;
                    end
                    for r = 0:nn
                        if (r==0)
                            N_zcomp = Lz;
                        else
                            N_zcomp = Lz/2;
                        end
                        N_pqr = N_xcomp * N_ycomp * N_zcomp ;
                        k_pqr = sqrt( (p*pi/Lx)^2+(q*pi/Ly)^2+(r*pi/Lz)^2 );
                        intsi = cos(r*pi*lpatch_z_42/Lz)* xcomp* ycomp;
                        impsum = impsum + intsi*intsi/(N_pqr*(kc^2-k_pqr^2));
                    end
                end
            end
            impedance(l) = -1i*rho*omega*impsum;

    
end
% 
% fileid = fopen('par.txt','a');
% fprintf(fileid,'%8.1f ',nn);
% fprintf(fileid,'%18.12f',log10(abs(impedance)));
% fprintf(fileid,' \n');
% fclose(fileid);
% end

toc
%matlabpool close
plot(freq,20*log10(abs(impedance/2e-5)))