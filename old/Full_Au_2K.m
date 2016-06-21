%Transmission Loss
%Double-leaf window incorperated with m TARs
%Determine the resonator location
clear
close all

nFull_Couple_Flag = 0 % Fully structural and acoustic coupling model
%nFull_Couple_Flag=0 % Not fully structural and acoustic coupling model

%Aproximates integral using Gauss-Legendre quadrature method
cti = cputime;

%Physical parameters
%Aluminum Plates with Lsx=0.83m, Lsy=1.2782m, thick=0.0141m
Lsx=1; %m; length
Lsy=Lsx*1.54;%0.83; %m; width
thick_i=Lsx*1.7/100;%0.0136;%0.003; %m; thickness
thick_r=thick_i; %m; thickness
ros_i=2700;%2373; % kg/m3  Aluminum 6061-T6 2700 kg/m^3 
ros_r=ros_i; %2700;%kg/m3
E_i=7.310e10;%60e9; %Pa, Young's modulus   Aluminum 6061-T6 7.310E+10 Pa 
E_r=E_i; %7.310e10;%Pa, Young's modulus
mu_si=0.33;%0.22; %Poisson ratio   Aluminum 6061-T6 0.3300 
mu_sr=mu_si; %0.33;%Poisson ratio   Aluminum 6061-T6 0.3300 
theta_c=0.0028; % acoustic damping ratio
theta_si=0.005;%0.0015; %structural damping ratio
theta_sr=theta_si;

mI=ros_i*thick_i; %kg/m2, Surface density
mR=ros_r*thick_r;  %kg/m2, Surface density
S_str=Lsx*Lsy; %m2, Area of the plates

%Incident angle of the incident plane wave
Pin=1;
inc_theta=pi/3;
inc_fa=pi/6;

%Concentrated mass at (xM,yM)
%M=mI*Lsx*Lsy*0.01;
%M=mI*Lsx*Lsy*0.1; %
M=mI*Lsx*Lsy*1;
%xM=Lsx/8;%Lsx/2;
%yM=Lsy/7;%Lsy/2;
xM=Lsx/2;
yM=Lsy/2;

%Concentrated center-point Force Q0 at (xF,yF)
F=1;
xF=Lsx/10;%Lsx/8;%Lsx/2;
yF=Lsy/14;%Lsy/7;%Lsy/2;

%F1=F*Lsx/D;

%Cavity
Lcx=Lsx;%m; length
Lcy=Lsy; %m; width
Lcz=0.1; %m; height
VC=Lcx*Lcy*Lcz; %m3, Volume of the enclosure

ro_air=1.21; %kg/m3
c_air=343.6; %m/s, speed of sound in air 20
ro_water=998.2; %kg/m3     998.2  at 20¡ãC
c_water=1482; %m/s, speed of sound in water, 1482 m/s at 20¡ãC

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Chang medium   1 of 5!!
ro1=ro_air; %         Incident plate side
c1=c_air; %m/s, speed of sound in air 20  
ro2=ro_air; %ro1;%kg/m3 %                          Cavity
c2=c_air; %c1;%m/s, speed of sound in air 20
ro3=ro_water;%ro_air; %kg/m3 water    998.2  at 20¡ãC    Radiation plate side 
c3=c_water;%c_air; %m/s, speed of sound in water, 1482 m/s at 20¡ãC

%Cavity
%Natural frequencies of cavity
w0C0=[];
pp0=0;
nL=20;%25;%10;%8;%6;
nM=20;%25;%10;%8;%6;
nN=3;%4;%6;%1;
J1=nL*nM*nN;
for l=1:1:nL
    kx=(l-1)*pi/Lcx;
    for m=1:1:nM
        ky=(m-1)*pi/Lcy;
        for n=1:1:nN
            kz=(n-1)*pi/Lcz;
            kr=sqrt(kx^2+ky^2+kz^2);
            Radia_Freq=c2*kr;
            %index vector
            pp0=1+pp0;
            %frequency vector
            w0C0=[w0C0 Radia_Freq];
            eigen_C0(pp0)=struct('index', pp0,'frequency',Radia_Freq,'l',l-1, 'm',m-1, 'n',n-1);
        end
    end
end

%Sorting frequency
[YC0, IC0]=sort(w0C0);

for j=1:J1
    C_l(j)=eigen_C0(IC0(j)).l;
    C_m(j)=eigen_C0(IC0(j)).m;
    C_n(j)=eigen_C0(IC0(j)).n;
    C_freq(j)=eigen_C0(IC0(j)).frequency;
end
kl=C_l*(pi/Lcx);
km=C_m*(pi/Lcy);
kn=C_n*(pi/Lcz);

%Modal mass of the cavity. It is 1 if the mode shape is normalized.
%for j=1:J1
%    if C_l(j)==0 
%        Al=1;
%    else
%        Al=0.5;
%    end
%    if C_m(j)==0
%        Am=1;
%    else
%        Am=0.5;
%    end
%    if  C_n(j)==0
%        An=1;
%    else
%        An=0.5;
%    end
%     Aj(j)=Al*Am*An*VC;
%end

%Factors for normalizing mode shapes
%Cavity
for j=1:J1
    if C_l(j)==0 
        Al=1;   %Lcx
    else
        Al=2;   %Lcx/2
    end
    if C_m(j)==0
        Am=1;  %Lcy
    else
        Am=2;  %Lcy/2
    end
    if  C_n(j)==0
        An=1;  %Lcz
    else
        An=2;  %Lcz/2
    end
     Aj(j)=sqrt(Al*Am*An/VC);
end

%Incident plate
%Radiating plate
%Natural frequencies of simply suported plates
w0S0i=[];
w0S0r=[];
pp1=0;
nP=10;%8;
nQ=10;%8;
J2=nP*nQ;
J3=J2;
for p=1:nP
    kxp=p*pi/Lsx;
    for q=1:nQ
        
        kyq=q*pi/Lsy;
        Struc_Freqi=(kxp^2+kyq^2)*sqrt(E_i*thick_i^2/(12*ros_i*(1-mu_si^2)));
        Struc_Freqr=(kxp^2+kyq^2)*sqrt(E_r*thick_r^2/(12*ros_r*(1-mu_sr^2)));

        %index vector
        pp1=1+pp1;
        %frequency vector
        w0S0i=[w0S0i Struc_Freqi];
        eigen_S0i(pp1)=struct('index', pp1,'frequency',Struc_Freqi, 'p',p, 'q',q);
        
        w0S0r=[w0S0r Struc_Freqr];
        eigen_S0r(pp1)=struct('index', pp1,'frequency',Struc_Freqr, 'p',p, 'q',q);
    end
end

%Sorting frequency
[YS0i, IS0i]=sort(w0S0i);
[YS0r, IS0r]=sort(w0S0r);

for j=1:J2
    Si_p(j)=eigen_S0i(IS0i(j)).p;
    Si_q(j)=eigen_S0i(IS0i(j)).q;
    Si_freq(j)=eigen_S0i(IS0i(j)).frequency;
    
    Sr_p(j)=eigen_S0r(IS0r(j)).p;
    Sr_q(j)=eigen_S0r(IS0r(j)).q;
    Sr_freq(j)=eigen_S0r(IS0r(j)).frequency;
end

kpi=Si_p*(pi/Lsx);
kqi=Si_q*(pi/Lsy);
kpr=Sr_p*(pi/Lsx);
kqr=Sr_q*(pi/Lsy);

%Modal mass of plates: For normalized mode shape, it is mI or mR.
%for n=1:J2  
%    MnI(n)=ros_i*thick_i*Lsx*Lsy/4;
%    MnR(n)=ros_r*thick_r*Lsx*Lsy/4;
%end
%Factors for normalizing mode shapes
%Plates
AI(1:J2)=2/sqrt(S_str); %Simmply supported
AR(1:J2)=2/sqrt(S_str); %Simply supported

%Normalized mode shapes
%The mode shapes at (xM, yM) and (xF,yF) of the incident plate
FaS0_Mi=AI.*sin(kpi*xM).*sin(kqi*yM);
FaS0_Fi=AI.*sin(kpi*xF).*sin(kqi*yF);

%Interaction between plates and cavity
%Coupling factors for not normalized modes
%for j=1:J1
%    for n=1:J2
%        if (Si_p(n) == C_l(j))
%            Clp=0;
%        else
%            Clp= Si_p(n)*((-1)^(C_l(j)+Si_p(n))-1)/(C_l(j)^2-Si_p(n)^2);
%        end
%        if (Si_q(n) == C_m(j))
%            Cmq=0;
%        else
%            Cmq=Si_q(n)*((-1)^(C_m(j)+Si_q(n))-1)/(C_m(j)^2-Si_q(n)^2);
%        end
%        CR(j,n)=S_str*Clp*Cmq/pi^2;
%        CI(j,n)= CR(j,n)*(-1)^C_n(j);
%    end
%end

%Coupling factors for normalized mode shapes of cavity and plates
for j=1:J1
    for n=1:J2
        if (Si_p(n) == C_l(j))
            Clp=0;
        else
            Clp= Si_p(n)*((-1)^(C_l(j)+Si_p(n))-1)/(C_l(j)^2-Si_p(n)^2);
        end
        if (Si_q(n) == C_m(j))
            Cmq=0;
        else
            Cmq=Si_q(n)*((-1)^(C_m(j)+Si_q(n))-1)/(C_m(j)^2-Si_q(n)^2);
        end
        CR(j,n)=S_str*Clp*Cmq/pi^2*(Aj(j)*AR(n));
        CI(j,n)=S_str*Clp*Cmq/pi^2*(-1)^C_n(j)*(Aj(j)*AI(n));
    end
end

%%%%
%Call modal transffer resistance matrix Eij
if nFull_Couple_Flag ==1
    load D:\Research_PolyUME\Me_Papers\Newdraft\Simulation\vibroacoustics\E_Au_airXY_2K.mat;
    load D:\Research_PolyUME\Me_Papers\Newdraft\Simulation\vibroacoustics\E_Au_waterXY_2K.mat;
else
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
    %%No fully coupling
    E_Au_airXY_2K=zeros(2000,J2,J2);
    E_Au_waterXY_2K=zeros(2000,J2,J2);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Geometric grid for the calculation of transmission power
    Ix=40;  %number in x-direction
    Iy=40;  %number in y-direction
    nPointNo=Ix*Iy; %total point number

    dx=Lsx/(2*Ix);
    dy=Lsy/(2*Iy);
    S_ele=4*dx*dy;  %area of each element

    ii=[1:Ix];
    xx=(2*ii-1)*dx;
    jj=[1:Iy];
    yy=(2*jj-1)*dy;
    [x,y]=meshgrid(xx,yy);
    %Distance between point ii1 and ii2
    for ii1=1:nPointNo
        for ii2=1:nPointNo
            rr(ii1,ii2)=sqrt((x(ii1)-x(ii2))^2+(y(ii1)-y(ii2))^2);  
        end
    end
    S_ele=4*dx*dy;  %area of each element
end

%Freq. Band of Interest
f=[1:2000];
nFrequency=length(f);

for idf=1:nFrequency
    w=2*pi*f(idf);
    k_in=w/c1;  %incident plate side
    k_c=w/c2;  %cavity
    k_ra=w/c3;  %Radiation plate side
        
    %Incident modal force
    kx=k_in*sin(inc_theta)*cos(inc_fa);
    ky=k_in*sin(inc_theta)*sin(inc_fa);
    for n=1:J2
        if (kx)^2==(kpi(n))^2
            Ixx=(-1i/2)*sign(sin(inc_theta)*cos(inc_fa));
        else
            Ixx=Si_p(n)*pi*(1-(-1)^Si_p(n)*(cos(kx*Lsx)-1i*sin(kx*Lsx)))/...
                ((Si_p(n)*pi)^2-(kx*Lsx)^2);
        end
        if (ky)^2==(kqi(n))^2
             Iyy=(-1i/2)*sign(sin(inc_theta)*sin(inc_fa));
        else
             Iyy=Si_q(n)*pi*(1-(-1)^Si_q(n)*(cos(ky*Lsy)-1i*sin(ky*Lsy)))/...
                ((Si_q(n)*pi)^2-(ky*Lsy)^2);
        end
        PIN(n)=2*Pin*S_str*AI(n)*Ixx*Iyy;  %Normalized mode shapes, plane wave
    end
    
    %%%%%%%%%%%%%
    %Concentrated center-point Force Q0 at (xe,ye)
    %PIN=F*FaS0_Fi;
    %PIN0=PIN0+PIN1;
 
    %Construct A matrix with resonator

    %The two plates have the same geometric and physical properties
    ZI=(theta_si*Si_freq.*(Si_freq/w)+1i*w*((Si_freq/w).^2-1))*mI;%(2*theta_si*Si_freq+i*w*((Si_freq/w).^2-1))*mI;  %Normalized mode shapes
    ZM=-1i*w*M*(FaS0_Mi)'*FaS0_Mi;   %Induced by the concentrated mass
    ZR=(theta_sr*Sr_freq.*(Sr_freq/w)+1i*w*((Sr_freq/w).^2-1))*mR;%(2*theta_sr*Sr_freq+i*w*((Sr_freq/w).^2-1))*mR;  %Normalized mode shapes
    YC=(2*theta_c*C_freq/w+1i*((C_freq/w).^2-1))*w/(ro2*c2^2);  %Normalized mode shapes

    %Cavity forcing term
    for j=1:J1
        for h=1:J1
            if h==j
                %Term of cavity
                a1=YC(h);
            else
                a1=0;
            end
            
            %No fully coupling
            %Term of incident plate
     %       a2=0;
     %       for r=1:J2
     %           a2=a2+CI(j,r)*CI(h,r)/KI(r);
     %       end
          
            %Term of radiation plate
     %       a3=0;
     %       for s=1:J3
     %           a3=a3+CR(j,s)*CR(h,s)/KR(s);
     %       end
          
            %Term of resonators
    %        a4=0;
    %        if Resonator ~= 0
    %            for m=1:MM
    %                a4=a4+FaC0_R(m,j)*FaC0_R(m,h)/Zr(m);
    %            end
    %        end
    %        AA(j,h)=a1+a2+a3;
            A(j,h)=a1;
        end
       
         %Terms of incident plates
         %CI(j,n);    %J1xJ2
         %CR(j,n);    %J1xJ2      
    end
   
    %Incident plate equation, considering fully coupling
    for m=1:J2
        %Term of incident plate
        for q=1:J2
           if  q==m
               a1=ZI(q);
           else
               a1=0;
           end
           %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
           %Chang medium   2 of 5!!
           BI(q,m)=a1+ZM(q,m)+E_Au_airXY_2K(idf,q,m);  %EI J2xJ2
        end
        %TCI(m,j)               %J2xJ1
    end
    BI=conj(BI'); 
    TCI=conj(CI');   %Transport

    %Radiation plate equation, considering fully coupling
    for m=1:J2
        %Term of radiation plate
        for q=1:J2
           if  q==m
               a1=ZR(q);
           else
               a1=0;
           end
           %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
           %Chang medium   3 of 5!!
           BR(q,m)=a1+E_Au_waterXY_2K(idf,q,m);   %ER J2xJ2
           %BR(q,m)=a1+E_Au_airXY_2K(idf,q,m);   %ER J2xJ2
        end
        %TCR(m,j)                %J2xJ1
    end
    BR=conj(BR');   %Transport
    TCR=conj(CR');   %Transport
   
    INVBI=inv(BI);
    INVBR=inv(BR);
    AA=A+CI*INVBI*TCI+CR*INVBR*TCR;
    BB=CI*INVBI*conj(PIN');
   
    %   PC0=AA\BB';
    %PC0=pinv(AA)*BB';
    %PC0=AA\conj(BB');  %J2 x 1
    PC0=AA\BB;  %J2 x 1

    %Modal velocity response of incident panel
    VI=INVBI*(conj(PIN')-TCI*PC0);
    VI=conj(VI');
    %Modal velocity response of radiating panel
    VR=INVBR*TCR*PC0; %J2 x 1
    VR=conj(VR');
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%No fully coupling
    
    %a1=conj(CI')*PC0;
    %a1=conj(a1');
    %VI=(PIN-a1)./KI;  %1xJ2

    %a1=conj(CR')*PC0;
    %a1=conj(a1');
    %VR=a1./KR; %1xJ2
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %pressure response at microphone 
    %p_C0(idf)=FaC0_M*PC0;

    %Transmission loss
    if nFull_Couple_Flag==1
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %Chang medium   4 of 5!!
        EE=squeeze(E_Au_waterXY_2K(idf,:,:));
        %EE=squeeze(E_Au_airXY_2K(idf,:,:));
        aab(idf)=real(VR*EE*VR'); %There is VR', but not conj(VR')
        TL(idf)=10*log10(S_str*cos(inc_theta)*(Pin)^2/(ro1*c1)/aab(idf));  
    else
        for ii=1:nPointNo 
             %Vibrating velocity at (xx, yy)
             %and the value of plate mode shape at (xx, yy)
             %Normalized mode shapes
             FaS0_Wi =AI.*sin(kpi*x(ii)).*sin(kqi*y(ii));
             FaS0_Wr =AR.*sin(kpr*x(ii)).*sin(kqr*y(ii));
             vi(ii)=VI*FaS0_Wi';
             vr(ii)=VR*FaS0_Wr';
        end

        %Reflection and Transmission (Transffer resistance matrix)
        for ii1=1:1:nPointNo
             for ii2=1:1:nPointNo
                 if (ii1==ii2)
                     RR(ii1,ii2)=ro3*w^2*S_ele^2/4/pi/c3;
                 else
                     RR(ii1,ii2)=(ro3*w^2*S_ele^2/4/pi/c3)*sin(k_ra*rr(ii1,ii2))/(k_ra*rr(ii1,ii2));
                 end
             end
        end
        %Transmission power
        PWR(idf)=vr*RR*vr';   %p*v', no tranport: conj(vr'). conjugate and transport:v'
 
        %Reflection power
        %PWI(idf)=vi*RR*vi';
    end
       
    %Averaged quadratic velocity at w     Normalized mode shapes
    Av=(VI*VI')/(2*S_str);  %W x conj(W)
    Av_vi(idf)=10*log10(Av/2.5e-15);
    Av=(VR*VR')/(2*S_str);   
    Av_vr(idf)=10*log10(Av/2.5e-15);
       
    %Averaged sound pressure level Lp at w in the enclsoure    Noemalized
    %mode shapes
    Av_p=(PC0'*PC0)/(2*VC);   %a=[1+2i 3+4i] =>a'=[1-2i;3-4i]=>conj(a')=[1+2i;3+4i]
    %conj(a)=[1-2i 3-4i]
    Av_Lp(idf)=10*log10(Av_p/0.00002^2);
end

if nFull_Couple_Flag~=1
    %Pure incident power
    PWI=(Pin)^2/(2*ro1*c1)*S_str*cos(inc_theta);
    
    %Transmission loss
    TL=10*log10(PWI./PWR); 
%    TL=10*log10(PWI./abs(PWR));
end

figure(1)
plot(f,TL, 'LineStyle', '-','LineWidth', 2, 'Color', [0 0 1])
xlabel('Frequency (Hz)')
ylabel('TL (dB)')
    
figure(2)
plot(f,Av_vi, 'LineStyle', '-','LineWidth', 2, 'Color', [0 0 1])
hold on    
plot(f,Av_vr, 'LineStyle', '--','LineWidth', 2, 'Color', [1 0 0])
xlabel('Frequency (Hz)')
ylabel('AQV (dB)')
legend('Incident plate','Radiation plate')
hold off
    
figure(3)
plot(f,Av_Lp, 'LineStyle', '-','LineWidth', 2, 'Color', [0 0 1])
xlabel('Frequency (Hz)')
ylabel('AP in enclosure (dB)')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Chang medium   5 of 5!!
%ss=['save D:\Research_PolyUME\Me_Papers\Newdraft\Simulation\vibroacoustics\data\Full_air_Au_air_Au_water_mass_2K_m20p_Distr f Av_Lp Av_vi Av_vr TL'];
%ss=['save C:\Research_PolyUME\Me_Papers\Newdraft\Simulation\vibroacoustics\TL_full_air_Au_air_Au_water_nomass f Av_Lp Av_vi Av_vr TL'];
%eval(ss)

ec_ti = cputime-cti


