clear;clc;
% VI. SEMI-INFINITE MEDIUM  -->  A.Radiated pressure  -->  FIG.7
% disp('Calculating Radiation Impedance Matrix Z')
rho_a = 1.293 ; % kg/m^3  density of air
c = 340  ;  % m/s
dij1 = 0.08;
dij2 = 0.84;

lx = 1.5; ly = 0.96; % panel size
nx = 19;  ny = 13; % x y patch count
patchamt = nx * ny ;

freq = 1:1:1000;
freqnum = length( freq );

%----------------------------------------------------------------------------------
tic
[Zrad,~] = f_radZ( lx,ly, nx,ny, rho_a,c,freq );
toc
disp('Inversing...')
YP = ones( patchamt, patchamt, freqnum );% initialization
for ii = 1:freqnum
    YP(:,:,ii) = inv(squeeze(Zrad(:,:,ii)));
end
tic
Zrad1 = ones( patchamt, patchamt, freqnum );% initialization
for ii = 1:freqnum
    Zrad1(:,:,ii) = f_radZ(lx,ly, nx,ny, rho_a,c,freq(ii) );
end
toc
fprintf('Complete!\nprinting...\n')
%----------------------------------------------------------------------------------


%----------------------------------------------------------------------------------
% delx = lx / nx;  dely = ly / ny; % patch size
% S_j = delx * dely; % patch aera
% ai = abs( sqrt(S_j/pi) ); % patch radius
% 
% % initialization
% patchnumber = 1:patchamt;
% patch = f_findcoor(patchnumber,nx,delx,dely);
% xdis = bsxfun(@minus,patch.x2.',patch.x2);% patch x distance
% ydis = bsxfun(@minus,patch.y2.',patch.y2);% patch y distance
% dis = sqrt( xdis.^2 + ydis.^2 );% patch distance
% 
% tic
% 
% Z = zeros( patchamt ); % initialization
% 
% omega = 2*pi*freq;
% k = omega / c;
% % parpool
% fprintf('current: %6.1f / ( start: %4.1f   step: %3.1f   end: %5.1f )\n',...
%     freq(1), freq(1), freq(2)-freq(1), freq(freqnum) )
% for ii = 1:freqnum
%     if (~mod(ii,100))
%         fprintf('current: %6.1f / ( start: %4.1f   step: %3.1f   end: %5.1f )\n',...
%             freq(ii), freq(1), freq(2)-freq(1), freq(freqnum) )
%     end
%     
% %     omega = 2*pi*freq(l);
% %     k = omega / c;
% %     Z(dis==0) = rho*c*(1-exp(-1i*k*ai));
% %     Z(dis~=0) = omega*rho*1i/2/pi*S_j* exp(-1i*k*dis(dis~=0)) ./ dis(dis~=0);
%     
%     Z(dis==0) = rho_a*c*(1-exp(-1i*k(ii)*ai));
%     Z(dis~=0) = omega(ii)*rho_a*1i/2/pi*S_j* exp(-1i*k(ii)*dis(dis~=0)) ./ dis(dis~=0);
% 
% %     Z(patchamt,patchamt) = rho*c*( 1-exp(-1i*k*ai) );
% %     for i = 1 : patchamt-1
% %         Z(i,i) = Z(patchamt,patchamt);
% %         for j = i+1 : patchamt
% %             distance = dis(i,j) ;
% %             Z(i,j) = omega*rho*1i*exp(-1i*k*distance)*S_j / (2*pi*distance);
% %             Z(j,i) = Z(i,j);
% %         end
% %     end
% %     
%     YP(:,:,ii) = inv(Z);
% end 
% toc
%-----------------------------------------------------------------------------------

% Yii = zeros( 1,freqnum );
Yii = squeeze(YP(1,1,:));
tr_db_i = 20 * log10( Yii );
% Yii2 = zeros( 1,freqnum );
% Yii2(:) = YP(150,150,:);
% tr_db_i2 = 20 * log10( Yii2 );

dis = f_distance( lx, ly, nx, ny );

dis_t = abs( dis - dij1 );
min_num1 = min( dis_t(:) );
fprintf('\n dij1 = %f , error = %f \n',dij1, min_num1 )
[xt,yt] = find( dis_t == min_num1 );
if ~isempty(xt)
%     Yij1 = zeros( 1,freqnum );
    Yij1 = squeeze( YP( xt(1),yt(1),:) );
    tr_db_t1 = 20 * log10( Yij1 );
end

dis_t = abs( dis - dij2 );
min_num2 = min( dis_t(:) );
fprintf('\n dij2 = %f , error = %f \n',dij2, min_num2 )
[xt,yt] = find( dis_t == min_num2 );
if ~isempty(xt)
%     Yij2 = zeros( 1,freqnum );
    Yij2 = squeeze( YP( xt(1),yt(1),:) );
    tr_db_t2 = 20 * log10( Yij2 );
end


figure(1);
plot(freq,tr_db_i,'b');
hold on;
% plot(freq,tr_db_i2,'m');
% figure(2);
plot(freq,(tr_db_t1),'m');
% figure(3);
plot(freq,tr_db_t2,'r');
legend('input PM','dij=0.08','dij=0.84')
xlabel('Frequency (HZ)');
ylabel('Magnitude (dB)');
%plot(freq,tmp_t)
