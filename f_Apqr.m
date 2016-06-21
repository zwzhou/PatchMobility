function [ A_pqr,idx_c ] = f_Apqr( lx,ly,lz, Xs,Ys,Zs,S_0, c,freq,rho_a,Nmodal,np,nq,nr  )
% return Apqr in function 41 42 and 43
%   lx,ly,lz: Cavity size
%   Xs,Ys,Zs: Source location
%   S_0: Source amplitude
%   Lx,Lz: panel size
%   c: sound speed
%   freq: frequency vector.  size: 1*freqNUM
%   Nmodal: total frequency order
%   np,nq,nr: frequency order in x,y and z axis.
%       !!  Nmodal must less than np*nq*nr  !!

fprintf('    ---> Calculating Apqr\n')

fprintf('    --->     step %3d/2  Modal Frequency\n',1)

% np = 350;
% nq = 250;
% nr = 150;
% Nmodal = 100000 ;
% [omegasqua,idx_c] = f_cavityfreqsqua(lx,ly,lz,c,Nmodal,np,nq,nr);
[kpqr_squa,idx_c] = f_cavityksqua(lx,ly,lz,Nmodal,np,nq,nr);

fprintf('    --->     step %3d/2  Apqr\n',2)

intphis = S_0 * cos(idx_c.x*pi*Xs/lx).*cos(idx_c.y*pi*Ys/ly).*cos(idx_c.z*pi*Zs/lz);

lxcomp(1:Nmodal,1) = lx;
lxcomp(idx_c.x~=0) = 0.5*lx;
lycomp(1:Nmodal,1) = ly;
lycomp(idx_c.y~=0) = 0.5*ly;
lzcomp(1:Nmodal,1) = lz;
lzcomp(idx_c.z~=0) = 0.5*lz;
% Npqr = repmat(lxcomp.*lycomp.*lzcomp,1,freqNUM);
Npqr = lxcomp.*lycomp.*lzcomp;

% k_pqr_squa = repmat(omegasqua/c^2,1,freqNUM);
% k_pqr_squa = omegasqua/c^2;

T_r = 10;
eta_r = 2.2./(freq*T_r);
% eta_r = 0.01;
cc = c* sqrt(1 + 1i*eta_r);
% kc_squa = repmat( (2*pi*freq./cc).^2, Nmodal,1);
kc_squa = (2*pi*freq./cc).^2;


% A_pqr   size: Nmodal*freqNUM
% A_pqr = bsxfun( @rdivide,intphis./Npqr,(bsxfun(@minus,kc_squa,kpqr_squa)) );
omega = 2*pi*freq;
A_pqr = bsxfun(@times,-1i*omega*rho_a,intphis./Npqr)./bsxfun(@minus,kc_squa,kpqr_squa) ;

disp('    ---> Complete')

end

