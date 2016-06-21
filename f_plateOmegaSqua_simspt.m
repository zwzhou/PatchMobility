function [ omegasqua, idx ] = f_plateOmegaSqua_simspt( lx,ly,D,m,Np,nx,ny )
% square circular frequency of a plate
%   lx,ly : plate size
%   nx,ny : order of x and z
%   D : D=Eh^3/[12(1-v^2)]
%   m : unit mass of the plate
%   Np : return order
%
%   idx : struct 
%       idx.x idx.y : index of the circular frequency
%   
[a,b] = meshgrid(1:nx,1:ny);
omega_tmp = ((a./lx).^2+(b./ly).^2).^2;
[omegasqua_t,i] = sort(reshape(omega_tmp,[],1));
tmp = (pi^4)*D/m;
omegasqua = tmp* omegasqua_t(1:Np);

a = reshape(a,[],1);
a = a(i);
idx.x = a(1:Np);

b = reshape(b,[],1);
b = b(i);
idx.y = b(1:Np);
end

