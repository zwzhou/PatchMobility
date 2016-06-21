function [ ksqua,idx ] = f_cavityksqua( lx,ly,lz,Nc,nx,ny,nz )
% square wave number of a cavity ; matrix: Nc*1
%   lx,ly,lz : cavity size
%   nx,ny,nz : order of x, y and z

%   Nc : return order of the cavity
%
%   ksqua : Nc*1 matrix; k squared
%   idx : struct 
%       idx.x idx.y idx.z : Nc*1 matrix; index of the circular frequency

[a,b,d] = meshgrid(0:nx,0:ny,0:nz);
k_tmp =  ((a./lx).^2+(b./ly).^2+(d./lz).^2);% all frequency
% pi^2/(2*pi) == pi/2
[ksqua_t,i] = sort(reshape(k_tmp,[],1));
ksqua = pi^2*ksqua_t(1:Nc);

a = reshape(a,[],1);
a = a(i);
idx.x = a(1:Nc);
b = reshape(b,[],1);
b = b(i);
idx.y = b(1:Nc);
d = reshape(d,[],1);
d = d(i);
idx.z = d(1:Nc);

end

