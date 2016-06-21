function [ patch ] = f_ploc( patchamt, nx, delx, dely )
% Input patch amount and output all patch coordinate value. output size: 1*patchamt
%   nx and ny: patch number at x and y axis;  delx and dely: patch size 
%   out put:
%
%           y ^
%             |          y2  __________
%             |             |          |
%             |             |  patch   |
%             |             |  number  |
%             |          y1 |__________|    
%             |            x1         x2   
%             |_____                             |
%             |nx+1 |                            | 
%             |_____|                       _____|
%             |  1  |                      |  nx |
%             |_____|______________________|_____|___\>
%             0                                        x
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                    20141229
patchnumber = 1:patchamt;
rownumber = ceil( patchnumber ./ nx ) ;
colnumber = patchnumber - (rownumber-1).*nx ;
patch.x2 = colnumber .* delx ;
patch.x1 = patch.x2 - delx ;
patch.y2 = rownumber .* dely ;
patch.y1 = patch.y2 - dely ;
end

