function [ x1, x2, y1, y2 ] = findcoor( patchnumber, nx, delx, dely )
%Input patch number and output its coordinate value
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
%tic
rownumber = ceil( patchnumber / nx ) ;
% sw = patchnumber - (rownumber-1)*nx ;
% switch ( sw )
%     case 0
%         colnumber = nx ;
%     otherwise
%         colnumber = sw;
% end
colnumber = patchnumber - (rownumber-1)*nx ;
% colnumber = mod( patchnumber, nx ) ;
        
x2 = colnumber * delx ;
x1 = x2 - delx ;
y2 = rownumber * dely ;
y1 = y2 - dely ;
%toc

% tic
% for i = 0:ny-1 
%     colnumber = patchnumber - i*nx ;
%     if ( colnumber <= nx && colnumber >= 0 )
%         rownumber = i+1 ;
%         x2 = colnumber * delx ;
%         x1 = x2 - delx ;
%         y2 = rownumber * dely ;
%         y1 = y2 - dely ;
%         toc
%         return
%     end
% end

end

