function [ coor ] = f_findcoor( patchNumber, nx, delx, dely )
%Input patch number and output its coordinate value
% return coor.x1 coor.x2 coor.y1 and coor.y2 coor.dis
% coor.dis: size patchNumber*patchNumber
% matrix size: all returned items identical to patchNumber: patchNmuber*1
% e.g. if patchNumber size is 1 by n (recommend), coor.x1 coor.x2 coor.y1 
% and coor.y2 size are 1 by n
%   nx and ny: patch number at x and y axis; ( ny not used )
%   delx and dely: patch size 
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
rowNumber = ceil( patchNumber ./ nx ) ;
% sw = patchnumber - (rownumber-1)*nx ;
% switch ( sw )
%     case 0
%         colnumber = nx ;
%     otherwise
%         colnumber = sw;
% end
colNumber = patchNumber - (rowNumber-1).*nx ;
% colnumber = mod( patchnumber, nx ) ;
        
coor.x2 = colNumber .* delx ;
coor.x1 = coor.x2 - delx ;
coor.y2 = rowNumber .* dely ;
coor.y1 = coor.y2 - dely ;

coor.xc = coor.x2 - 0.5 * delx;
coor.yc = coor.y2 - 0.5 * dely;
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

xdis = bsxfun(@minus,coor.x2.',coor.x2);% patch x distance
ydis = bsxfun(@minus,coor.y2.',coor.y2);% patch y distance
coor.dis = sqrt( xdis.^2 + ydis.^2 );% patch distance


end

