function [ coor ] = f_findcoorL( patchnumber, nx, delx, dely, Ax, Ay )
%Input patch number and output its coordinate value
%   nx and ny: patch number at x and y axis;  delx and dely: patch size 
%   out put:
%
%    y ^ 
%      |
%      |
%      |
%      |____________________________________________________
%      |                                                    |
%      |                                                    |
%      |                                                    |           
%      |       __________________________________           |           
%      |      |                                  |          |           
%      |      |          y2  __________          |          |           
%      |      |             |          |         |          |           
%      |      |             |  patch   |         |          |           
%      |      |             |  number  |         |          |           
%      |      |          y1 |__________|         |          |           
%      |      |            x1         x2         |          |           
%      |      |_____                             |          |           
%      |      |nx+1 |                            |          |           
%      |      |_____|                       _____|          |           
%      |      |  1  |                      |  nx |          |           
%      |      |_____|______________________|_____| _ _      |           
%      |  Ax                                        |       |           
%      |______|                                     | Ay    |           
%      |                                            |       |           
%      |____________________________________________|_______|_____>     
%      0                                                          x     
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%tic
rownumber = ceil( patchnumber ./ nx ) ;
% sw = patchnumber - (rownumber-1)*nx ;
% switch ( sw )
%     case 0
%         colnumber = nx ;
%     otherwise
%         colnumber = sw;
% end
colnumber = patchnumber - (rownumber-1).*nx ;
% colnumber = mod( patchnumber, nx ) ;
        
coor.x2 = Ax + colnumber .* delx ;
coor.x1 = coor.x2 - delx ;
coor.y2 = Ay + rownumber .* dely ;
coor.y1 = coor.y2 - dely ;
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

