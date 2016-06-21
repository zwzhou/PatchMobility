function [ row, column ] = findlocation( patchnumber, nx, ny )
%Input patch number and output its location
%   
for i = 0:ny-1 
    colnumber = patchnumber - i*nx ;
    if ( colnumber <= nx && colnumber >= 0 )
        row = i+1 ;
        column = colnumber ;
        return
    end

end

end

