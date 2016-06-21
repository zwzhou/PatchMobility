function [ frequency ] = platefreq( lx, ly, h, fn )
% caculate frequency of a plate
%   Detailed explanation goes here

rho = 7.800e-9;% Mkg/cm^3
nu = 0.3; % Poisson ratio
E = 2.06e5; % Pa   Young's modulus
D = E*h^3/(12*(1-nu^2));
nl=fn;   nm=fn;
total = (nl)*(nm);
freq = ones(total,1);
cl = ones(total,1);
cm = ones(total,1);
i = 0;
fr = sqrt(D/(rho*h)) * pi^2;
for l=1:nl
    for m=1:nm
        i = i+1;
        freq(i) = fr*((l/lx)^2+(m/ly)^2);
        cl(i) = l;
        cm(i) = m;
    end
end
[freq0,ix] = sort(freq/(2*pi));
frequency = freq0(1:fn);
sl = cl(ix);
sm = cm(ix);
fprintf(' order   frequency    nx    ny  \n');
for j=1:fn
    %disp([j frequency(j) sl(j) sm(j) sn(j)])
    if ~mod(j,15)
        fprintf('\n order   frequency    nx    ny  \n');
    end
    fprintf('%5i   %10.6f %5i %5i\n',j,frequency(j), sl(j), sm(j) );
end
[maxl,indexl] = max(sl(1:fn));
[maxm,indexm] = max(sm(1:fn));
fprintf('\n max nx = %3i , in order %i \n',maxl,indexl);
fprintf(' max ny = %3i , in order %i \n',maxm,indexm);


end