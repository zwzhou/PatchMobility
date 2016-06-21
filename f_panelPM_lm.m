function [ yp ] = f_panelPM_lm( lx,ly,h, nx,ny, rho_p, freq, omegac_pq_squa, int_patch )
% return patch mobility martix of a panel: patchAmt*patchAmt
%   lx, ly & h: length, width and thickness of the plate
%   nx & ny: patch number along length and width
%   rho_p: density of the plate
%   freq: frequency number
%   omegac_pq_squa: modal frequency, size: Nmodal*1
%   int_patch: intral over patch. size: Nmodal*patchAmt
%      see: test_panelPM.m
% use in loop to save memory

panel_area = lx * ly;
M_pq = rho_p*h*panel_area/4 ;

patchAmt = nx * ny;
delarea = panel_area / patchAmt;

Nmodal = length(omegac_pq_squa);
% modal frequency
% [omega_c_squa,idx_p]  size: Nmodal*1
% [omegac_pq_squa,idx_p] = f_plateOmegaSqua_simspt(lx,ly,Dc,rho_p*h,Nmodal,np,nq);


% freqNum = length(freq);


omega = 2*pi*freq;
omega_squa = omega^2;% 

% integral patch i and j

	
	% yp1 = bsxfun(@minus,omegac_pq_squa,omega_squa);% size: Nmodal*1*freqNUM
	% yp2 = bsxfun(@rdivide,intij,yp1);% size: Nmodal*patchamt*freqNUM
	% % equal to the two lines above.  size: Nmodal*patchamt*freqNUM
	% % yp2 = bsxfun(@rdivide,intij,bsxfun(@minus,omegac_pq_squa,omega_squa));
	% yp3 = reshape(yp2,Nmodal,patchamt*freqNUM);% size: Nmodal*(patchamt*freqNUM); two dimensional maxtrix
	
	% equal to the three lines above.
	% size: Nmodal*(patchamt*freqNUM); two dimensional maxtrix
% 	yp3 = reshape(...
%         bsxfun(@rdivide,int_patch, bsxfun(@minus,omegac_pq_squa,omega_squa) ),...
%         Nmodal,patchAmt*freqNum);
% 	yp = bsxfun( @times, 1i/delarea^2/M_pq*omega ,...
%         reshape(int_patch.'*yp3,patchAmt,patchAmt,freqNum) );
	% clear yp3
	


end

