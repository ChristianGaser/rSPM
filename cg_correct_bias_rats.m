function cg_correct_bias_rats
%
% $Id: cg_correct_bias_rats.m 30 2011-09-08 20:00:39Z gaser $

if nargin==0
    P = spm_select(Inf,'image','Select images');
end
V = spm_vol(P);
s0 = [6 6 6];

spm_progress_bar('Init',size(P,1),'remove NaN','volumes completed')
for i = 1:size(P,1)
	vol = spm_read_vols(V(i));
	
    VOX = sqrt(sum(V(i).mat(1:3,1:3).^2));
    s  = s0./VOX;                        % voxel anisotropy
    s1 = s/sqrt(8*log(2));              % FWHM -> Gaussian parameter

    x  = round(6*s1(1)); x = -x:x; x = spm_smoothkern(s(1),x,1); x  = x/sum(x);
    y  = round(6*s1(2)); y = -y:y; y = spm_smoothkern(s(2),y,1); y  = y/sum(y);
    z  = round(6*s1(3)); z = -z:z; z = spm_smoothkern(s(3),z,1); z  = z/sum(z);

    i2  = (length(x) - 1)/2;
    j2  = (length(y) - 1)/2;
    k2  = (length(z) - 1)/2;

    vol2 = zeros(size(vol));
    spm_conv_vol(vol,vol2,x,y,z,-[i2,j2,k2]);

    vol = vol./(vol2 + eps);

    [pth,nm,xt,vr] = spm_fileparts(V(i).fname);
    V2 = V(i);
    V2.fname = fullfile(pth,['m' nm xt vr]);
    V2.dt(1) = 16;
    spm_write_vol(V2,vol);
	spm_progress_bar('Set',i);
end
spm_progress_bar('Clear');