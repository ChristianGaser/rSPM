function [vol_name, vol_left, vol_right] = cg_volume_paxinos(P)
% Write Out Paxinos ROI's of Jacobian Determinant Images For Longitudinal Data.
% FORMAT [vol_name, vol_left, vol_right] = cg_volume_paxinos_long(data)
%
% % $Id: cg_volume_paxinos_long.m 38 2014-04-09 09:01:05Z gaser $

if isempty(P), return; end;

tbs = spm('tbs');
for i=1:length(tbs)
	if strcmp(lower(tbs(i).name),'rspm')
		rSPMdir = tbs(i).dir; 
	end
end

atlas_label   = fullfile(rSPMdir,'Paxinos_labeled.nii');
brainmask     = fullfile(rSPMdir,'Brainmask-Paxinos-avg176.nii');
paxinos_atlas = fullfile(rSPMdir,'Paxinos_label.txt');

% get bounding box of Paxinos labels
[bb, vox] = bbvox_from_V(spm_vol(atlas_label));

fid = fopen(paxinos_atlas);
if fid > 0
	atlas = textscan(fid,'%s%d%d','Delimiter','\t');
else
	error(sprintf('Could not open file %s.',paxinos_atlas));
end
fclose(fid);

Vlabel = spm_vol(atlas_label);
Vmask = spm_vol(brainmask);
V = spm_vol(P);

label = zeros(Vlabel.dim,'uint8');
JD = zeros(Vlabel.dim);

vol_left  = 0;
vol_right = 0;
mid = round(Vlabel.dim(1)/2);
for j=1:Vlabel.dim(3),
	Mi  = spm_matrix([0 0 j]);
	M1  = Vlabel.mat\Vmask.mat\Mi;
	msk = spm_slice_vol(Vmask,M1,Vlabel.dim(1:2),[1 0]);
	M2  = Vlabel.mat\V.mat\Mi;
	img = spm_slice_vol(V,M2,Vlabel.dim(1:2),[1 0]);
	img = img.*(msk > 0) + 1;
	JD(:,:,j) = img;
	vol_left  = vol_left  + sum(sum(img(1:mid,:)));
	vol_right = vol_right + sum(sum(img(mid+1:Vlabel.dim(1),:)));
end;

vx_vol = abs(prod(vox));
vol_name = strvcat('Hemispheres');
vol_left  = vx_vol*vol_left;
vol_right = vx_vol*vol_right;
 
for j=1:Vlabel.dim(3),
	img = spm_slice_vol(Vlabel,spm_matrix([0 0 j]),Vlabel.dim(1:2),1);
	label(:,:,j) = uint8(img);
end;

for i = 1:length(atlas{1})
	fprintf('.');
	vol_name = strvcat(vol_name, char(atlas{1}(i)));
	vol_left_subject  = 0;
	vol_right_subject = 0;
	for j=1:Vlabel.dim(3)
		vol_left_subject  = vol_left_subject  + vx_vol*sum(sum(JD(:,:,j).*(label(:,:,j)==atlas{2}(i))));
		vol_right_subject = vol_right_subject + vx_vol*sum(sum(JD(:,:,j).*(label(:,:,j)==atlas{3}(i))));
	end
	vol_left  = [vol_left;  vol_left_subject];
	vol_right = [vol_right; vol_right_subject];
end

return;


%_______________________________________________________________________
function [bb,vx] = bbvox_from_V(V)
vx = sqrt(sum(V.mat(1:3,1:3).^2));
if det(V.mat(1:3,1:3))<0, vx(1) = -vx(1); end;

o  = V.mat\[0 0 0 1]';
o  = o(1:3)';
bb = [-vx.*(o-1) ; vx.*(V.dim(1:3)-o)]; 
return;
%_______________________________________________________________________

