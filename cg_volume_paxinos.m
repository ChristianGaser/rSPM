function [vol_name, vol_left, vol_right] = cg_volume_paxinos(prm)
% Write Out Jacobian Determinant Images.
% FORMAT [vol_name, vol_left, vol_right] = cg_volume_paxinos(prm)
%
% @(#)cg_volume_paxinos.m	v1.02 Christian Gaser 2008/08/27

if isempty(prm), return; end;

if ischar(prm)
	pth = fileparts(prm);
	prm = load(prm);
else
	pth = '';
end;

if isempty(prm.Tr)
	error('No non-linear component found');
end

def_flags = struct('interp',1,'vox',NaN,'bb',NaN,'wrap',[0 0 0],'preserve',1);

tbs = spm('tbs');
for i=1:length(tbs)
	if strcmp(lower(tbs(i).name),'rspm')
		rSPMdir = tbs(i).dir; 
	end
end

atlas_label   = fullfile(rSPMdir,'Paxinos_labeled.nii');
brainmask     = fullfile(rSPMdir,'Brainmask-Paxinos.nii');
paxinos_atlas = fullfile(rSPMdir,'Paxinos_label.txt');

flags = def_flags;

% get bounding box of Paxinos labels
[flags.bb, flags.vox] = bbvox_from_V(spm_vol(atlas_label));

[x,y,z,mat] = get_xyzmat(prm,flags.bb,flags.vox);

fid = fopen(paxinos_atlas);
if fid > 0
	atlas = textscan(fid,'%s%d%d','Delimiter','\t');
else
	error(sprintf('Could not open file %s.',paxinos_atlas));
end
fclose(fid);

VO = cg_calc_jacdet(prm.VF,prm,x,y,z,mat,flags,pth);

Vmask = spm_vol(brainmask);

vol_left  = 0;
vol_right = 0;
mid = round(VO.dim(1)/2);
for j=1:VO.dim(3),
	Mi  = spm_matrix([0 0 j]);
	M1  = VO.mat\Vmask.mat\Mi;
	msk = spm_slice_vol(Vmask,M1,VO.dim(1:2),[1 0]);
	img = VO.dat(:,:,j).*(msk > 0);
	vol_left  = vol_left  + sum(sum(img(1:mid,:)));
	vol_right = vol_right + sum(sum(img(mid+1:VO.dim(1),:)));
end;

vx_vol = abs(prod(flags.vox));
vol_name = strvcat('Hemispheres');
vol_left  = vx_vol*vol_left;
vol_right = vx_vol*vol_right;

Vlabel = spm_vol(atlas_label);
label = uint8(0);
label(Vlabel.dim(1),Vlabel.dim(2),Vlabel.dim(3)) = uint8(0);
 
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
		vol_left_subject  = vol_left_subject  + vx_vol*sum(sum(VO.dat(:,:,j).*(label(:,:,j)==atlas{2}(i))));
		vol_right_subject = vol_right_subject + vx_vol*sum(sum(VO.dat(:,:,j).*(label(:,:,j)==atlas{3}(i))));
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

%_______________________________________________________________________
function [x,y,z,mat] = get_xyzmat(prm,bb,vox)
% The old voxel size and origin notation is used here.
% This requires that the position and orientation
% of the template is transverse.  It would not be
% straitforward to account for templates that are
% in different orientations because the basis functions
% would no longer be seperable.  The seperable basis
% functions mean that computing the deformation field
% from the parameters is much faster.

% bb  = sort(bb);
% vox = abs(vox);

msk       = find(vox<0);
bb        = sort(bb);
bb(:,msk) = flipud(bb(:,msk));

% Adjust bounding box slightly - so it rounds to closest voxel.
% Comment out if not needed.  I chose not to change it because
% it would lead to being bombarded by questions about spatially
% normalised images not having the same dimensions.
bb(:,1) = round(bb(:,1)/vox(1))*vox(1);
bb(:,2) = round(bb(:,2)/vox(2))*vox(2);
bb(:,3) = round(bb(:,3)/vox(3))*vox(3);

M   = prm.VG(1).mat;
vxg = sqrt(sum(M(1:3,1:3).^2));
if det(M(1:3,1:3))<0, vxg(1) = -vxg(1); end;
ogn = M\[0 0 0 1]';
ogn = ogn(1:3)';

% Convert range into range of voxels within template image
x   = (bb(1,1):vox(1):bb(2,1))/vxg(1) + ogn(1);
y   = (bb(1,2):vox(2):bb(2,2))/vxg(2) + ogn(2);
z   = (bb(1,3):vox(3):bb(2,3))/vxg(3) + ogn(3);

og  = -vxg.*ogn;

% Again, chose whether to round to closest voxel.
of  = -vox.*(round(-bb(1,:)./vox)+1);
%of = bb(1,:)-vox;

M1  = [vxg(1) 0 0 og(1) ; 0 vxg(2) 0 og(2) ; 0 0 vxg(3) og(3) ; 0 0 0 1];
M2  = [vox(1) 0 0 of(1) ; 0 vox(2) 0 of(2) ; 0 0 vox(3) of(3) ; 0 0 0 1];
mat = prm.VG(1).mat*inv(M1)*M2;

LEFTHANDED = true;
if (LEFTHANDED && det(mat(1:3,1:3))>0) || (~LEFTHANDED && det(mat(1:3,1:3))<0),
	Flp = [-1 0 0 (length(x)+1); 0 1 0 0; 0 0 1 0; 0 0 0 1];
	mat = mat*Flp;
	x   = flipud(x(:))';
end;
return;
%_______________________________________________________________________


