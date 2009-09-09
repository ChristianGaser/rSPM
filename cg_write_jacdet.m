function VO = cg_write_jacdet(prm,flags,extras)
% Write Out Jacobian Determinant Images.
% FORMAT VO = cg_write_jacdet(V,matname,flags)
% V         - Images to transform (filenames or volume structure).
% matname   - Transformation information (filename or structure).
% flags     - flags structure, with fields...
%           interp   - interpolation method (0-7)
%           wrap     - wrap edges (e.g., [1 1 0] for 2D MRI sequences)
%           vox      - voxel sizes (3 element vector - in mm)
%                      Non-finite values mean use template vox.
%           bb       - bounding box (2x3 matrix - in mm)
%                      Non-finite values mean use template bb.
%
% Images are written prefixed by "j0".
%
% Non-finite vox or bounding box suggests that values should be derived
% from the template image.
%
% Don't use interpolation methods greater than one for data containing
% NaNs.
%
%_______________________________________________________________________
% Copyright (C) 2005 Wellcome Department of Imaging Neuroscience

% John Ashburner
% $Id: cg_write_jacdet.m 771 2007-03-20 14:19:41Z john $
%
% @(#)cg_write_jacdet.m	v1.02 Christian Gaser 2008/06/16

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

[def_flags.bb, def_flags.vox] = bbvox_from_V(prm.VG(1));

if nargin < 3,
	flags = def_flags;
else
	fnms = fieldnames(def_flags);
	for i=1:length(fnms),
		if ~isfield(flags,fnms{i}),
			flags.(fnms{i}) = def_flags.(fnms{i});
		end;
	end;
end;

if ~all(isfinite(flags.vox(:))), flags.vox = def_flags.vox; end;
if ~all(isfinite(flags.bb(:))),  flags.bb  = def_flags.bb;  end;

[x,y,z,mat] = get_xyzmat(prm,flags.bb,flags.vox);

if nargout==0,
		cg_calc_jacdet(prm.VF,prm,x,y,z,mat,flags,pth);
else
		VO = cg_calc_jacdet(prm.VF,prm,x,y,z,mat,flags,pth);
end;

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



