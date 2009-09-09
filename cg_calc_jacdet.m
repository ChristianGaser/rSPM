function VO = cg_calc_jacdet(V,prm,x,y,z,mat,flags,pth)

% @(#)cg_calc_jacdet.m	v1.02 Christian Gaser 2008/06/16

[X,Y] = ndgrid(x,y);
Tr = prm.Tr;
BX = spm_dctmtx(prm.VG(1).dim(1),size(Tr,1),x-1);
BY = spm_dctmtx(prm.VG(1).dim(2),size(Tr,2),y-1);
BZ = spm_dctmtx(prm.VG(1).dim(3),size(Tr,3),z-1);
if flags.preserve,
	DX = spm_dctmtx(prm.VG(1).dim(1),size(Tr,1),x-1,'diff');
	DY = spm_dctmtx(prm.VG(1).dim(2),size(Tr,2),y-1,'diff');
	DZ = spm_dctmtx(prm.VG(1).dim(3),size(Tr,3),z-1,'diff');
end;
d  = [flags.interp*[1 1 1]' flags.wrap(:)];

for i=1:numel(V),
	VO     = make_hdr_struct(V(i),x,y,z,mat);
	detAff = det(prm.VF.mat*prm.Affine/prm.VG(1).mat);

	if flags.preserve || nargout>0,
		Dat = single(0);
		Dat(VO.dim(1),VO.dim(2),VO.dim(3)) = single(0);
	else
		VO  = spm_create_vol(VO);
	end;

	for j=1:length(z),   % Cycle over planes
		% Nonlinear deformations
		%----------------------------------------------------------------------------
		tx = get_2Dtrans(Tr(:,:,:,1),BZ,j);
		ty = get_2Dtrans(Tr(:,:,:,2),BZ,j);
		tz = get_2Dtrans(Tr(:,:,:,3),BZ,j);

		if ~flags.preserve,
			if nargout>0,
				Dat(:,:,j) = single(dat);
			else
				VO = spm_write_plane(VO,dat,j);
			end;
		else
			j11 = DX*tx*BY' + 1; j12 = BX*tx*DY';     j13 = BX*get_2Dtrans(Tr(:,:,:,1),DZ,j)*BY';
			j21 = DX*ty*BY';     j22 = BX*ty*DY' + 1; j23 = BX*get_2Dtrans(Tr(:,:,:,2),DZ,j)*BY';
			j31 = DX*tz*BY';     j32 = BX*tz*DY';     j33 = BX*get_2Dtrans(Tr(:,:,:,3),DZ,j)*BY' + 1;

			% The determinant of the Jacobian reflects relative volume changes.
			%------------------------------------------------------------------
			dat       = (j11.*(j22.*j33-j23.*j32) - j21.*(j12.*j33-j13.*j32) + j31.*(j12.*j23-j13.*j22)) * detAff;
			Dat(:,:,j) = single(dat);
		end;
	end;
	if nargout==0,
		if flags.preserve,
			VO = spm_write_vol(VO,Dat);
		end;
	else
		VO.pinfo  = [1 0]';
		VO.dt     = [spm_type('float32') spm_platform('bigend')];
		VO.dat    = Dat;
	end;
end;
return;
%_______________________________________________________________________

%_______________________________________________________________________
function VO   = make_hdr_struct(V,x,y,z,mat)
ind = findstr(V.fname,'_sn.mat');
if ~isempty(ind)
	[pth,nm,xt,vr] = fileparts(V.fname);
	V.fname        = fullfile(pth,[nm(ind:end,:) '.img' vr]);
end
VO            = V;
VO.fname      = prepend(V.fname,'jw');
VO.mat        = mat;
VO.dim(1:3)   = [length(x) length(y) length(z)];
VO.pinfo      = [1 0]';
VO.dt         = [spm_type('float32') spm_platform('bigend')];
VO.descrip    = 'rSPM - Jacobian determinant';
return;
%_______________________________________________________________________

%_______________________________________________________________________
function T2 = get_2Dtrans(T3,B,j)
d   = [size(T3) 1 1 1];
tmp = reshape(T3,d(1)*d(2),d(3));
T2  = reshape(tmp*B(j,:)',d(1),d(2));
return;
%_______________________________________________________________________

%_______________________________________________________________________
function PO = prepend(PI,pre)
[pth,nm,xt,vr] = fileparts(deblank(PI));
PO             = fullfile(pth,[pre nm xt vr]);
return;
%_______________________________________________________________________
