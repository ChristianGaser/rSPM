function out = spm_hdw(job)
% Warp a pair of same subject images together
%
% Very little support can be provided for the warping routine, as it
% involves optimising a very nonlinear objective function.
% Also, don't ask what the best value for the regularisation is.
%_______________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% modified version of spm_warp from
% John Ashburner
% $Id: cg_hdw.m 27 2011-06-21 15:00:15Z gaser $

for i=1:numel(job.subj),
    out(i).files = cell(numel(job.subj(i).mov)-1,1);
    for j=2:numel(job.subj(i).mov)
        [pth,nam,ext,num] = spm_fileparts(job.subj(i).mov{j});
        if isempty(char(job.warp_opts.weight))
            out(i).files{j-1} = fullfile(pth,['jd_y_', nam, ext, num]);
        else
            out(i).files{j-1} = fullfile(pth,['jd_my_', nam, ext, num]);
        end
        run_warping(job.subj(i).mov{j},job.subj(i).mov{1},job.warp_opts,job.bias_opts);
    end
end;
%_______________________________________________________________________

%_______________________________________________________________________
function run_warping(PF,PG,warp_opts,bias_opts)

VG  = spm_vol(PG);
VF  = spm_vol(PF);
if isempty(char(warp_opts.subsamp))
  VG0 = VF;
else
  VG0 = spm_vol(char(warp_opts.subsamp));
end
mask  = warp_opts.weight;
reg = warp_opts.reg;
nit = warp_opts.nits;

% Load the image pair as 8 bit.
%-----------------------------------------------------------------------
bo            = bias_opts;
[VG,sf]       = loadfloat(VG,VG0,mask); % Template
VF      = bias_correction(VF,VG,VG0,mask,bo.nits,bo.fwhm,bo.reg,bo.lmreg,sf);

% Try loading pre-existing deformation fields.  Otherwise, create
% deformation fields from uniform affine transformations.
%-----------------------------------------------------------------------
[pth,nme,ext,num] = spm_fileparts(VF.fname);
ofname = fullfile(pth,['y_' nme '.nii']);
Def    = cell(3,1);
fprintf('Generating uniform affine transformation field\n');
Def{1} = single(1:VG.dim(1))';
Def{1} = Def{1}(:,ones(VG.dim(2),1),ones(VG.dim(3),1));
Def{2} = single(1:VG.dim(2));
Def{2} = Def{2}(ones(VG.dim(1),1),:,ones(VG.dim(3),1));
Def{3} = reshape(single(1:VG.dim(3)),1,1,VG.dim(3));
Def{3} = Def{3}(ones(VG.dim(1),1),ones(VG.dim(2),1),:);
spm_affdef(Def{:},VF.mat\VG.mat);
inverse = 0;
if inverse	
	% inverse field
  iDef    = cell(3,1);
  fprintf('Generating inverse affine transformation field\n');
  iDef{1} = single(1:VG.dim(1))';
  iDef{1} = iDef{1}(:,ones(VG.dim(2),1),ones(VG.dim(3),1));
  iDef{2} = single(1:VG.dim(2));
  iDef{2} = iDef{2}(ones(VG.dim(1),1),:,ones(VG.dim(3),1));
  iDef{3} = reshape(single(1:VG.dim(3)),1,1,VG.dim(3));
  iDef{3} = iDef{3}(ones(VG.dim(1),1),ones(VG.dim(2),1),:);
 spm_affdef(iDef{:},VG.mat\VF.mat);
end

% Voxel sizes
%-----------------------------------------------------------------------
vxg = sqrt(sum(VG.mat(1:3,1:3).^2))';if det(VG.mat(1:3,1:3))<0, vxg(1) = -vxg(1); end;
vxf = sqrt(sum(VF.mat(1:3,1:3).^2))';if det(VF.mat(1:3,1:3))<0, vxf(1) = -vxf(1); end;

% Do warping
%-----------------------------------------------------------------------
fprintf('Warping (iterations=%d regularisation=%g)\n', nit, reg);

cg_warp(VG.single,VF.single,Def{:},[vxg vxf],[nit,reg,1,0]);

if inverse
	% inverse warping
	cg_warp(VF.single,VG.single,iDef{:},[vxg vxf],[nit,reg,1,0]);

	% combine both by averaging
	try
		[iDef{1},iDef{2},iDef{3}] = spm_invdef(iDef{1},iDef{2},iDef{3},VG.dim(1:3),VG.mat\VF.mat,VF.mat\VG.mat);
		for i=1:3
		    Def{i} = (Def{i} + iDef{i})/2;
		end
	catch
		warning('spm_invdef failed: No inverse warping is used.')
	end
end

% Convert mapping from voxels to mm
%-----------------------------------------------------------------------
spm_affdef(Def{:},VF.mat);

% Write the deformations
%-----------------------------------------------------------------------
save_def(Def,VG.mat,ofname,mask)

return;
%_______________________________________________________________________

function [VO,sf] = loadfloat(V,VG0,mask)
% Load data from file indicated by V into an array of floating values.

if nargin<3, mask=''; end
if nargin<2, VG0=V; end

if isempty(VG0), VG0 = V; end

if ~isempty(char(mask))
  fprintf('Apply mask to %s\n',V.fname);
  Vm = spm_vol(char(mask));
end

M = VG0.mat;

spm_progress_bar('Init',V.dim(3),...
                ['Computing max/min of ' spm_str_manip(V.fname,'t')],...
                'Planes complete');
mx = -Inf;
for p=1:V.dim(3),
    img = spm_slice_vol(V,spm_matrix([0 0 p]),V.dim(1:2),1);
    mx  = max([max(img(:)) mx]);
    spm_progress_bar('Set',p);
end;
spm_progress_bar('Init',V.dim(3),...
        ['Loading ' spm_str_manip(V.fname,'t')],...
        'Planes loaded');
sf = 1;

udat = zeros(VG0.dim(1:3),'single');
for p=1:VG0.dim(3),
        M1 = M\V.mat\spm_matrix([0 0 p]);
        img = spm_slice_vol(V,M1,VG0.dim(1:2),1);
        udat(:,:,p) = single(img*sf);
        % load mask image
        if ~isempty(char(mask))
            Mm = M\Vm.mat\spm_matrix([0 0 p]);
            tmp_mask = spm_slice_vol(Vm,Mm,VG0.dim(1:2),1);
            udat(:,:,p) = udat(:,:,p).*single(tmp_mask>0);
        end
        spm_progress_bar('Set',p);
end;
spm_progress_bar('Clear');
VO = VG0;
VO.fname = V.fname;
VO.single = udat;
return;
%_______________________________________________________________________

%_______________________________________________________________________
function t = transf(B1,B2,B3,T)
d2 = [size(T) 1];
t1 = reshape(reshape(T, d2(1)*d2(2),d2(3))*B3', d2(1), d2(2));
t  = B1*t1*B2';
return;
%_______________________________________________________________________

%_______________________________________________________________________
function VO = bias_correction(VF,VG,VG0,mask,nits,fwhm,reg1,reg2,sf)
% This function is intended for doing bias correction prior
% to high-dimensional intra-subject registration.  Begin by
% coregistering, then run the bias correction, and then do
% the high-dimensional warping.  A version of the first image
% is returned out, that has the same bias as that of the second
% image.  This allows warping of the corrected first image, to
% match the uncorrected second image.
%
% Ideally, it would all be combined into the same model, but
% that would require quite a lot of extra work.  I should really
% sort out a good convergence criterion, and make it a proper
% Levenberg-Marquardt optimisation.

if isempty((VG0)), VG0 = VF; end

if ~isempty(char(mask))
  fprintf('Apply mask to %s\n',VF.fname);
  Vm = spm_vol(char(mask));
end

M = VG0.mat;

vx      = sqrt(sum(VG0.mat(1:3,1:3).^2));
d       = VG0.dim(1:3);
sd      = vx(1)*VG0.dim(1)/fwhm; d3(1) = ceil(sd*2); krn_x   = exp(-(0:(d3(1)-1)).^2/sd.^2)/sqrt(vx(1));
sd      = vx(2)*VG0.dim(2)/fwhm; d3(2) = ceil(sd*2); krn_y   = exp(-(0:(d3(2)-1)).^2/sd.^2)/sqrt(vx(2));
sd      = vx(3)*VG0.dim(3)/fwhm; d3(3) = ceil(sd*2); krn_z   = exp(-(0:(d3(3)-1)).^2/sd.^2)/sqrt(vx(3));

Cbias   = kron(krn_z,kron(krn_y,krn_x)).^(-2)*prod(d)*reg1;
Cbias   = sparse(1:length(Cbias),1:length(Cbias),Cbias,length(Cbias),length(Cbias));
B3bias  = spm_dctmtx(d(3),d3(3));
B2bias  = spm_dctmtx(d(2),d3(2));
B1bias  = spm_dctmtx(d(1),d3(1));
lmRb    = speye(size(Cbias))*prod(d)*reg2;
Tbias   = zeros(d3);

ll = Inf;
spm_chi2_plot('Init','Bias Correction','- Log-likelihood','Iteration');
for subit=1:nits,

    % Compute objective function and its 1st and second derivatives
    Alpha = zeros(prod(d3),prod(d3)); % Second derivatives
    Beta  = zeros(prod(d3),1); % First derivatives
    oll   = ll;
    ll    = 0.5*Tbias(:)'*Cbias*Tbias(:);

    for z=1:VG0.dim(3),
        M1 = M\VF.mat\spm_matrix([0 0 z]);
        f1o = spm_slice_vol(VF,M1,VG0.dim(1:2),0);
        f2o = double(VG.single(:,:,z))/sf;
        % load mask image
        if ~isempty(char(mask))
            Mm = M\Vm.mat\spm_matrix([0 0 z]);
            tmp_mask = spm_slice_vol(Vm,Mm,VG0.dim(1:2),1);
            f1o = f1o.*double(tmp_mask>0);
            f2o = f2o.*double(tmp_mask>0);
        end
        msk = (f1o==0) & (f2o==0);
        f1o(msk) = 0;
        f2o(msk) = 0;
        ro       = transf(B1bias,B2bias,B3bias(z,:),Tbias);
        msk      = abs(ro)>0.01; % false(d(1:2));

        % Use the form based on an integral for bias that is
        % far from uniform.
        f1  = f1o(msk);
        f2  = f2o(msk);
        r   = ro(msk);
        e   = exp(r);
        t1  = (f2.*e-f1);
        t2  = (f1./e-f2);
        ll  = ll + 1/4*sum(sum((t1.^2-t2.^2)./r));
        wt1 = zeros(size(f1o));
        wt2 = zeros(size(f1o));
        wt1(msk) = (2*(t1.*f2.*e+t2.*f1./e)./r + (t2.^2-t1.^2)./r.^2)/4;
        wt2(msk) = ((f2.^2.*e.^2-f1.^2./e.^2+t1.*f2.*e-t2.*f1./e)./r/2 ...
                 - (t1.*f2.*e+t2.*f1./e)./r.^2 + (t1.^2-t2.^2)./r.^3/2);

        % Use the simple symmetric form for bias close to uniform
        f1  = f1o(~msk);
        f2  = f2o(~msk);
        r   = ro(~msk);
        e   = exp(r);
        t1  = (f2.*e-f1);
        t2  = (f1./e-f2);
        ll  = ll + (sum(t1.^2)+sum(t2.^2))/4;
        wt1(~msk) = (t1.*f2.*e-t2.*f1./e)/2;
        wt2(~msk) = ((f2.*e).^2+t1.*f2.*e + (f1./e).^2+t2.*f1./e)/2;

        b3    = B3bias(z,:)';
        Beta  = Beta  + kron(b3,spm_krutil(wt1,B1bias,B2bias,0));
        Alpha = Alpha + kron(b3*b3',spm_krutil(wt2,B1bias,B2bias,1));
    end;
    spm_chi2_plot('Set',ll/prod(d));


    if subit > 1 && ll>oll,
        % Hasn't improved, so go back to previous solution
        Tbias = oTbias;
        ll    = oll;
        lmRb  = lmRb*10;
    else
        % Accept new solution
        oTbias = Tbias;
        Tbias  = Tbias(:);
        Tbias  = Tbias - (Alpha + Cbias + lmRb)\(Beta + Cbias*Tbias);
        Tbias  = reshape(Tbias,d3);
    end;
end;

udat = zeros(VG0.dim(1:3),'single');
for z=1:VG0.dim(3),
    M1 = M\VF.mat\spm_matrix([0 0 z]);
    f1 = spm_slice_vol(VF,M1,VG0.dim(1:2),1);
    r  = transf(B1bias,B2bias,B3bias(z,:),Tbias);
    f1 = f1./exp(r);
    udat(:,:,z) = single(f1*sf);
    % load mask image
    if ~isempty(char(mask))
        Mm = M\Vm.mat\spm_matrix([0 0 z]);
        tmp_mask = spm_slice_vol(Vm,Mm,VG0.dim(1:2),1);
        udat(:,:,z) = udat(:,:,z).*single(tmp_mask>0);
    end
end;

VO = VG0;
VO.fname = VF.fname;
VO.single = udat;

return;

%_______________________________________________________________________
function save_def(Def,mat,fname,mask)
% Save a deformation field as an image

dim   = [size(Def{1},1) size(Def{1},2) size(Def{1},3) 1 3];
dtype = 'FLOAT32-BE';
off   = 0;
scale = 1;
inter = 0;
dat   = file_array(fname,dim,dtype,off,scale,inter);

N      = nifti;
N.dat  = dat;
N.mat  = mat;
N.mat0 = mat;
N.mat_intent  = 'Aligned';
N.mat0_intent = 'Aligned';
N.intent.code = 'VECTOR';
N.intent.name = 'Mapping';
N.descrip     = 'Deformation field';
create(N);
N.dat(:,:,:,1,1) = Def{1};
N.dat(:,:,:,1,2) = Def{2};
N.dat(:,:,:,1,3) = Def{3};

% Write Jacobian determinants
[pth,nme,ext,num] = spm_fileparts(fname);
jfname = fullfile(pth,['jd_' nme '.nii']);
if nargin > 3
  if ~isempty(char(mask))
    jfname = fullfile(pth,['jd_m' nme '.nii']);
  end
end
dt     = spm_def2det(Def{:},N.mat);
dt = double(dt);
dt(dt~=0) = dt(dt~=0) - 1;		
dt(isnan(dt)) = 0;
dat    = file_array(jfname,dim(1:3),dtype,off,scale,inter);
N      = nifti;
N.dat  = dat;
N.mat  = mat;
N.mat0 = mat;
N.mat_intent  = 'Aligned';
N.mat0_intent = 'Aligned';
N.intent.code = 'NONE';
N.intent.name = 'Det';
N.descrip     = 'Jacobian Determinant - 1';
create(N);
N.dat(:,:,:) = dt;
return;
%_______________________________________________________________________
