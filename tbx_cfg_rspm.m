function opts = cg_config_rSPM
% Configuration file for normalise jobs
%_______________________________________________________________________
% Copyright (C) 2005 Wellcome Department of Imaging Neuroscience

% based on John Ashburners version of tbx_cfg_rspm.m
% $Id$

cg_rSPM_defaults

%addpath(fileparts(which(mfilename)));

addpath(fullfile(spm('dir'),'toolbox','HDW'));

%------------------------------------------------------------------------

subsamp = cfg_files;
subsamp.name = 'Reference Image for subsampling';
subsamp.tag  = 'subsamp';
subsamp.filter = 'image';
subsamp.ufilter = '^.*';
subsamp.num  = [0 1];
subsamp.def     = @(val)cg_rSPM_get_defaults('hdw.subsamp', val{:});
subsamp.help   = {[...
'This is the reference image, which is used to subsample data to a smaller size']};
%------------------------------------------------------------------------

weight = cfg_files;
weight.name = 'Mask image';
weight.tag  = 'weight';
weight.filter = 'image';
weight.ufilter = '^.*';
weight.num  = [0 1];
weight.def     = @(val)cg_rSPM_get_defaults('hdw.weight', val{:});
weight.help   = {[...
'This is the mask image to limit warpings to values inside mask.']};
%------------------------------------------------------------------------

mov = cfg_files;
mov.name = 'Images to warp';
mov.tag  = 'mov';
mov.filter = 'image';
mov.ufilter = '^rw.*\.nii$';
mov.num  = [1 Inf];
mov.help   = {[...
'These are the images, which are warped to the first image.']};
%------------------------------------------------------------------------

subj = cfg_branch;
subj.name = 'Subject';
subj.tag = 'subj';
subj.val = {mov};
subj.help = {[...
'Images of the same subject, which are to be registered together.  Prior to nonlinear high-dimensional warping, the images should be rigidly registered with each other.']};

%------------------------------------------------------------------------

esubjs         = cfg_repeat;
esubjs.tag     = 'esubjs';
esubjs.name    = 'Data';
esubjs.values  = {subj };
esubjs.num     = [1 Inf];
esubjs.help = {[...
'Specify pairs of images to match together.']};

%------------------------------------------------------------------------

nits = cfg_entry;
nits.name = 'Iterations for bias correction';
nits.tag  = 'nits';
nits.strtype = 'n';
nits.num  = [1 1];
nits.def     = @(val)cg_rSPM_get_defaults('hdw.nits_bias', val{:});
nits.help = {'Number of iterations for the bias correction.'};

%------------------------------------------------------------------------

biasfwhm = cfg_menu;
biasfwhm.name = 'Bias FWHM';
biasfwhm.tag  = 'fwhm';
biasfwhm.labels = {...
'3mm cutoff','4mm cutoff','5mm cutoff','6mm cutoff','7mm cutoff',...
'8mm cutoff','9mm cutoff','10mm cutoff','11mm cutoff','12mm cutoff',...
'13mm cutoff','14mm cutoff','15mm cutoff','No correction'};
biasfwhm.values = {3,4,5,6,7,8,9,10,11,12,13,14,15,Inf};
biasfwhm.def     = @(val)cg_rSPM_get_defaults('hdw.biasfwhm', val{:});
biasfwhm.help = {[...
'FWHM of Gaussian smoothness of bias. If your intensity nonuniformity is very smooth, then choose a large FWHM. This will prevent the algorithm from trying to model out intensity variation due to different tissue types. The model for intensity nonuniformity is one of i.i.d. Gaussian noise that has been smoothed by some amount, before taking the exponential. Note also that smoother bias fields need fewer parameters to describe them. This means that the algorithm is faster for smoother intensity nonuniformities.']};

%------------------------------------------------------------------------

biasreg = cfg_menu;
biasreg.name = 'Bias regularisation';
biasreg.tag = 'reg';
biasreg.labels = {...
'no regularisation','extremely light regularisation',...
'very light regularisation','light regularisation',...
'medium regularisation','heavy regularisation',...
'very heavy regularisation','extremely heavy regularisation'};
biasreg.values = {0, 1e-9, 1e-8, 1e-7, 1e-6, 1e-5, 1e-4, 1e-3};
biasreg.def     = @(val)cg_rSPM_get_defaults('hdw.biasreg', val{:});
biasreg.help = {[...
'We know a priori that intensity variations due to MR physics tend to be spatially smooth, whereas those due to different tissue types tend to contain more high frequency information. A more accurate estimate of a bias field can be obtained by including prior knowledge about the distribution of the fields likely to be encountered by the correction algorithm. For example, if it is known that there is little or no intensity non-uniformity, then it would be wise to penalise large values for the intensity nonuniformity parameters. This regularisation can be placed within a Bayesian context, whereby the penalty incurred is the negative logarithm of a prior probability for any particular pattern of nonuniformity.']};

%------------------------------------------------------------------------

lmreg = cfg_menu;
lmreg.name = 'Levenberg-Marquardt regularisation';
lmreg.tag = 'lmreg';
lmreg.labels = {...
'no regularisation','extremely light regularisation',...
'very light regularisation','light regularisation',...
'medium regularisation','heavy regularisation',...
'very heavy regularisation','extremely heavy regularisation'};
lmreg.values = {0, 1e-9, 1e-8, 1e-7, 1e-6, 1e-5, 1e-4, 1e-3};
lmreg.def     = @(val)cg_rSPM_get_defaults('hdw.lmreg', val{:});
lmreg.help = {[...
'Levenberg-Marquardt regularisation keeps the bias correction part stable. Higher values means more stability, but slower convergence.']};

%------------------------------------------------------------------------

bias_opts = cfg_branch;
bias_opts.name = 'Bias Correction Options';
bias_opts.tag = 'bias_opts';
bias_opts.val = {nits,biasfwhm,biasreg,lmreg};
bias_opts.help = {[...
'MR images are usually corrupted by a smooth, spatially varying artifact that modulates the intensity of the image (bias). These artifacts, although not usually a problem for visual inspection, can impede automated processing of the images.'],...
[...
'Before registering the images, an approximate bias correction is estimated for the moved image. This is based on minimising the difference between the images an a symmetric way. Prior to registering the images, they should be rigidly aligned together.  The bias correction is estimated once for these aligned images.']};

%------------------------------------------------------------------------

warpreg = cfg_entry;
warpreg.name = 'Warping regularisation';
warpreg.tag  = 'reg';
warpreg.strtype = 'e';
warpreg.num  = [1 1];
warpreg.def     = @(val)cg_rSPM_get_defaults('hdw.warpreg', val{:});
warpreg.help = {[...
'There is a tradeoff between the smoothness of the estimated warps, and the difference between the registered images.  Higher values mean smoother warps, at the expense of a lower mean squared difference between the images.']};

nits = cfg_entry;
nits.name = 'Iterations for warping';
nits.tag  = 'nits';
nits.strtype = 'n';
nits.num  = [1 1];
nits.def     = @(val)cg_rSPM_get_defaults('hdw.nits_reg', val{:});
nits.help = {'Number of iterations for the warping.'};

warp_opts = cfg_branch;
warp_opts.name = 'Warping Options';
warp_opts.tag = 'warp_opts';
warp_opts.val = {subsamp,weight,nits,warpreg};
warp_opts.help = {'There are a couple of user-customisable warping options.'};

%------------------------------------------------------------------------

hdw = cfg_exbranch;
hdw.name = 'High-Dimensional Warping';
hdw.tag  = 'hdw';
hdw.val  = {esubjs bias_opts warp_opts};
hdw.prog = @cg_hdw;
hdw.vout = @vout_hdw;
hdw.help = {
'This option provides a Bayesian method for three dimensional registration of brain images/* \cite{ashburner00a} */. A finite element approach is used to obtain a maximum a posteriori (MAP) estimate of the deformation field at every voxel of a template volume.  The priors used by the MAP estimate penalize unlikely deformations and enforce a continuous one-to-one mapping.  The deformations are assumed to have some form of symmetry, in that priors describing the probability distribution of the deformations should be identical to those for the inverses (i.e., warping brain A to brain B should not be different probablistically from warping B to A).  A gradient descent algorithm is used to estimate the optimum deformations.'
''
'Deformation fields are written with the same name as the moved image, but with "y_" prefixed on to the filename.  Jacobian determinant images are also written (prefixed by "jy_").'};

%------------------------------------------------------------------------

%------------------------------------------------------------------------

bb = cfg_entry;
bb.name = 'Bounding box';
bb.tag  = 'bb';
bb.num  = [2 3];
bb.strtype = 'e';
bb.def     = @(val)cg_rSPM_get_defaults('write.bb', val{:});
bb.help = {[...
'The bounding box (in mm) of the volume which is to be written (relative to the anterior commissure).']};

%------------------------------------------------------------------------

vox = cfg_entry;
vox.name = 'Voxel sizes';
vox.tag  = 'vox';
vox.num  = [1 3];
vox.strtype = 'e';
vox.def     = @(val)cg_rSPM_get_defaults('write.vox', val{:});
vox.help = {'The voxel sizes (x, y & z, in mm) of the written normalised images.'};

%------------------------------------------------------------------------

interp = cfg_menu;
interp.name = 'Interpolation';
interp.tag  = 'interp';
interp.labels = {'Nearest neighbour','Trilinear','2nd Degree B-spline',...
'3rd Degree B-Spline ','4th Degree B-Spline ','5th Degree B-Spline',...
'6th Degree B-Spline','7th Degree B-Spline'};
interp.values = {0,1,2,3,4,5,6,7};
interp.def     = @(val)cg_rSPM_get_defaults('write.interp', val{:});
interp.help = {...
['The method by which the images are sampled when being written in a different space.'],...
['    Nearest Neighbour:     - Fastest, but not normally recommended.'],...
['    Bilinear Interpolation:     - OK for PET, or realigned fMRI.'],...
['    B-spline Interpolation:     - Better quality (but slower) interpolation/* \cite{thevenaz00a}*/, especially with higher degree splines.  Do not use B-splines when       there is any region of NaN or Inf in the images. '],...
};

%------------------------------------------------------------------------

wrap = cfg_menu;
wrap.name = 'Wrapping';
wrap.tag  = 'wrap';
wrap.labels = {'No wrap','Wrap X','Wrap Y','Wrap X & Y','Wrap Z',...
'Wrap X & Z','Wrap Y & Z','Wrap X, Y & Z'};
wrap.values = {[0 0 0],[1 0 0],[0 1 0],[1 1 0],[0 0 1],[1 0 1],[0 1 1],[1 1 1]};
wrap.def     = @(val)cg_rSPM_get_defaults('write.wrap', val{:});
wrap.help = {...
'These are typically:',...
['    No wrapping: for PET or images that have already                   been spatially transformed. '],...
['    Wrap in  Y: for (un-resliced) MRI where phase encoding                   is in the Y direction (voxel space).']};

%------------------------------------------------------------------------

roptions = cfg_branch;
roptions.name = 'Writing Options';
roptions.tag  = 'roptions';
roptions.val  = {bb,vox,interp,wrap};
roptions.help = {'Various options for writing normalised images.'};

%------------------------------------------------------------------------

matname = cfg_files;
matname.name = 'Parameter File';
matname.tag  = 'matname';
matname.num  = [1 1];
matname.filter  = 'mat';
matname.ufilter = '.*_sn\.mat$';
matname.help = {[...
'Select the ''_sn.mat'' file containing the spatial normalisation parameters for that subject.']};

%------------------------------------------------------------------------

data = cfg_files;
data.name    = 'Files';
data.tag     = 'data';
data.filter  = 'image';
data.num = [1 Inf];
data.help = {'List of subjects.'};

%------------------------------------------------------------------------

data = cfg_files;
data.name    = 'Parameter files';
data.tag     = 'data';
data.filter  = 'mat';
data.ufilter = '.*_sn\.mat$';
data.num = [1 Inf];
data.help = {'List of subjects.'};

writjac = cfg_exbranch;
writjac.name = 'Write Jacobian determinant';
writjac.tag  = 'writejac';
writjac.val  = {data,roptions};
writjac.prog = @writejac;
%writjac.vout = @vout_writejac;
writjac.help   = {[...
'Use previously estimated warps (stored in imagename''_sn.mat'' files) to write Jacobian determinant (volume changes).']};

%------------------------------------------------------------------------

outpt = cfg_entry;
outpt.name = 'Output Filename';
outpt.tag  = 'output';
outpt.strtype = 's';
outpt.num  = [1 Inf];
outpt.val  = {'volumes.csv'};
outpt.help  = {[...
'The volumes will written to current working directory unless a valid full pathname is given.']};

data = cfg_files;
data.name    = 'Parameter files';
data.tag     = 'data';
data.filter  = 'mat';
data.ufilter = '.*_sn\.mat$';
data.num = [1 Inf];
data.help = {'List of subjects.'};

calcpaxinos = cfg_exbranch;
calcpaxinos.name = 'Calculate volume of Paxinos atlas';
calcpaxinos.tag  = 'calcpaxinos';
calcpaxinos.val  = {data,outpt};
calcpaxinos.prog = @volume_paxinos;
calcpaxinos.help = {[...
'Use previously estimated warps (stored in imagename''_sn.mat'' files) to calculate volumes of Paxinos atlas.']};
%------------------------------------------------------------------------

opts = cfg_choice;
opts.name = 'rSPM';
opts.tag  = 'rspm';
opts.values = {writjac,calcpaxinos,hdw};
opts.help = {
'Help text needed!'};
%------------------------------------------------------------------------

return;
%------------------------------------------------------------------------
 
%------------------------------------------------------------------------
function volume_paxinos(varargin)
job    = varargin{1};

fid = fopen(job.output,'w');

spm_progress_bar('Init',length(job.data),'Calculate volumes','files completed');
for i=1:length(job.data),
	[pth, name] = fileparts(strvcat(job.data{i}));
	fprintf('Calculate %s',name(1:end-3));
	[vol_name, vol_left{i}, vol_right{i}] = cg_volume_paxinos(strvcat(job.data{i}));
	% print names for left and right hemisphere
	if i==1
		fprintf(fid,'Name');
		for j=1:size(vol_name,1)
			fprintf(fid,';%s (L);%s (R)',deblank(vol_name(j,:)),deblank(vol_name(j,:)));
		end
		fprintf(fid,'\n');
	end
	fprintf(fid,'%s',name(1:end-3));
	for j=1:size(vol_name,1)
		fprintf(fid,';%7.5f;%7.5f',vol_left{i}(j), vol_right{i}(j));
	end
	fprintf(fid,'\n');
	fprintf('\n');
	spm_progress_bar('Set',i);
end;

spm_progress_bar('Clear');
fclose(fid);

return;
%------------------------------------------------------------------------

%------------------------------------------------------------------------
function writejac(varargin)
job    = varargin{1};
o      = job.roptions;
rflags = struct(...
	'bb',      o.bb,...
	'vox',     o.vox,...
	'interp',  o.interp,...
	'wrap',    o.wrap);

for i=1:length(job.data),
	cg_write_jacdet(strvcat(job.data{i}),rflags);
end;
return;
%------------------------------------------------------------------------

%------------------------------------------------------------------------
function dep = vout_hdw(job)
for k=1:numel(job.subj)
    dep(k)            = cfg_dep;
    dep(k).sname      = sprintf('Jacobian determinant (Subj %d)',k);
    dep(k).src_output = substruct('()',{k},'.','files');
    dep(k).tgt_spec   = cfg_findspec({{'filter','image','strtype','e'}});
end;
%------------------------------------------------------------------------

%------------------------------------------------------------------------
function dep = vout_writejac(job)
for k=1:numel(job.data)
    dep(k)            = cfg_dep;
    dep(k).sname      = sprintf('Norm Params File (Subj %d)',k);
    dep(k).src_output = substruct('()',{k},'.','files');
    dep(k).tgt_spec   = cfg_findspec({{'filter','image','strtype','e'}});
end;
%------------------------------------------------------------------------

%------------------------------------------------------------------------
function vf = vfiles_writejac(varargin)
job = varargin{1};
vf  = {};
for i=1:length(job.data),
    res = job.subj(i).matname;
    vf1 = cell(1,length(res));
    for j=1:length(res),
        [pth,nam,ext,num] = spm_fileparts(res{j});
        vf1{j} = fullfile(pth,['jw', nam(1:end-3), '.img', num]);
    end;
    vf = {vf{:} vf1{:}};
end;
