function cg_preprocess_rats
%
% $Id$

spmdir = spm('dir');

for j = 1:1000,
	P{j} = spm_select(Inf,'image',['Select images, subject ' num2str(j)]);
	if isempty(P{j}), break; end;
end;

% use only affine nomalization for longitudinal data (# of files > 1)
if size(P{1},1) > 1
  cutoff = Inf;
  fprintf('Processing of longitudinal data\n');
else
  cutoff = 2;
  fprintf('Processing of cross-sectional data\n');
end

m = size(P,2) - 1;

com = spm_input('Use center-of-mass for initial realignment ?','+1','m','yes|no',[1 0],1);

% Step 1
% Use registration to T2 template to get better starting estimates for the next steps
for j=1:m
	n = size(P{j},1);
	for i=1:n
	  file = deblank(P{j}(i,:));
    matlabbatch{1}.spm.spatial.coreg.estimate.ref = {fullfile(spmdir,'toolbox','rSPM','T2-Paxinos-avg176.nii')};  
    if com
      cg_set_com(file);
    end
    matlabbatch{1}.spm.spatial.coreg.estimate.source = {file};
    matlabbatch{1}.spm.spatial.coreg.estimate.other = {''};
    matlabbatch{1}.spm.spatial.coreg.estimate.eoptions.cost_fun = 'nmi';
    matlabbatch{1}.spm.spatial.coreg.estimate.eoptions.sep = 0.6;
    matlabbatch{1}.spm.spatial.coreg.estimate.eoptions.tol = [0.1 0.1 0.1 0.005 0.005 0.005 0.05 0.05 0.05 0.005 0.005 0.005];
    matlabbatch{1}.spm.spatial.coreg.estimate.eoptions.fwhm = [7 7];
	
	  spm_jobman('run', matlabbatch);
  end
end

clear matlabbatch

% Step 2:
% 1st realignement
% normalization
% 2nd realignement
for j=1:m
	n = size(P{j},1);
	C = cell(n,1);
	for i=1:n
		C{i} = deblank(P{j}(i,:));
	end
	
	% 1st registration
	matlabbatch{1}.spm.spatial.realign.estimate.data = {C};
	matlabbatch{1}.spm.spatial.realign.estimate.eoptions.quality = 0.9;
	matlabbatch{1}.spm.spatial.realign.estimate.eoptions.sep = 0.4;
	matlabbatch{1}.spm.spatial.realign.estimate.eoptions.fwhm = 0.5;
	matlabbatch{1}.spm.spatial.realign.estimate.eoptions.rtm = 0;
	matlabbatch{1}.spm.spatial.realign.estimate.eoptions.interp = 2;
	matlabbatch{1}.spm.spatial.realign.estimate.eoptions.wrap = [0 0 0];
	matlabbatch{1}.spm.spatial.realign.estimate.eoptions.weight = {''};
	
	% normalization
	matlabbatch{2}.spm.spatial.normalise.estwrite.subj.source(1) = cfg_dep;
	matlabbatch{2}.spm.spatial.normalise.estwrite.subj.source(1).tname = 'Source Image';
	matlabbatch{2}.spm.spatial.normalise.estwrite.subj.source(1).tgt_spec = {};
	matlabbatch{2}.spm.spatial.normalise.estwrite.subj.source(1).sname = 'Realign: Estimate: Realigned Images (Sess 1)';
	matlabbatch{2}.spm.spatial.normalise.estwrite.subj.source(1).src_exbranch = substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1});
	matlabbatch{2}.spm.spatial.normalise.estwrite.subj.source(1).src_output = substruct('.','sess', '()',{1}, '.','cfiles');
	matlabbatch{2}.spm.spatial.normalise.estwrite.subj.wtsrc = '';
	matlabbatch{2}.spm.spatial.normalise.estwrite.subj.resample(1) = cfg_dep;
	matlabbatch{2}.spm.spatial.normalise.estwrite.subj.resample(1).tname = 'Images to Write';
	matlabbatch{2}.spm.spatial.normalise.estwrite.subj.resample(1).tgt_spec = {};
	matlabbatch{2}.spm.spatial.normalise.estwrite.subj.resample(1).sname = 'Realign: Estimate: Realigned Images (Sess 1)';
	matlabbatch{2}.spm.spatial.normalise.estwrite.subj.resample(1).src_exbranch = substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1});
	matlabbatch{2}.spm.spatial.normalise.estwrite.subj.resample(1).src_output = substruct('.','sess', '()',{1}, '.','cfiles');
	matlabbatch{2}.spm.spatial.normalise.estwrite.eoptions.template = {fullfile(spmdir,'toolbox','rSPM','T2-Paxinos-avg176.nii')};
	matlabbatch{2}.spm.spatial.normalise.estwrite.eoptions.weight = {fullfile(spmdir,'toolbox','rSPM','Brainmask-Paxinos-avg176.nii')};
	matlabbatch{2}.spm.spatial.normalise.estwrite.eoptions.smosrc = 0.8;
	matlabbatch{2}.spm.spatial.normalise.estwrite.eoptions.smoref = 0.8;
	matlabbatch{2}.spm.spatial.normalise.estwrite.eoptions.regtype = 'none';
	matlabbatch{2}.spm.spatial.normalise.estwrite.eoptions.cutoff = cutoff;
	matlabbatch{2}.spm.spatial.normalise.estwrite.eoptions.nits = 24;
	matlabbatch{2}.spm.spatial.normalise.estwrite.eoptions.reg = 0.5;
	matlabbatch{2}.spm.spatial.normalise.estwrite.roptions.preserve = 0;
	matlabbatch{2}.spm.spatial.normalise.estwrite.roptions.bb = [-8.5 -16 -12; 8.5 8 1];
	matlabbatch{2}.spm.spatial.normalise.estwrite.roptions.vox = [0.2 0.2 0.2];
	matlabbatch{2}.spm.spatial.normalise.estwrite.roptions.interp = 1;
	matlabbatch{2}.spm.spatial.normalise.estwrite.roptions.wrap = [0 0 0];
	matlabbatch{2}.spm.spatial.normalise.estwrite.roptions.prefix = 'w';
	
	% 2nd registration
	matlabbatch{3}.spm.spatial.realign.estwrite.data{1}(1) = cfg_dep;
	matlabbatch{3}.spm.spatial.realign.estwrite.data{1}(1).tname = 'Session';
	matlabbatch{3}.spm.spatial.realign.estwrite.data{1}(1).tgt_spec = {};
	matlabbatch{3}.spm.spatial.realign.estwrite.data{1}(1).sname = 'Normalise: Estimate & Write: Normalised Images (Subj 1)';
	matlabbatch{3}.spm.spatial.realign.estwrite.data{1}(1).src_exbranch = substruct('.','val', '{}',{2}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1});
	matlabbatch{3}.spm.spatial.realign.estwrite.data{1}(1).src_output = substruct('()',{1}, '.','files');
	matlabbatch{3}.spm.spatial.realign.estwrite.eoptions.quality = 1;
	matlabbatch{3}.spm.spatial.realign.estwrite.eoptions.sep = 0.4;
	matlabbatch{3}.spm.spatial.realign.estwrite.eoptions.fwhm = 0.5;
	matlabbatch{3}.spm.spatial.realign.estwrite.eoptions.rtm = 0;
	matlabbatch{3}.spm.spatial.realign.estwrite.eoptions.interp = 2;
	matlabbatch{3}.spm.spatial.realign.estwrite.eoptions.wrap = [0 0 0];
	matlabbatch{3}.spm.spatial.realign.estwrite.eoptions.weight = {''};
	matlabbatch{3}.spm.spatial.realign.estwrite.roptions.which = [2 0];
	matlabbatch{3}.spm.spatial.realign.estwrite.roptions.interp = 4;
	matlabbatch{3}.spm.spatial.realign.estwrite.roptions.wrap = [0 0 0];
	matlabbatch{3}.spm.spatial.realign.estwrite.roptions.mask = 0;
	matlabbatch{3}.spm.spatial.realign.estwrite.roptions.prefix = 'r';

  if isinf(cutoff)
    % warp between time points
	  matlabbatch{4}.spm.tools.rspm.hdw.subj.mov(1) = cfg_dep;
	  matlabbatch{4}.spm.tools.rspm.hdw.subj.mov(1).tname = 'Images to warp';
	  matlabbatch{4}.spm.tools.rspm.hdw.subj.mov(1).tgt_spec{1}(1).name = 'filter';
	  matlabbatch{4}.spm.tools.rspm.hdw.subj.mov(1).tgt_spec{1}(1).value = 'image';
	  matlabbatch{4}.spm.tools.rspm.hdw.subj.mov(1).tgt_spec{1}(2).name = 'strtype';
	  matlabbatch{4}.spm.tools.rspm.hdw.subj.mov(1).tgt_spec{1}(2).value = 'e';
    matlabbatch{4}.spm.tools.rspm.hdw.subj.mov(1).sname = 'Realign: Estimate & Reslice: Resliced Images (Sess 1)';
    matlabbatch{4}.spm.tools.rspm.hdw.subj.mov(1).src_exbranch = substruct('.','val', '{}',{3}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1});
    matlabbatch{4}.spm.tools.rspm.hdw.subj.mov(1).src_output = substruct('.','sess', '()',{1}, '.','rfiles');
	  matlabbatch{4}.spm.tools.rspm.hdw.bias_opts.nits = 8;
	  matlabbatch{4}.spm.tools.rspm.hdw.bias_opts.fwhm = 10;
	  matlabbatch{4}.spm.tools.rspm.hdw.bias_opts.reg = 1e-06;
	  matlabbatch{4}.spm.tools.rspm.hdw.bias_opts.lmreg = 1e-06;
	  matlabbatch{4}.spm.tools.rspm.hdw.warp_opts.subsamp = {''};
	  matlabbatch{4}.spm.tools.rspm.hdw.warp_opts.weight = {fullfile(spmdir,'toolbox','rSPM','Brainmask-Paxinos-avg176.nii')};
	  matlabbatch{4}.spm.tools.rspm.hdw.warp_opts.nits = 30;
	  matlabbatch{4}.spm.tools.rspm.hdw.warp_opts.reg = 20;

  	% final non-linear normalization if previously only affine normalization was applied
	  matlabbatch{5}.spm.spatial.normalise.estwrite.subj.source(1) = cfg_dep;
	  matlabbatch{5}.spm.spatial.normalise.estwrite.subj.source(1).tname = 'Source Image';
    matlabbatch{5}.spm.spatial.normalise.estwrite.subj.source(1).tgt_spec{1}(1).name = 'filter'; 
    matlabbatch{5}.spm.spatial.normalise.estwrite.subj.source(1).tgt_spec{1}(1).value = 'image';
    matlabbatch{5}.spm.spatial.normalise.estwrite.subj.source(1).tgt_spec{1}(2).name = 'strtype';
    matlabbatch{5}.spm.spatial.normalise.estwrite.subj.source(1).tgt_spec{1}(2).value = 'e';
	  matlabbatch{5}.spm.spatial.normalise.estwrite.subj.source(1).sname = 'Realign: Estimate & Reslice: Resliced Images (Sess 1)';
	  matlabbatch{5}.spm.spatial.normalise.estwrite.subj.source(1).src_exbranch = substruct('.','val', '{}',{3}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1});
	  matlabbatch{5}.spm.spatial.normalise.estwrite.subj.source(1).src_output = substruct('.','sess', '()',{1}, '.','rfiles');
	  matlabbatch{5}.spm.spatial.normalise.estwrite.subj.wtsrc = '';
	  matlabbatch{5}.spm.spatial.normalise.estwrite.subj.resample(1) = cfg_dep;
	  matlabbatch{5}.spm.spatial.normalise.estwrite.subj.resample(1).tname = 'Images to Write';
	  matlabbatch{5}.spm.spatial.normalise.estwrite.subj.resample(1).tgt_spec = {};
	  matlabbatch{5}.spm.spatial.normalise.estwrite.subj.resample(1).sname = 'Realign: Estimate & Reslice: Resliced Images (Sess 1)';
	  matlabbatch{5}.spm.spatial.normalise.estwrite.subj.resample(1).src_exbranch = substruct('.','val', '{}',{3}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1});
	  matlabbatch{5}.spm.spatial.normalise.estwrite.subj.resample(1).src_output = substruct('.','sess', '()',{1}, '.','rfiles');
    matlabbatch{5}.spm.spatial.normalise.estwrite.subj.resample(2) = cfg_dep;
    matlabbatch{5}.spm.spatial.normalise.estwrite.subj.resample(2).tname = 'Images to Write';
    matlabbatch{5}.spm.spatial.normalise.estwrite.subj.resample(2).tgt_spec = {};
    matlabbatch{5}.spm.spatial.normalise.estwrite.subj.resample(2).sname = 'High-Dimensional Warping: Jacobian determinant (Subj 1)';
    matlabbatch{5}.spm.spatial.normalise.estwrite.subj.resample(2).src_exbranch = substruct('.','val', '{}',{4}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1});
    matlabbatch{5}.spm.spatial.normalise.estwrite.subj.resample(2).src_output = substruct('()',{1}, '.','files');
	  matlabbatch{5}.spm.spatial.normalise.estwrite.eoptions.template = {fullfile(spmdir,'toolbox','rSPM','T2-Paxinos-avg176.nii')};
	  matlabbatch{5}.spm.spatial.normalise.estwrite.eoptions.weight = {fullfile(spmdir,'toolbox','rSPM','Brainmask-Paxinos-avg176.nii')};
	  matlabbatch{5}.spm.spatial.normalise.estwrite.eoptions.smosrc = 0.8;
	  matlabbatch{5}.spm.spatial.normalise.estwrite.eoptions.smoref = 0.8;
	  matlabbatch{5}.spm.spatial.normalise.estwrite.eoptions.regtype = 'none';
	  matlabbatch{5}.spm.spatial.normalise.estwrite.eoptions.cutoff = 2;
	  matlabbatch{5}.spm.spatial.normalise.estwrite.eoptions.nits = 24;
	  matlabbatch{5}.spm.spatial.normalise.estwrite.eoptions.reg = 0.5;
	  matlabbatch{5}.spm.spatial.normalise.estwrite.roptions.preserve = 0;
	  matlabbatch{5}.spm.spatial.normalise.estwrite.roptions.bb = [-8.5 -16 -12; 8.5 8 1];
	  matlabbatch{5}.spm.spatial.normalise.estwrite.roptions.vox = [0.2 0.2 0.2];
	  matlabbatch{5}.spm.spatial.normalise.estwrite.roptions.interp = 1;
	  matlabbatch{5}.spm.spatial.normalise.estwrite.roptions.wrap = [0 0 0];
	  matlabbatch{5}.spm.spatial.normalise.estwrite.roptions.prefix = 'w';
	end
	spm_jobman('run', matlabbatch);
end

function cg_set_com(P)

V = spm_vol(P);
n = size(P,1);

com_reference = [0 -4.5 -13.5];

for i=1:n
  fprintf('Correct center-of-mass for %s\n',V(i).fname);
  Affine = eye(4);
  vol = spm_read_vols(V(i));
  avg = mean(vol(:));
  avg = mean(vol(find(vol>avg)));
  
	% don't use background values
	[x,y,z] = ind2sub(size(vol),find(vol>avg));
	com = V(i).mat(1:3,:)*[mean(x) mean(y) mean(z) 1]';
	com = com';

	M = spm_get_space(V(i).fname);
	Affine(1:3,4) = (com - com_reference)';
  spm_get_space(V(i).fname,Affine\M);
end
