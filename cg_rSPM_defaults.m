function cg_rspm_defaults
% Sets the defaults for rSPM
% FORMAT cg_rspm_defaults
%_______________________________________________________________________
%
% This file is intended to be customised for the site.
%
% Care must be taken when modifying this file
%_______________________________________________________________________
% @(#)cg_rspm_m	v1.02 Christian Gaser 2008/06/16

global rspm

tbs = spm('tbs');
for i=1:length(tbs)
	if strcmp(lower(tbs(i).name),'rspm')
		rSPMdir = tbs(i).dir; 
	end
end

% Estimation options
%=======================================================================
rspm.estimate.smosrc  = 0.8;
rspm.estimate.smoref  = 0.8;
rspm.estimate.regtype = 'none';
rspm.estimate.weight  = {fullfile(rSPMdir,'Brainmask-Paxinos.nii')};
rspm.estimate.template= {fullfile(rSPMdir,'T2-Paxinos-avg36.nii')};
rspm.estimate.cutoff  = 2;
rspm.estimate.nits    = 24;
rspm.estimate.reg     = 0.5;
rspm.estimate.wtsrc   = 0;

% Writing options
%=======================================================================
rspm.write.preserve   = 0;
rspm.write.bb         = [[-8.5 -16 -12];[8.5 8 1]];
rspm.write.vox        = [0.2 0.2 0.2];
rspm.write.interp     = 1;
rspm.write.wrap       = [0 0 0];

% HDW options
%=======================================================================
rspm.hdw.subsamp      = {''};
rspm.hdw.weight       = {fullfile(rSPMdir,'Brainmask-Paxinos.nii')};
rspm.hdw.nits_bias    = 8;
rspm.hdw.biasfwhm     = 6;
rspm.hdw.biasreg      = 1e-6;
rspm.hdw.lmreg        = 1e-6;
rspm.hdw.warpreg      = 10;
rspm.hdw.nits_reg     = 50;

% DICOM Import defaults
%=======================================================================
dicom.root    = 'flat';
dicom.format  = 'nii';
dicom.icedims = 0;
