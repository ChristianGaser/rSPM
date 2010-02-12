function cg_rspm_defaults
% Sets the defaults for rSPM
% FORMAT cg_rspm_defaults
%_______________________________________________________________________
%
% This file is intended to be customised for the site.
%
% Care must be taken when modifying this file
%_______________________________________________________________________
% @(#)cg_rspm_defaults.m	v1.02 Christian Gaser 2008/06/16

global defaults

tbs = spm('tbs');
for i=1:length(tbs)
	if strcmp(lower(tbs(i).name),'rspm')
		rSPMdir = tbs(i).dir; 
	end
end

% Estimation options
%=======================================================================
defaults.rspm.estimate.smosrc  = 0.8;
defaults.rspm.estimate.smoref  = 0.8;
defaults.rspm.estimate.regtype = 'none';
defaults.rspm.estimate.weight  = {fullfile(rSPMdir,'Brainmask-Paxinos.nii')};
defaults.rspm.estimate.template= {fullfile(rSPMdir,'T2-Paxinos-avg36.nii')};
defaults.rspm.estimate.cutoff  = 2;
defaults.rspm.estimate.nits    = 24;
defaults.rspm.estimate.reg     = 0.5;
defaults.rspm.estimate.wtsrc   = 0;

% Writing options
%=======================================================================
defaults.rspm.write.preserve   = 0;
defaults.rspm.write.bb         = [[-8.5 -16 -12];[8.5 8 1]];
defaults.rspm.write.vox        = [0.2 0.2 0.2];
defaults.rspm.write.interp     = 1;
defaults.rspm.write.wrap       = [0 0 0];

% HDW options
%=======================================================================
defaults.rspm.hdw.subsamp      = {''};
defaults.rspm.hdw.weight       = {fullfile(rSPMdir,'Brainmask-Paxinos.nii')};
defaults.rspm.hdw.nits_bias    = 8;
defaults.rspm.hdw.biasfwhm     = 6;
defaults.rspm.hdw.biasreg      = 1e-6;
defaults.rspm.hdw.lmreg        = 1e-6;
defaults.rspm.hdw.warpreg      = 10;
defaults.rspm.hdw.nits_reg     = 50;

% DICOM Import defaults
%=======================================================================
defaults.dicom.root    = 'flat';
defaults.dicom.format  = 'nii';
defaults.dicom.icedims = 0;
