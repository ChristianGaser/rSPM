function spm_rSPM
% rSPM Toolbox wrapper to call rSPM functions
%_______________________________________________________________________
% $Id: spm_rSPM.m 41 2016-05-20 10:29:36Z gaser $

addpath(fileparts(which(mfilename)));

SPMid = spm('FnBanner',mfilename,'v1.02');
[Finter,Fgraph,CmdLine] = spm('FnUIsetup','Rat SPM Toolbox');
spm_help('!ContextHelp',mfilename);
spm_help('!Disp','rSPM.man','',Fgraph,'Morphometry toolbox for rats');

fig = spm_figure('GetWin','Interactive');
h0  = uimenu(fig,...
	'Label',	'rSPM',...
	'Separator',	'on',...
	'Tag',		'rSPM',...
	'HandleVisibility','on');
h1  = uimenu(h0,...
	'Label',	'DICOM Import',...
	'Separator',	'off',...
	'Tag',		'Convert DICOM data to nifti',...
	'CallBack','spm_jobman(''interactive'','''',''matlabbatch.spm.util.dicom'');',...
	'HandleVisibility','on');
h11  = uimenu(h0,...
	'Label',	'Correct bias',...
	'Separator',	'off',...
	'Tag',		'Correct bias',...
	'CallBack','cg_correct_bias_rats;',...
	'HandleVisibility','on');
h2  = uimenu(h0,...
	'Label',	'Preprocess data',...
	'Separator',	'off',...
	'Tag',		'Preprocess data',...
	'CallBack','cg_preprocess_rats;',...
	'HandleVisibility','on');
h3  = uimenu(h0,...
	'Label',	'Cross-sectional data',...
	'Separator',	'off',...
	'Tag',		'Cross-sectional data',...
	'HandleVisibility','on');
h31  = uimenu(h3,...
	'Label',	'Write Jacobian determinant',...
	'Separator',	'off',...
	'Tag',		'Write segmentations',...
	'CallBack','spm_jobman(''interactive'','''',''matlabbatch.spm.tools.rspm.writejac'');',...
	'HandleVisibility','on');
h32  = uimenu(h3,...
	'Label',	'Calculate volumes of Paxinos atlas',...
	'Separator',	'off',...
	'Tag',		'Calculate volumes',...
	'CallBack','spm_jobman(''interactive'','''',''matlabbatch.spm.tools.rspm.calcpaxinos'');',...
	'HandleVisibility','on');
h4  = uimenu(h0,...
	'Label',	'Longitudinal data',...
	'Separator',	'off',...
	'Tag',		'Longitudinal data',...
	'HandleVisibility','on');
h41  = uimenu(h4,...
	'Label',	'Calculate volumes of Paxinos atlas',...
	'Separator',	'off',...
	'Tag',		'Calculate volumes',...
	'CallBack','spm_jobman(''interactive'','''',''matlabbatch.spm.tools.rspm.calcpaxinos_long'');',...
	'HandleVisibility','on');
h5  = uimenu(h0,...
	'Label',	'Show results (confidence plots)',...
	'Separator',	'off',...
	'Tag',		'Show results (confidence plots)',...
	'CallBack','cg_confplot_spm;',...
	'HandleVisibility','on');
h6  = uimenu(h0,...
	'Label',	'Check for updates',...
	'Separator',	'off',...
	'Tag',		'Check for updates',...
	'CallBack','cg_rSPM_update(1);',...
	'HandleVisibility','on');
	