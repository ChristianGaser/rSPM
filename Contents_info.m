% ______________________________________________________________________
% Copyright (C) 2010 Christian Gaser christian.gaser@uni-jena.de
%        ___  ____  __  __
%  ___  / __)(  _ \(  \/  )  
% (___) \__ \ )___/ )    (   Rat Statistical Parametric Mapping
% (_)   (___/(__)  (_/\/\_)  http://dbm.neuro.uni-jena.de/vbm/
%
% ______________________________________________________________________
% $Id$
%  http://dbm.neuro.uni-jena.de
%
% General files
%   INSTALL.txt           - installation instructions
%   rSPM.man              - notes on rSPM toolbox
%   Howto.txt             - step-by-step manual
%
% rSPM functions
%   spm_rSPM.m            - Toolbox wrapper to call functions
%   tbx_cfg_rspm.m        - configure rSPM
%   cg_volume_paxinos.m   - calculate local volume of labeled Paxinos atlas
%   cg_write_jacdet.m     - write jacaobian determinant
%   cg_calc_jacdet.m      - calculate jacaobian determinant
%   cg_rSPM_defaults.m    - sets the defaults for rSPM
%   cg_preprocess_rats.m  - preprocess rat data
%
% Modified SPM functions
%   spm_orthviews.m         
%   spm_sections.m
%   spm_image.m
%   spm_dicom_convert.m
%
% Images
%   T2-Paxinos-avg176.nii - T2-weighted template based on average of 176 rats
%   Brainmask-Paxinos-avg176.nii - brainmask of Paxinos template
%   Ref0.4mm.nii          - reference image for subsampling to 0.4mm
%   Paxinos_labeled.nii   - labeling according to Paxinos atlas
%   Paxinos_label.txt     - name and ID of labels
