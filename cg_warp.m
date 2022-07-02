function cg_warp(G,F,Defx,Defy,Defz,vx,param)
% FORMAT out = ornlmMex(in, v, f, h)
% 
% G            - reference image
% F            - image to warp
% Deformations - deformation field
% vx           - array of voxel size reference/image
% param        - array of [nit,reg,epsilon,meth]
% nit          - number of warping iterations
% reg          - warping regularization (lambda)
% epsilon      - ???
% meth         - ???
%
%
% Christian Gaser
% $Id: cg_warp.m 27 2011-06-21 15:00:15Z gaser $

rev = '$Rev: 27 $';

disp('Compiling cg_warp.c')

pth = fileparts(which(mfilename));
p_path = pwd;
cd(pth);
mex -O cg_warp.c 
cd(p_path);

cg_warp(G,F,Defx,Defy,Defz,vx,param)

return
