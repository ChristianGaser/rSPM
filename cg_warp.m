function cg_warp(G,F,Deformations,vx,param)
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
% $Id: cg_warp.m 224 2009-12-02 23:39:15Z gaser $

rev = '$Rev: 224 $';

disp('Compiling cg_warp.c')

pth = fileparts(which(mfilename));
p_path = pwd;
cd(pth);
mex -O cg_warp.c 
cd(p_path);

cg_warp(G,F,Deformations,vx,param)

return