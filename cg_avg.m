function cg_avg
%
% $Id: cg_avg.m 30 2011-09-08 20:00:39Z gaser $

global SWD

SCCSid = '1.0';
[Finter,Fgraph,CmdLine] = spm('FnUIsetup','cg_avg',0);
SPMid = spm('FnBanner',mfilename,SCCSid);
spm_help('!ContextHelp',[mfilename,'.m'])

spm2 = 0;
if strcmp(spm('ver'),'SPM2'), spm2 = 1; end

if ~spm2
	P = spm_select(Inf,'image','Select images');
else	
	P = spm_get(Inf,'w*.img','Select images');
end
if isempty(P), error('no input images specified'), end
n = size(P,1);

% propose outputname
[tmp, name]=spm_str_manip(spm_str_manip(P,'t'),'C');
Q = ['avg_' name.s name.e];
pos = strfind(name.e,',1');
if ~isempty(pos)
  name.e = name.e(1:pos-1);
end
Q = spm_input('Output filename',1,'s',Q);

% make function string like (i1+i2+i3)/3
f='(i1';
for i=2:n, f = [f '+i' num2str(i)]; end
f = [f ')/' num2str(n)];

norm = spm_input('Global normalisation ?','+1','yes|no',[1 0],2);

% mask image
mask = spm_input('Masking ?','+1','yes|no',[1 0],2);
if mask
    if ~spm2
        Pm = spm_select(1,'mask','Select mask image');
    else	
		Pm = spm_get(1,'.img','select mask',[SWD '/apriori/']);
	end
	P = strvcat(P,Pm);
	f = ['i' num2str(n+1) '.*' f];
end

hold = 1;
type = 4;
dmtx = 0;
mask = 0;

spm('FigName','ImCalc: working',Finter,CmdLine);
spm('Pointer','Watch')


%-Map input files
%-----------------------------------------------------------------------
Vi = spm_vol(char(P));
if isempty(Vi), error('no input images specified'), end


if norm
	gm=zeros(size(Vi,1),1);
	disp('Calculating globals...');
	for i=1:size(Vi,1), gm(i) = spm_global(Vi(i)); end
	gm_all = mean(gm);
	for i=1:size(Vi,1)
		Vi(i).pinfo(1:2,:) = gm_all*Vi(i).pinfo(1:2,:)/gm(i);
	end
end

%-Work out filename for output image
%------------------------------------------------------------------
Vo = struct(	'fname',	Q,...
		'dim',		Vi(1).dim(1:3),...
		'mat',		Vi(1).mat,...
		'descrip',	'avg');

if ~spm2
	Vo.dt = [spm_type('int16') spm_platform('bigend')];
else
	Vo.dim = [Vo.dim(1:3) type];
end

%-Call spm_imcalc to handle computations
%------------------------------------------------------------------
args = {{dmtx,mask,hold}};
Vo   = spm_imcalc(Vi,Vo,f,args{:});

%-End
%------------------------------------------------------------------
spm('Pointer');
spm('FigName','ImCalc: done',Finter,CmdLine);
