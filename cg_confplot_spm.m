% [signal_change0, xyz] = cg_confplot_spm(SPM,xSPM,hReg,scale,names, Ic)
%
% SPM, xSPM, hReg	- parameters saved in workspace
% scale    - scale factor for precent signal change of data
% name			- optional names of columns given as {'name1','name2'...}
% Ic			  - number of contrast (usually 1  for effects of interest)		
%
% signal_change0 - scaled beta to obtain percent signal change
% xyz			       - coordinates of local cluster maximum			

try
    [xyz,i] = spm_XYZreg('NearestXYZ',spm_XYZreg('GetCoords',hReg),xSPM.XYZmm);
    spm_XYZreg('SetCoords',xyz,hReg);
catch
    [hReg xSPM SPM] = spm_results_ui('Setup');
    [xyz,i] = spm_XYZreg('NearestXYZ',spm_XYZreg('GetCoords',hReg),xSPM.XYZmm);
    spm_XYZreg('SetCoords',xyz,hReg);
end

CI    = 1.6449;					% = spm_invNcdf(1 - 0.05);

%-Colour specifications
%-----------------------------------------------------------------------
Col   = [0 0 0; .8 .8 .8; 1 0 0];

Cplot = 'Parameter estimates';

%-Specify VOI
%-----------------------------------------------------------------------
xY.def    = spm_input('VOI definition...',1,'b',...
			{'sphere','box','cluster','voxel'},[],3);
Q       = ones(1,size(xSPM.XYZmm,2));

if ~exist('scale','var')
    scale = spm_input('Scaling factor','+1', 'm',['1 (no scaling)|100 (for precent signal change)'],[1 100],2);
end

switch xY.def

	case 'sphere'
	%---------------------------------------------------------------
  xY.spec = spm_input('VOI radius (mm)','!+0','r',0,1,[0,Inf]);
	d     = [xSPM.XYZmm(1,:) - xyz(1);
		 xSPM.XYZmm(2,:) - xyz(2);
		 xSPM.XYZmm(3,:) - xyz(3)];
	Q     = find(sum(d.^2) <= xY.spec^2);
    XYZstr = sprintf(' averaged in sphere (radius %d mm)', xY.spec);

	case 'box'
	%---------------------------------------------------------------
	xY.spec = spm_input('box dimensions [x y z] {mm}',...
			'!+0','r','0 0 0',3);
	Q     = find(all(abs(xSPM.XYZmm - xyz*Q) <= xY.spec(:)*Q/2));
    XYZstr = sprintf(' averaged in box dimensions (%3.2f %3.2f %3.2f)', xY.spec);

	case 'cluster'
	%---------------------------------------------------------------
	[x i] = spm_XYZreg('NearestXYZ',xyz,xSPM.XYZmm);
	A     = spm_clusters(xSPM.XYZ);
	Q     = find(A == A(i));
    XYZstr = sprintf(' averaged in cluster');

	case 'voxel'
	%---------------------------------------------------------------
	d     = [xSPM.XYZmm(1,:) - xyz(1);
		 xSPM.XYZmm(2,:) - xyz(2);
		 xSPM.XYZmm(3,:) - xyz(3)];
	d2 = sum(d.^2);
	Q = find(d2==min(d2));
    XYZstr = sprintf(' in voxel');
end

    XYZ     = xSPM.XYZ(:,Q); 		% coordinates

	%-Parameter estimates:   beta = xX.pKX*xX.K*y;
	%-Residual mean square: ResMS = sum(R.^2)/xX.trRV
	%---------------------------------------------------------------
	
	beta0  = spm_get_data(SPM.Vbeta, XYZ);
	beta   = mean(beta0,2);
	
	try
        y      = spm_get_data(SPM.xY.VY, XYZ);
  catch
       warning('No raw data found');
  end
	ResMS  = spm_get_data(SPM.VResMS,XYZ);
	ResMS  = mean(ResMS,2);
	Bcov   = ResMS*SPM.xX.Bcov;
	Bcov   = Bcov;

	% determine which contrast
	%---------------------------------------------------------------
  if ~exist('Ic','var')
      Ic    = spm_input('Which contrast?','!+1','m',{SPM.xCon.name});
	end
  
	TITLE = {Cplot XYZstr};

	% find contrast and related colu,ms in design matrix
	%--------------------------------------------------------------	
  c = SPM.xCon(Ic).c';
  [ind_x, ind_y] = find(c~=0);
  ind_y = unique(ind_y);
  X = SPM.xX.X;
  X = X(:,ind_y);
  n_effects = size(X,2);
  
  if ~exist('names','var')
      define_names = spm_input('Define names?',1,'yes|use numbers',[1 0],1);
      if define_names
          names = [];
          for i=1:n_effects
              new_name = spm_input(['Name for parameter ' num2str(i)],1,'s');
              names = strvcat(names,new_name);
          end
      else
          names = num2str((1:n_effects)');
      end
  end

	% compute contrast of parameter estimates and 90% C.I.
	%--------------------------------------------------------------	
	signal_change0 = SPM.xCon(Ic).c'*beta0;
	signal_change  = SPM.xCon(Ic).c'*beta;
	CI    = CI*sqrt(diag(SPM.xCon(Ic).c'*Bcov*SPM.xCon(Ic).c));

	% GUI figure
	%--------------------------------------------------------------
	h10 = figure(10);
	clf
	set(h10,'Position',[0 800 150 292],'MenuBar','none','NumberTitle','off');
  hNewButton = uicontrol(h10,...
      'Position',[20 250 110 20],...
      'Callback','cg_confplot_spm',...
      'Interruptible','on',...
      'Style','Pushbutton',...
      'String','Plot',...
      'Backgroundcolor',[1 .5 .5]);
  hClearButton = uicontrol(h10,...
      'position',[20 210 110 20],...
      'Callback','clear names Ic scale',...
      'Interruptible','on',...
      'Style','Pushbutton',...
      'string','Reset variables',...
      'backgroundcolor',[1 .5 .5]);
  hCloseButton = uicontrol(h10,...
      'position',[20 170 110 20],...
      'Callback','close(10,11,12)',...
      'Interruptible','on',...
      'Style','Pushbutton',...
      'string','Close windows',...
      'backgroundcolor',[1 .5 .5]);

	% % signal change plot
	%--------------------------------------------------------------
	h11 = figure(11);
	set(h11,'Position',[150 800 400 270],'NumberTitle','off','Toolbar','none');
	cla
	hold on

	% estimates
	%--------------------------------------------------------------
	h     = bar(signal_change');
	set(h,'FaceColor',Col(2,:));

	% standard error
	%--------------------------------------------------------------
	for j = 1:length(signal_change)
		line([j j],([CI(j) 0 - CI(j)] + signal_change(j)),...
			    'LineWidth',2,'Color',Col(3,:))
	end

	title(TITLE,'FontSize',14,'FontWeight','bold')
	ylabel('parameter estimate','FontSize',12)
	set(gca,'XLim',[0.4 (length(signal_change) + 0.6)],'XTick',1:length(signal_change));
  if exist('names','var')
		if size(names,1) == length(signal_change)
			set(gca,'XTickLabel',names);
		end
	end
	hold off

	% prepare raw values for boxplot
	%--------------------------------------------------------------
  y2 = cell(1,n_effects);

  if scale == 1
      y_label = 'raw signal change';
  else
      y_label = 'percent signal change';
  end
  
  try
      for i=1:n_effects
        y2{i} = scale*y(find(X(:,i)==1),:);
      end
	    h12 = figure(12);
	    set(h12,'Position',[550 800 400 270],'NumberTitle','off','Toolbar','none');
	    cla
	    cg_boxplot(y2);

	    TITLE = {'Boxplot of raw data ' XYZstr};
	    title(TITLE,'FontSize',14,'FontWeight','bold')
	    ylabel(y_label,'FontSize',12)
	    set(gca,'XLim',[0.4 (length(signal_change) + 0.6)],'XTick',1:length(signal_change));
      if exist('names','var')
		    if size(names,1) == length(signal_change)
			    set(gca,'XTickLabel',names);
		    end
	    end
  catch
  
  end
 
return
	% bar chart
	%--------------------------------------------------------------
	h = figure(13);
	set(h,'Position',[500 800 400 270],'NumberTitle','off');
	cla
	for i=1:length(signal_change)
		beta1{i} = beta0(i,:);
		signal_change1{i} = signal_change0(i,:);
	end
	
	cg_boxplot(signal_change1);

	TITLE = {'Boxplot of percent signal change' XYZstr};
	title(TITLE,'FontSize',14,'FontWeight','bold')
	ylabel('percent signal change','FontSize',12)
	set(gca,'XLim',[0.4 (length(signal_change) + 0.6)],'XTick',1:length(signal_change));
  if exist('names','var')
		if size(names,1) == length(signal_change)
			set(gca,'XTickLabel',names);
		end
	end
