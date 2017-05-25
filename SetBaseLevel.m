function [SC]=SetBaseLevel(DEM,FD,A,S,method,varargin)
	% Function to adjust base level of streams within a network via either an elevation or drainage area condition.
	% 
	% Reqiured Inputs:
	% 	DEM - DEM Grid Object (assumes unconditioned DEM)
	% 	FD - FLOW object
	% 	A - GRID object of flow accumulations
	%	S - STREAM object
	%	method - controls how stream network base level is adjusted. Options for control of base level are:
	%		'elevation' - extract streams only above a given elevation (provided by the user using the 'min_elevation' parameter) to ensure that base level
	%			elevation for all streams is uniform. If the provided elevation is too low (i.e. some outlets of the unaltered stream network are above this
	%			elevation) then a warning will be displayed, but the code will still run.
	%		'drain_area' - extract streams only below a given maximum drainage area (provided by the user using the 'max_drainage_area' parameter) to ensure
	%			that the outlets of all extracted streams have the same drainage areas. If the provided maximum drainage area is too large (i.e. some outlets
	%			have drainage areas smaller than this maximum) then a warning will be displayed, but the code will still run.
	%		'max_out_elevation' - uses the maximum elevation of all stream outlets to extract streams only above this elevation
	%		'min_out_drain_area' - uses the minimum drainage area of all stream outlets to extract streams only below this drainage area
	%
	% Optional Inputs:
	%	complete_network_only [true] - if true (default) the code will only populate portions of the stream network that are complete all the way to the desired
	%									end condition (e.g. if method is 'elevation' and 'min_elevation' is 500 meters, streams with outlets along the edge of the
	%									DEM that have elevations above 500 meters will be excluded). If false, streams draining out of the edge of the DEM but whose 
	%									outlets are either above the minimum elevation or the below the maximum drainage area will be included in the output STREAMobj
	%	min_elevation [] - parameter to set minimum elevation for base level if 'method' is set to 'elevation' (if no value is provided, the user
	%						will be prompted with a view of the topography to interactively select a base level)
	%	max_drainage_area [] - parameter to set maximum drainage area for base level if 'method' is set to 'drain_area' (if no value is provided, the user
	%						will be prompted with a view of the topography to interactively select a base level)
	%
	% Example:
	%	[SN]=SetBaseLevel(DEM,FD,A,S,'max_out_elevation');
	%	[SN]=SetBaseLevel(DEM,FD,A,S,'drain_area','max_drainage_area',1e8);
	%
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% Function Written by Adam M. Forte - Last Revised Winter 2017 %
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


	% Parse Inputs
	p = inputParser;         
	p.FunctionName = 'SetBaseLevel';
	addRequired(p,'DEM',@(x) isa(x,'GRIDobj'));
	addRequired(p,'FD', @(x) isa(x,'FLOWobj'));
	addRequired(p,'A', @(x) isa(x,'GRIDobj'));
	addRequired(p,'S',@(x) isa(x,'STREAMobj'));
	addRequired(p,'method',@(x) ischar(validatestring(x,{'elevation','drain_area','max_out_elevation','min_out_drain_area'})));

	addParamValue(p,'complete_network_only',true,@(x) islogical(x) & isscalar(x));
	addParamValue(p,'min_elevation',[],@(x) isnumeric(x));
	addParamValue(p,'max_drainage_area',[],@(x) isnumeric(x));

	parse(p,DEM,FD,A,S,method,varargin{:});
	DEM=p.Results.DEM;
	FD=p.Results.FD;
	A=p.Results.A;
	S=p.Results.S;
	method=p.Results.method;

	cno=p.Results.complete_network_only;
	me=p.Results.min_elevation;
	ma=p.Results.max_drainage_area;

	%% Initiate graphical picker if no values for either min drainage area or min elevation are provided
	if strcmp(method,'drain_area') & isempty(ma)
		DA=A.*(A.cellsize^2);
		lDA=log10(DA);
		f1=figure(1);
		set(f1,'Units','normalized','Position',[0.1 0.1 0.75 0.75]);
		hold on
		imageschs(DEM,lDA);
		title('Zoom to desired view and press enter')
		pause()
		title('Pick point on stream to select base level drainage area');
		[x,y]=ginput(1);
		hold off
		close(f1)
		[xn,yn]=snap2stream(S,x,y);
		da_ix=coord2ind(DEM,xn,yn);
		ma=DA.Z(da_ix);
		disp(['Selected drainage area is : ' num2str(ma) ' m^2']);
	elseif strcmp(method,'elevation') & isempty(me)
		f1=figure(1);
		set(f1,'Units','normalized','Position',[0.1 0.1 0.75 0.75]);
		hold on
		imageschs(DEM,DEM);
		title('Zoom to desired view and press enter')
		pause()
		title('Pick point on DEM to select base level elevation');
		[x,y]=ginput(1);
		hold off
		close(f1)
		el_ix=coord2ind(DEM,x,y);
		me=DEM.Z(el_ix);
		disp(['Selected elevation is : ' num2str(me) ' m']);
	end

	%% Main switch between methods
	switch method
	case 'elevation'
		st_el=getnal(S,DEM);
		idx=st_el>=me;
		IX=S.IXgrid(idx);
		W=GRIDobj(DEM);
		W.Z(IX)=1;
		W.Z=logical(W.Z);
		SC=STREAMobj(FD,W);
		% Check to see if all outlets meet the condition
		coix=streampoi(SC,'outlets','ix');
		coel=DEM.Z(coix);
		max_coel=max(coel);
		if sum(coel>me)~=0 & ~cno
			warning(['One or more stream outlets are above the provided elevation, maximum outlet elevation is ' num2str(max_coel)]);
		elseif sum(coel>me)~=0 & cno
			[xo,yo]=getoutline(DEM,true);
			% Control for incosistent output of getoutline
			sz=size(xo);
			if sz(1)==1 & sz(2)>1
				[oxy]=[xo' yo'];
			elseif sz(2)==1 & sz(1)>1
				[oxy]=[xo yo];
			end
			[coxy]=streampoi(SC,'outlets','xy');
			idx2=coel>me & ismember(coxy,oxy,'rows'); % Find streams with outlets greater than min elevation AND along boundary of DEM
			coix(idx2)=[];
			W=GRIDobj(DEM);
			W.Z(coix)=1;
			W.Z=logical(W.Z);
			SC=modify(SC,'upstreamto',W);			
		end

	case 'drain_area'
		st_da=getnal(S,DA);
		idx=st_da<=ma;
		IX=S.IXgrid(idx);
		W=GRIDobj(DEM);
		W.Z(IX)=1;
		W.Z=logical(W.Z);
		SC=STREAMobj(FD,W);	
		% Check to see if all outlets meet the condition
		coix=streampoi(S,'outlets','ix');
		DA=A.*(A.cellsize^2);
		coda=DA.Z(coix);
		min_coda=min(coda);
		if sum(coda<ma)~=0 & ~cno
			warning(['One or more steam outlets have drainage areas less than provided maximum drainage area, minimum outlet drainage area is ' num2str(min_coda)]);
		elseif sum(coda<ma)~=0 & cno
			[xo,yo]=getoutline(DEM,true);
			% Control for incosistent output of getoutline
			sz=size(xo);
			if sz(1)==1 & sz(2)>1
				[oxy]=[xo' yo'];
			elseif sz(2)==1 & sz(1)>1
				[oxy]=[xo yo];
			end
			[coxy]=streampoi(SC,'outlets','xy');
			idx2=coda<ma & ismember(coxy,oxy,'rows'); % Find streams with outlets less than max draiange area AND along boundary of DEM
			coix(idx2)=[];
			W=GRIDobj(DEM);
			W.Z(coix)=1;
			W.Z=logical(W.Z);
			SC=modify(SC,'upstreamto',W);			
		end
	case 'max_out_elevation'
		coix=streampoi(S,'outlets','ix');
		coel=DEM.Z(coix);
		max_coel=max(coel);
		st_el=getnal(S,DEM);
		idx=st_el>=max_coel;
		IX=S.IXgrid(idx);
		W=GRIDobj(DEM);
		W.Z(IX)=1;
		W.Z=logical(W.Z);
		SC=STREAMobj(FD,W);
	case 'min_out_drain_area'
		coix=streampoi(S,'outlets','ix');
		DA=A.*(A.cellsize^2);
		coda=DA.Z(coix);
		min_coda=min(coda);
		st_da=getnal(S,DA);
		idx=st_da<=min_coda;
		IX=S.IXgrid(idx);
		W=GRIDobj(DEM);
		W.Z(IX)=1;
		W.Z=logical(W.Z);
		SC=STREAMobj(FD,W);	
	end
end