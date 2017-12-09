function [SC]=SetOutlet(DEM,FD,A,S,method,varargin)
	% Function to adjust outlet elevation of streams within a network and ensure network is complete for chi analysis.
	%
	% If you use the result of this code in a publication, please cite Forte, A.M. & Whipple, K.X., In Review, Criteria and Tools for Determining
	% Drainage Divide Stability, submitted to EPSL. And while it's in review, check out the supporting text in preprint form at https://eartharxiv.org/anr29
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
	%		'max_out_elevation' - uses the maximum elevation of all stream outlets to extract streams only above this elevation
	%		'complete_only' - does not check elevations, only removes portions of network that are incomplete watersheds (influenced by edge pixels)
	%
	% Optional Inputs:
	%	complete_network_only [true] - if true (default) the code will only populate portions of the stream network that are complete all the way to the desired
	%									end condition (e.g. if method is 'elevation' and 'min_elevation' is 500 meters, streams with outlets along the edge of the
	%									DEM that have elevations above 500 meters will be excluded along with drainage basins which are influenced by pixels along the
	%									edge of the DEM). If false, streams draining out of the edge of the DEM but whose outlets are above the minimum elevation 
	%									and incomplete drainage networks will be included in the output STREAMobj. Cannot be set to false if method is 'complete_only'
	%	min_elevation [] - parameter to set minimum elevation for base level if 'method' is set to 'elevation' (if no value is provided, the user
	%						will be prompted with a view of the topography to interactively select a base level)
	%
	% Example:
	%	[SN]=SetOutlet(DEM,FD,A,S,'max_out_elevation');
	%	[SN]=SetOutlet(DEM,FD,A,S,'elevation','min_elevation',200);
	%	[SN]=SetOutlet(DEM,FD,A,S,'complete_only')
	%
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% Function Written by Adam M. Forte - Last Revised Winter 2017 %
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


	% Parse Inputs
	p = inputParser;         
	p.FunctionName = 'SetOutlet';
	addRequired(p,'DEM',@(x) isa(x,'GRIDobj'));
	addRequired(p,'FD', @(x) isa(x,'FLOWobj'));
	addRequired(p,'A', @(x) isa(x,'GRIDobj'));
	addRequired(p,'S',@(x) isa(x,'STREAMobj'));
	addRequired(p,'method',@(x) ischar(validatestring(x,{'elevation','max_out_elevation','complete_only'})));

	addParamValue(p,'complete_network_only',true,@(x) islogical(x) & isscalar(x));
	addParamValue(p,'min_elevation',[],@(x) isnumeric(x));


	parse(p,DEM,FD,A,S,method,varargin{:});
	DEM=p.Results.DEM;
	FD=p.Results.FD;
	A=p.Results.A;
	S=p.Results.S;
	method=p.Results.method;

	cno=p.Results.complete_network_only;
	me=p.Results.min_elevation;

	if ~cno & strcmp(method,'complete_network_only')
		error('Cannot set method to complete_only and set complete_network_only to false');
	end

	if cno % Remove portions of the DEM that are not complete networks
		% Find nodes influenced by edge (from TopoToolbox blog)
		IXE = GRIDobj(DEM,'logical');
		IXE.Z(:,:) = true;
		IXE.Z(2:end-1,2:end-1) = false;
		IXE = influencemap(FD,IXE);
		% Rest is mine
		% Find drainage basins and all outlets
		[DB,oixi]=drainagebasins(FD);
		% Find where these share pixels other than the edge
		db=DB.Z; db=db(3:end-2,3:end-2); % Move 2 pixels in to be conservative
		ixe=IXE.Z; ixe=ixe(3:end-2,3:end-2); % Move 2 pixels in to be conservative
		dbL=db(ixe);
		% Compile list of drainage basins that are influenced by edge pixels
		dbL=unique(dbL);
		% Index list of outlets based on this
		idxi=ismember(DB.Z(oixi),dbL);
		oixi(idxi)=[];
		% Remove drainage basins based on this
		mask=dependencemap(FD,oixi);
		DEM.Z(~mask.Z)=NaN;
	end


	%% Initiate graphical picker if no values for either min drainage area or min elevation are provided
	if strcmp(method,'elevation') & isempty(me)
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
	case 'complete_only'
		IX=S.IXgrid;
		W=GRIDobj(DEM);
		W.Z(IX)=1;
		W.Z=logical(W.Z);
		IDX=isnan(DEM);
		W.Z(IDX.Z)=false;
		SC=STREAMobj(FD,W);
	end
