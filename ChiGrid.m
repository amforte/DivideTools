function [ChiOBJ]=ChiGrid(DEM,FD,varargin)
	% Function to produce a continuous chi grid. Has the similar options for controling outlet elevations and network
	% completeness as the companion code 'SetOutlet' which is designed to modify a stream network.
	%
	% If you use this the result of this code in a publication, please cite Forte, A.M. & Whipple, K.X., In Review, Criteria and Tools for Determining
	% Drainage Divide Stability, submitted to EPSL.
	%
	% Reqiured Inputs:
	% 	DEM - DEM Grid Object (assumes unconditioned DEM)
	% 	FD - FLOW object
	%
	% Optional Inputs:
	%	adjust_outlet_elevation [false] - if true will only calculate chi above a given minimum elevation, that you must provide as an argument to the 
	%									'min_elevation' parameter.
	%	min_elevation [] - parameter to set minimum elevation for outlet elevation if 'method' is set to 'elevation'
	%	complete_network_only [true] - if true (default) the code will only populate portions of the stream network that are complete
	%	chi_ref_area [1] - reference area for calculating chi, setting this value to 1 will ensure that slope of the chi-z relationship is equivalent to 
	%			to ksn, but for this function, this value doesn't matter too much
	%	theta_ref [0.5] - reference concavity for calculating chi
	%	
	% Example:
	%	[CHI]=ChiGrid(DEM,FD);
	%	[CHI]=ChiGrid(DEM,FD,'adjust_outlet_elevation',true,'min_elevation',500);
	%
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% Function Written by Adam M. Forte - Last Revised Winter 2017 %
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	% Parse Inputs
	p = inputParser;         
	p.FunctionName = 'ChiGrid';
	addRequired(p,'DEM',@(x) isa(x,'GRIDobj'));
	addRequired(p,'FD', @(x) isa(x,'FLOWobj'));

	addParamValue(p,'file_name_prefix',[],@(x) ischar(x));
	addParamValue(p,'theta_ref',0.5,@(x) isscalar(x) && isnumeric(x));
	addParamValue(p,'chi_ref_area',1,@(x) isscalar(x) && isnumeric(x));
	addParamValue(p,'adjust_outlet_elevation',false,@(x) isscalar(x) && islogical(x));
	addParamValue(p,'complete_networks_only',true,@(x) isscalar(x) && islogical(x));
	addParamValue(p,'min_elevation',[],@(x) isnumeric(x));

	parse(p,DEM,FD,varargin{:});
	DEM=p.Results.DEM;
	FD=p.Results.FD;

	fnp=p.Results.file_name_prefix;
	mn=p.Results.theta_ref;
	a0=p.Results.chi_ref_area;
	abl=p.Results.adjust_outlet_elevation;
	me=p.Results.min_elevation;
	cno=p.Results.complete_networks_only;

	if abl & isempty(me)
		error('Must provide a minimum eleavtion if Adjust Base Level option is true');
	end

	if cno
		% Find nodes influenced by edge (From Topotoolbox blog)
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
	
	% Extract info from FLOWobj		
	if ~abl
		W=~isnan(DEM);
		I=W.Z(FD.ix);
		ix=double(FD.ix(I));
		ixc=double(FD.ixc(I));
		% Recalculate indices
		IX        = zeros(FD.size,'uint32');
		IX(W.Z)   = 1:nnz(W.Z);
		ix      = double(IX(ix));
		ixc     = double(IX(ixc));
		I          = ixc == 0;
		ix(I)    = [];
		ixc(I)   = [];
	elseif abl
		W=DEM>me;
		I=W.Z(FD.ix);
		ix=double(FD.ix(I));
		ixc=double(FD.ixc(I));
		% Recalculate indices
		IX        = zeros(FD.size,'uint32');
		IX(W.Z)   = 1:nnz(W.Z);
		ix      = double(IX(ix));
		ixc     = double(IX(ixc));
		I          = ixc == 0;
		ix(I)    = [];
		ixc(I)   = [];
	end

	% Generate coordinate list
	IXgrid=find(W.Z);
	[x,y]=ind2coord(DEM,IXgrid);

	% Distance between two nodes throughout grid
	d = nan(size(x));
	dedge = sqrt((x(ixc)-x(ix)).^2 + (y(ixc)-y(ix)).^2);
	d(:) = 0;
	d(ix) = dedge;

	% Cumulative trapezoidal numerical integration of draiange area grid
	DA=flowacc(FD).*(DEM.cellsize^2);
	da=DA.Z(IXgrid);
	c = zeros(size(da));
	da = ((a0) ./da).^mn;
	for r = numel(ix):-1:1;
	    c(ix(r)) = c(ixc(r)) + (da(ixc(r))+(da(ix(r))-da(ixc(r)))/2)*d(ix(r));
	end

	% Generate and populate full chi grid
	ChiOBJ=GRIDobj(DEM);
	ChiOBJ.Z(ChiOBJ.Z==0)=NaN;
	ChiOBJ.Z(IXgrid)=c;


	if abl
		% Remove portions of networks that drain to the edge and are above min elevation
		nrc = numel(x);
		M  = sparse(double(ix),double(ixc),true,nrc,nrc);
		V = (sum(M,2)==0) & (sum(M,1)'~=0);
		oix   = find(V);
		oix   = oix(:);
		oix  = IXgrid(oix);
		[cx,cy]=ind2coord(DEM,oix);
		coxy=[cx cy];

		out_els=DEM.Z(oix);

		[xo,yo]=getoutline(DEM,true);

		sz=size(xo);
		if sz(1)==1 & sz(2)>1
			[oxy]=[xo' yo'];
		elseif sz(2)==1 & sz(1)>1
			[oxy]=[xo yo];
		end

		idx=out_els>me & ismember(coxy,oxy,'rows');
		oix(idx)=[]; % Remove outlets that are above the minimum elevation and along the edge of the DEM

		DM=dependencemap(FD,oix);
		DM.Z=double(DM.Z);
		DM.Z(DM.Z==0)=NaN;
		ChiOBJ=ChiOBJ.*DM;
	end

	if ~isempty(fnp)
		GRIDobj2ascii(ChiOBJ,[fnp '_chigrid.txt']);
	end
end

