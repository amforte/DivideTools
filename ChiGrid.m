function [ChiOBJ]=ChiGrid(DEM,FD,varargin)
	% Function to produce a continuous chi grid. Has the similar options for controling outlet elevations and network
	% completeness as the companion code 'SetOutlet' which is designed to modify a stream network.
	%
    % If you use the result of this code in a publication, please cite Forte, A.M. & Whipple, K.X., 2018, Criteria and Tools for Determining
    % Drainage Divide Stability, Earth and Planetary Science Letters, v.493, p.102-112, DOI:10.1016/j.epsl.2018.04.026
	%
	% Reqiured Inputs:
	% 	DEM - DEM Grid Object (assumes unconditioned DEM)
	% 	FD - FLOW object
	%
	% Optional Inputs:
	%	file_name_prefix [] - name of output ascii file of the chi grid, if no input is provided the 
	%			code will not save an output, but the chi grid (as a GRIDobj) will still be produced
	%			in the workspace.
	%	min_elevation [] - parameter to set minimum elevation for calculating chi
	%	complete_networks_only [true] - if true (default) the code will only populate portions of the stream network that are complete. Generally, this
	%			option should probably be left as true (i.e. chi will not be accurate if drainage area is not accurate), but this can be overly agressive
	%			on certain DEMs and when used in tandem with 'min_elevation', it can be slow to calculate as it requires recalculation of the FLOWobj
	%	chi_ref_area [1] - reference area for calculating chi, setting this value to 1 will ensure that slope of the chi-z relationship is equivalent to 
	%			to ksn, but for this function, this value doesn't matter too much
	%	theta_ref [0.5] - reference concavity for calculating chi
	%	
	% Example:
	%	[CHI]=ChiGrid(DEM,FD);
	%	[CHI]=ChiGrid(DEM,FD,'adjust_outlet_elevation',true,'min_elevation',500);
	%
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% Function Written by Adam M. Forte - Last Revised Spring 2018 %
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	% Parse Inputs
	p = inputParser;         
	p.FunctionName = 'ChiGrid';
	addRequired(p,'DEM',@(x) isa(x,'GRIDobj'));
	addRequired(p,'FD', @(x) isa(x,'FLOWobj'));

	addParamValue(p,'file_name_prefix',[],@(x) ischar(x));
	addParamValue(p,'theta_ref',0.5,@(x) isscalar(x) && isnumeric(x));
	addParamValue(p,'chi_ref_area',1,@(x) isscalar(x) && isnumeric(x));
	addParamValue(p,'complete_networks_only',true,@(x) isscalar(x) && islogical(x));
	addParamValue(p,'min_elevation',[],@(x) isnumeric(x));

	parse(p,DEM,FD,varargin{:});
	DEM=p.Results.DEM;
	FD=p.Results.FD;

	fnp=p.Results.file_name_prefix;
	mn=p.Results.theta_ref;
	a0=p.Results.chi_ref_area;
	me=p.Results.min_elevation;
	cno=p.Results.complete_networks_only;

	if isempty(me)
		abl=false;
	else
		abl=true;
	end

	if cno && ~abl
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
		W=dependencemap(FD,oixi);
		% DEM.Z(~mask.Z)=NaN;
		% % Extract info from FLOWobj		
		% W=~isnan(DEM);
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
	elseif cno && abl
		% Recalculate flow directions after removing portions below min elevation
		IX=DEM>me;
		DEM.Z(IX.Z==false)=NaN;
		FD=FLOWobj(DEM,'preprocess','carve');
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
		% Find offending drainage basins and recalculate indicies
		W=dependencemap(FD,oixi);
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
	elseif ~cno && ~abl
		% Extract info from FLOWobj	
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
	elseif ~cno && abl
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


	if ~isempty(fnp)
		GRIDobj2ascii(ChiOBJ,[fnp '_chigrid.txt']);
	end
end