function [ChiOBJ]=ChiGrid(DEM,FD,varargin)

	% Parse Inputs
	p = inputParser;         
	p.FunctionName = 'ChiGrid';
	addRequired(p,'DEM',@(x) isa(x,'GRIDobj'));
	addRequired(p,'FD', @(x) isa(x,'FLOWobj'));

	addParamValue(p,'file_name_prefix','batch',@(x) ischar(x));
	addParamValue(p,'theta_ref',0.45,@(x) isscalar(x) && isnumeric(x));
	addParamValue(p,'a0',1,@(x) isscalar(x) && isnumeric(x));
	addParamValue(p,'complete_network_only',false,@(x) islogical(x) & isscalar(x));
	addParamValue(p,'adjust_base_level',false,@(x) isscalar(x) && islogical(x));
	addParamValue(p,'min_elevation',[],@(x) isnumeric(x));

	parse(p,DEM,FD,varargin{:});
	DEM=p.Results.DEM;
	FD=p.Results.FD;

	file_name_prefix=p.Results.file_name_prefix;
	mn=p.Results.theta_ref;
	a0=p.Results.a0;
	abl=p.Results.adjust_base_level;
	me=p.Results.min_elevation;
	cno=p.Results.complete_network_only;

	if cno & ~abl
		warning('Complete Netowrk Only option is only valid if Adjust Base Level option is also true and a minimum elevation is provided');
	elseif abl & isempty(me)
		error('Must provide a minimum eleavtion if Adjust Base Level option is true');
	end
	
	disp('Generating ChiGRID...')

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

	% Generate and populate chi grid
	ChiOBJ=GRIDobj(DEM);
	ChiOBJ.Z(ChiOBJ.Z==0)=NaN;
	ChiOBJ.Z(IXgrid)=c;

	if cno & abl
		% Complete networks only
		disp('Removing incomplete networks...')
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
		% Control for incosistent output of getoutline
		sz=size(xo);
		if sz(1)==1 & sz(2)>1
			[oxy]=[xo' yo'];
		elseif sz(2)==1 & sz(1)>1
			[oxy]=[xo yo];
		end

		idx=out_els>me & ismember(coxy,oxy,'rows');
		oix(idx)=[];

		DM=dependencemap(FD,oix);
		DM.Z=double(DM.Z);
		DM.Z(DM.Z==0)=NaN;
		ChiOBJ=ChiOBJ.*DM;
	end
end


