function [OUT]=DivideStability(DEM,FD,varargin)
	%%BETA%%BETA%%BETA%%BETA%%BETA%%BETA%%BETA%%BETA%%BETA%%BETA%%BETA%%
	%
	% Function is designed to produce various proposed metrics for drainage divide stability. Produces a shapefile of the river network with
	% fields for (1) channel elevation, (2) mean upstream gradient, (3) mean upstream local relief, and (4) chi. It also outputs a structure
	% containing these values at channel heads. Structure is for use in the companion code 'AcrossDivide'. See Whipple et al., 2016, Timescales of
	% landscape response to divide migration and drainage capture: Implications for role of divide mobility in landscape evolution; JGR-ES for further
	% discussion of the use of these metrics
	%
	% Required Inputs:
	% 		DEM - GRIDobj of the digital elevation model of your area loaded into the workspace
	% 		FD - FLOWobj of the flow direction of your area loaded into the workspace
	% Optional Inputs:
	% 		ref_area [1e6] - minimum accumulation area to define streams in meters squared and the reference area for which across divide values are
	%			computed.
	% 		rlf_rad [500] - radius for calculating local relief. For this purpose, the relief radius should be smaller or equal to the mean distance 
	%			from the divide to the 'channel head', if provided radius is larger than assumed distance (square root of reference area) then a warning
	%			will appear, but code will still run.
	% 		verbose [false] - flag to provide output of progress through code, can be useful for larger input to files to make sure the code has not hung.
	% 		shape_name ['div_stabil'] - string that will serve as the name for the output shape file, do not add .shp, this will be done automatically
	%		chi_ref_area [1] - reference area for calculating chi, setting this value to 1 will ensure that slope of the chi-z relationship is equivalent to 
	%			to ksn, but for this function, this value doesn't matter too much
	%		theta_ref [0.5] - reference concavity for calculating chi
	% Outputs:
	%		OUT - Structure with the following fields:
	%			chXY - xy coordinates of channel heads
	%			chIX - linear indices of channel heads
	%			E - upstream elevation at channel heads
	%			G - mean upstream gradient at channel heads
	%			R - mean upstream local relief at channel heads
	%			C - chi value at channel heads
	%			Stream - STREAMobj
	%			Ref_Area - reference area used to generate the STREAMobj
	%		Shapefile has four fields:
	%			chan_elev - normalized stream elevation [varies between 0 and 1]
	%			slope - normalized mean upstream gradient [varies between 0 and 1]
	%			relief - normalized mean upstream local relief [varies between 0 and 1]
	%			chi - values of chi along stream network
	%			
	% Examples:
	%		[AREA_OUT]=DivideStability(DEM,FD); 
	%		[AREA_OUT]=DivideStability(DEM,FD,'verbose',true,'ref_area',1e7,'rlf_rad',500); 
	%
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% Function Written by Adam M. Forte - Last Revised Fall 2016 %
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	p=inputParser;
	p.FunctionName = 'DivideStability';
	addRequired(p,'DEM',@(x) isa(x,'GRIDobj'));
	addRequired(p,'FD',@(x) isa(x,'FLOWobj'));

	addParamValue(p,'ref_area',1e6,@(x) isscalar(x) && isnumeric(x));
	addParamValue(p,'rlf_rad',500,@(x) isscalar(x) && isnumeric(x));
	addParamValue(p,'verbose',false,@(x) isscalar(x));
	addParamValue(p,'shape_name','div_stabil',@(x) ischar(x));
	addParamValue(p,'chi_ref_area',1,@(x) isscalar(x) && isnumeric(x));
	addParamValue(p,'theta_ref',0.5,@(x) isscalar(x) && isnumeric(x));

	parse(p,DEM,FD,varargin{:});
	DEM=p.Results.DEM;
	FD=p.Results.FD;

	ref_area=p.Results.ref_area;
	rlf_rad=p.Results.rlf_rad;
	verbose=p.Results.verbose;
	shape_name=p.Results.shape_name;
	chi_ref_area=p.Results.chi_ref_area;
	theta_ref=p.Results.theta_ref;

	if rlf_rad>sqrt(ref_area);
		warning('Radius for calculating local relief is larger than mean hilllslope length from reference area, consider decreasing relief radius or increasing reference area to avoid significant across divide smearing of local relief values')
	end

	% Make new STREAMobj based on input paramaters
	if verbose
		w1=waitbar(0/13,'Building STREAMobj');
	end
	S=STREAMobj(FD,'unit','mapunits','minarea',ref_area); 

	% Calculate Mean Slopes Upstream
	if verbose
		waitbar(1/12,w1,'Building Upstream Slope Raster');
	end
	G=gradient8(DEM);
	if verbose
		waitbar(2/12,w1,'Building Upstream Slope Raster');
	end	
	UpG=upslopestats(FD,G,'mean'); 

	% Calculate Mean Relief Upstream
	if verbose
		waitbar(3/12,w1,'Building Upstream Relief Raster');
	end	
	R=localtopography(DEM,rlf_rad); 
	if verbose
		waitbar(4/12,w1,'Building Upstream Relief Raster');
	end	
	UpR=upslopestats(FD,R,'mean');

	% Calculate chi
	if verbose
		waitbar(5/12,w1,'Calculating chi');
	end	
	A=flowacc(FD);
	DA=A.*(DEM.cellsize^2);
	c=chitransform(S,DA,'a0',chi_ref_area,'mn',theta_ref);
	CHI=GRIDobj(DEM);
	CHI.Z(S.IXgrid)=c;

	% Normalize GRIDs
	if verbose
		waitbar(6/12,w1,'Normalizing Grids');
	end	
	[e,g,r]=getnal(S,DEM,UpG,UpR);
	if verbose
		waitbar(7/12,w1,'Normalizing Grids');
	end	
	minE=nanmin(e);
	minG=nanmin(g);
	minR=nanmin(r);

	UpEN=DEM-minE;
	UpGN=UpG-minG;
	UpRN=UpR-minR;

	if verbose
		waitbar(8/12,w1,'Normalizing Grids');
	end	
	[e,g,r]=getnal(S,DEM,UpG,UpR);
	if verbose
		waitbar(9/12,w1,'Normalizing Grids');
	end	
	maxE=nanmax(e);
	maxG=nanmax(g);
	maxR=nanmax(r);

	UpEN=UpEN/maxE;
	UpGN=UpGN/maxG;
	UpRN=UpRN/maxR;

	% Saving Shapefile
	if verbose
		waitbar(10/12,w1,'Building Shapefile');
	end	
	MS=STREAMobj2mapstruct(S,'seglength',DEM.cellsize*3,'attributes',...
	{'chan_elev' UpEN @max 'slope' UpGN @max 'relief' UpRN @max 'chi' CHI @max});
	if verbose
		waitbar(11/12,w1,'Saving Shapefile');
	end	
	fileName=[shape_name '.shp'];
	shapewrite(MS,fileName); 
	if verbose
		waitbar(12/12,w1,'Saving Shapefile');
		close(w1);
	end	

	% Process GRIDs and extract values at channel heads for use in AcrossDivide tool
	chXY=streampoi(S,'channelheads','xy'); chX=chXY(:,1); chY=chXY(:,2);
	chIX=streampoi(S,'channelheads','ix');

	OUT=struct;
	OUT.chXY=chXY;
	OUT.chIX=chIX;
	OUT.E=DEM.Z(chIX);
	OUT.G=UpG.Z(chIX);
	OUT.R=UpR.Z(chIX);
	OUT.C=CHI.Z(chIX);
	OUT.Stream=S;
	OUT.Ref_Area=ref_area;

end
