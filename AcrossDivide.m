function [head_vals]=AcrossDivide(DEM,FD,DS_OUT,varargin)
    %%BETA%%BETA%%BETA%%BETA%%BETA%%BETA%%BETA%%BETA%%BETA%%BETA%%BETA%%
    %
    % Function is designed to evaluate across divide differences in various proposed metrics for drainage divide stability. Allows for several different
    % strategies for defining a divide of interest and then evaluates across divide differences in four metrics: (1) mean upstream elevation, 
    % (2) mean upstream gradient, (3) mean upstream local relief, and (4) chi. Seleting a divide of interest depends on the definition of the drainage basins that define this
    % divide, see details of the 'outlet_method' parameter for more discussion. Produces a plot of across divide values and a table of values. Title of subplots will also specify
    % which direction the divide is predicted to move based on that metric. See Whipple et al., 2016, Timescales of landscape response to divide migration and drainage capture: 
    % Implications for role of divide mobility in landscape evolution; JGR-ES for further discussion of the use of these metrics.
    %
    % Required Inputs:
    %       DEM - GRIDobj of the digital elevation model of your area loaded into the workspace
    %       FD - FLOWobj of the flow direction of your area loaded into the workspace
    %       DS_OUT - output structure from the 'DivideStability' function
    % Optional Inputs:
    %       outlet_method ['auto_outlets'] - switch to specify how drainage basin outlets are selected, this will control how divides are defined, see descriptions below:
    %           'auto_outlets' - default method. Defines different drainage basins on the basis of stream outlets. This will work well on smaller areas, but for larger 
    %               DEMs, this may produce extremely large basins and thus may not provide fine enough detail for selecting portions of a drainage divide.
    %           'streamorder' - defines drainage basins based on outlets but also junctions of streams above a provided minimum stream order (default is 3, modify with
    %               'minimum_order' parameter).
    %           'pick_new_outlets' - will open a gui interface to pick the outlets of drainage basins that define the divide of interest. To speed plotting, the original stream
    %               network is downsampled based on stream order, the default is to only plot and snap selected pour points to third order or greater drainages. The minimum stream
    %               order used for this downsampling can be modified with the optional 'minimum_order' parameter. This downsampling does not effect the channel heads as the location of
    %               these is still computed from the original full stream newtork provided to the code
    %           'picked_outlets' - short circuits the gui selection of drainage basins that happens with other 'outlet-method' options and allows the user to provide a list of pour points
    %               that define the drainage basins which share the divide of interest. This option requires that the user also supply an input for the optional 'river_mouths' option. This
    %               input can either be (1) a mx3 array with columns x, y, and pairs of numbers defining particular divides of (2) the path to a valid shapefile. Coordinates for either should be the same 
    %		    as the projection used in the DEM file. If a shapefile, the file must be a point file and must only have one value column. For either type of input, the third columns of numbers 
    %               groups basins into being on the same side of a divide, e.g. all the provided pour points with a 1 in the third column are all on one side of a divide of interest, all the pour
    %               points with a 2 in the third column are on the opposite of the divide. Can provide any number of divide pairs, but must always define both sides of a divide so the maximum value 
    %		    in this third column must be even, e.g. pour points with numbers 1,2,3,4,5,6,7,8 would define 4 divide pairs 1-2, 3-4, 5-6, 7-8.
    %       divide_buffer ['moderate'] - switch to control distance between channel heads within selected basins that are considered to define the divide, options are 
    %               'conservative', 'moderate', and 'lax' with increasing acceptable distance
    %       wl_method ['std_dev'] - switch to determine restrictiveness of criteria for defining whether a divide is stable, 'std_dev' uses the standard deviation to determine if the means of
    %               a give metric overlap, alternatively 'std_err' uses the standard error, which is generally smaller than the standard deviation.
    %       minimum_order [3] - minimum stream order for either defining which confluences should be used to define drainage basins ('outlet_method','streamorder') or how much to downsample the
    %           stream network for plotting purposes ('outlet_method','pick_new_outlets').
    %       river_mouths - mx3 array with columns x,y, and a 1 or 2 indicating which side of the divide that particular pour point is on. Required input for 'outlet_method','picked_outlets'
    %       save_shape_files [false] - flag to output a shapefile of the drainage basins used to define the divide and the values at selected channel heads(true) or to not save them (false - default)
    %       out_file_name ['basinsDivide'] - name for shapefiles and plots if optional save functions are used, DO NOT include a file extension in this name.
    %       save_plot [false] - save the divide stability plot, if true will use the name provided for 'drainage_basin_name'.
    %       save_heads [false] - save the channel head values as a mat file, will use the name provided for 'drainage_basin_name' as a prefix.
    %       extract_profiles [false] - option to extract and save out channel profiles for channels downstream of channel heads used to define the divide,
    %           output is necessary as an input for 'DivideProfiles' function.
    %       plot_style [histograms] - method of displaying divide output metrics, either 'points' or 'histograms'
    %       chi_ref_area [1] - reference area used for calculating chi profiles if 'extract_profiles' is true.
    %       theta_ref [0.5] - reference concavity used for calculating chi profiles if 'extract_profiles' is true.
    %       display_map [true] - display map of selected channel heads (optional argument only recognized if 'outlet_method' is 'picked_outlets' to supress display of map figure)
    % Outputs:
    %       head_vals - mx8 array with columns x and y coordinates of relevant channel heads, mean upstream elevation, mean upstream gradient, mean upstream relief, and chi at these channel heads
    %           and a divide identifying number (e.g. 1, 2, 3, 4,...) to indicate which channel heads are grouped together on one or the other side of a divide an individual channel ID number
    %           
    % Examples:
    %       [channel_head_values]=AcrossDivide(DEM,FD,DivStabil_OUT);
    %       [channel_head_values]=AcrossDivide(DEM,FD,DivStabil_OUT,'outlet_method','streamorder','minimum_order',4);
    %       [channel_head_values]=AcrossDivide(DEM,FD,DivStabil_OUT,'outlet_method','picked_outlets','river_mouths',array_of_river_mouths);
    %       [channel_head_values]=AcrossDivide(DEM,FD,DivStabil_OUT,'outlet_method','picked_outlets','river_mouths','name_of_shape_file'); % DO NOT APPEND THE .shp TO THE NAME OF THE FILE
    %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Function Written by Adam M. Forte - Last Revised Spring 2017 %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    p=inputParser;
    p.FunctionName = 'AcrossDivide';
    addRequired(p,'DEM',@(x) isa(x,'GRIDobj'));
    addRequired(p,'FD',@(x) isa(x,'FLOWobj'));
    addRequired(p,'DS_OUT',@(x) isa(x,'struct'));

    addParamValue(p,'outlet_method','auto_outlets',@(x) ischar(validatestring(x,{'auto_outlets','picked_outlets','streamorder','pick_new_outlets'})));
    addParamValue(p,'river_mouths',[],@(x) (isnumeric(x) && size(x,2)==3) | ischar(x));
    addParamValue(p,'minimum_order',3,@(x) isscalar(x) && isnumeric(x));
    addParamValue(p,'divide_buffer','moderate',@(x) ischar(validatestring(x,{'conservative','moderate','lax'})));
    addParamValue(p,'wl_method','std_dev',@(x) ischar(validatestring(x,{'std_dev','std_err'})));
    addParamValue(p,'save_shape_files',false,@(x) islogical(x));
    addParamValue(p,'out_file_name','basinsDivide',@(x) ischar(x));
    addParamValue(p,'plot_style','histograms',@(x) ischar(validatestring(x,{'points','histograms'})));
    addParamValue(p,'save_plot',false,@(x) islogical(x));
    addParamValue(p,'display_map',true,@(x) islogical(x));
    addParamValue(p,'save_heads',false,@(x) islogical(x));
    addParamValue(p,'extract_profiles',false,@(x) islogical(x));
    addParamValue(p,'chi_ref_area',1,@(x) isscalar(x) && isnumeric(x));
    addParamValue(p,'theta_ref',0.5,@(x) isscalar(x) && isnumeric(x));

    parse(p,DEM,FD,DS_OUT,varargin{:});
    DEM=p.Results.DEM;
    FD=p.Results.FD;
    DS_OUT=p.Results.DS_OUT;

    outlet_method=p.Results.outlet_method;
    rm=p.Results.river_mouths;
    mo=p.Results.minimum_order;
    dbd=p.Results.divide_buffer;
    wl_method=p.Results.wl_method;
    sv_db=p.Results.save_shape_files;
    db_nm=p.Results.out_file_name;
    plot_style=p.Results.plot_style;
    save_plot=p.Results.save_plot;
    dm=p.Results.display_map;
    sv_hv=p.Results.save_heads;
    ep=p.Results.extract_profiles;
    a0=p.Results.chi_ref_area;
    mn=p.Results.theta_ref;


    % Perform some checks on the river mouths input (if valid)
    if strcmp(outlet_method,'picked_outlets') & isempty(rm)
        error('Please provide a m x 3 array of river mouths or the valid path to a shapefile');
    elseif strcmp(outlet_method,'picked_outlets') & ischar(rm)
        disp('Reading shapefile')
        try
            rm_name=[rm '.shp'];
            rm_shp=shaperead(rm_name);
            fn=fieldnames(rm_shp);
            rm=horzcat([rm_shp.X]',[rm_shp.Y]',[rm_shp.(fn{4})]');
        catch
            error('Error reading shapefile, make sure the name provided did NOT include .shp, the shapefile is a point file with a single value column, and you have a license for the Mapping Toolbox');
        end
    end

    % Load in necessary results from DivideStability function
	S=DS_OUT.Stream;
    ra=DS_OUT.Ref_Area;
	chXY=DS_OUT.chXY; chX=chXY(:,1); chY=chXY(:,2);
	chIX=DS_OUT.chIX;
    E=DS_OUT.E;
    R=DS_OUT.R;
    G=DS_OUT.G;
    C=DS_OUT.C;

    % Switch between different methods
    switch outlet_method
    %%%%%%%%%%%%%%%%%%
    case 'auto_outlets'
        disp('Finding Ridgelines...')
    	[RL,DB]=FindRidgeLines(FD,S,'auto_outlets');
        disp('Completed')

    	ch_db=DB.Z(chIX);

        f1=figure(1);
        clf
        set(f1,'Units','normalized','Position',[0.1 0.1 0.75 0.75]);
        hold on
        imageschs(DEM,RL);
        hold off

        mg=msgbox('Zoom or pan to area of interest and then press enter');
        uiwait(mg);

        pause()

        hold on
        plot(S,'-w');
        scatter(chX,chY,5,'w','filled');
        hold off

        str='Click within drainage basins defining one side of the divide of interest, presss "Enter/Return", click within drainage basins on other side of the divide and then press "Enter/Return" again.';
        mg=msgbox(str);
        uiwait(mg);

        % Select and paint channel heads on one side of divide
        [x1,y1]=ginput;
        ix1=coord2ind(DEM,x1,y1);
        db1=DB.Z(ix1);
        chIDX1=ismember(ch_db,db1);
        chIX1=chIX(chIDX1);
        chXY1=chXY(chIDX1,:); chX1=chXY1(:,1); chY1=chXY1(:,2);
        aE1=E(chIDX1); aG1=G(chIDX1); aR1=R(chIDX1); aC1=C(chIDX1);

        hold on
        scatter(chX1,chY1,20,'y','filled');
        hold off

        % Select and paint channel heads on the other side of the divide
        [x2,y2]=ginput;
        ix2=coord2ind(DEM,x2,y2);
        db2=DB.Z(ix2);
        chIDX2=ismember(ch_db,db2);
        chIX2=chIX(chIDX2);
        chXY2=chXY(chIDX2,:);  chX2=chXY2(:,1); chY2=chXY2(:,2);
        aE2=E(chIDX2); aG2=G(chIDX2); aR2=R(chIDX2); aC2=C(chIDX2);

        hold on
        scatter(chX2,chY2,20,'g','filled');
        hold off  
    %%%%%%%%%%%%%%%%%%%         
    case 'streamorder'
        disp('Finding Ridgelines...')
        [RL,DB]=FindRidgeLines(FD,S,'streamorder','minimum_order',mo);
        disp('Completed')

        ch_db=DB.Z(chIX);

        f1=figure(1);
        clf
        set(f1,'Units','normalized','Position',[0.1 0.1 0.75 0.75]);
        hold on
        imageschs(DEM,RL);
        hold off

        mg=msgbox('Zoom or pan to area of interest and then press enter');
        uiwait(mg);

        pause()

        hold on
        plot(S,'-w');
        scatter(chX,chY,5,'w','filled');
        hold off

        str='Click within drainage basins defining one side of the divide of interest, presss "Enter/Return", click within drainage basins on other side of the divide and then press "Enter/Return" again.';
        mg=msgbox(str);
        uiwait(mg);

        % Select and paint channel heads on one side of divide
        [x1,y1]=ginput;
        ix1=coord2ind(DEM,x1,y1);
        db1=DB.Z(ix1);
        chIDX1=ismember(ch_db,db1);
        chIX1=chIX(chIDX1);
        chXY1=chXY(chIDX1,:); chX1=chXY1(:,1); chY1=chXY1(:,2);
        aE1=E(chIDX1); aG1=G(chIDX1); aR1=R(chIDX1); aC1=C(chIDX1);

        hold on
        scatter(chX1,chY1,20,'y','filled');
        hold off

        % Select and paint channel heads on the other side of the divide
        [x2,y2]=ginput;
        ix2=coord2ind(DEM,x2,y2);
        db2=DB.Z(ix2);
        chIDX2=ismember(ch_db,db2);
        chIX2=chIX(chIDX2);
        chXY2=chXY(chIDX2,:);  chX2=chXY2(:,1); chY2=chXY2(:,2);
        aE2=E(chIDX2); aG2=G(chIDX2); aR2=R(chIDX2); aC2=C(chIDX2);

        hold on
        scatter(chX2,chY2,20,'g','filled');
        hold off
    %%%%%%%%%%%%%%%%%%%%%%
    case 'pick_new_outlets'
        f1=figure(1);
        clf 
        set(f1,'Units','normalized','Position',[0.1 0.1 0.75 0.75]);
        hold on
        imageschs(DEM,DEM);
        hold off

        mg=msgbox('Zoom or pan to area of interest and then press enter');
        uiwait(mg);

        pause()

        hold on
        crit=['>=' num2str(mo)];
        Sd=modify(S,'streamorder',crit);
        plot(Sd,'-w');
        hold off

        mg=msgbox('Select outlets of basins that define the drainage divide of interest');
        uiwait(mg);

        [xb,yb]=ginput;
        [xb,yb]=snap2stream(Sd,xb,yb);
        rm=[xb yb];

        disp('Finding Ridgelines...')
        [RL,DB]=FindRidgeLines(FD,S,'outlets','river_mouths',rm);
        % DB=GridExpand(DB,DEM,nan);
        disp('Completed')

        f1=figure(1);
        clf
        set(f1,'Units','normalized','Position',[0.1 0.1 0.75 0.75]);
        hold on
        imageschs(DEM,RL);
        hold off

        ch_db=DB.Z(chIX);

        mg=msgbox('Zoom or pan to area of interest and then press enter');
        uiwait(mg);

        pause()

        hold on
        plot(S,'-w');
        scatter(chX,chY,5,'w','filled');
        hold off

        str='Click within drainage basins defining one side of the divide of interest, presss "Enter/Return", click within drainage basins on other side of the divide and then press "Enter/Return" again.';
        mg=msgbox(str);
        uiwait(mg);

        % Select and paint channel heads on one side of divide
        [x1,y1]=ginput;
        ix1=coord2ind(DEM,x1,y1);
        db1=DB.Z(ix1);
        chIDX1=ismember(ch_db,db1);
        chIX1=chIX(chIDX1);
        chXY1=chXY(chIDX1,:); chX1=chXY1(:,1); chY1=chXY1(:,2);
        aE1=E(chIDX1); aG1=G(chIDX1); aR1=R(chIDX1); aC1=C(chIDX1);

        hold on
        scatter(chX1,chY1,20,'y','filled');
        hold off

        % Select and paint channel heads on the other side of the divide
        [x2,y2]=ginput;
        ix2=coord2ind(DEM,x2,y2);
        db2=DB.Z(ix2);
        chIDX2=ismember(ch_db,db2);
        chIX2=chIX(chIDX2);
        chXY2=chXY(chIDX2,:);  chX2=chXY2(:,1); chY2=chXY2(:,2);
        aE2=E(chIDX2); aG2=G(chIDX2); aR2=R(chIDX2); aC2=C(chIDX2);

        hold on
        scatter(chX2,chY2,20,'g','filled');
        hold off
    %%%%%%%%%%%%%%%%%%%%%
    case 'picked_outlets'
        disp('Snapping picked outlets to streams...')
        max_rm=max(rm(:,3));

        if rem(max_rm,2)~=0
            error('Must be numbers defining either side of the divide, so maximum ID number should be even')
        end

        [rmx,rmy]=snap2stream(S,rm(:,1),rm(:,2));
        num_divides=max_rm/2;

        odds=[1:2:(num_divides*2)-1];
        evens=[2:2:(num_divides*2)];

        % Generate empty cells for outputs
        chX1=cell(num_divides,1); chY1=cell(num_divides,1);
        aE1=cell(num_divides,1); aG1=cell(num_divides,1); aR1=cell(num_divides,1); aC1=cell(num_divides,1);
        chX2=cell(num_divides,1); chY2=cell(num_divides,1);
        aE2=cell(num_divides,1); aG2=cell(num_divides,1); aR2=cell(num_divides,1); aC2=cell(num_divides,1);
        S1=cell(num_divides,1); S2=cell(num_divides,1);

        for ii=1:num_divides
            idx1=rm(:,3)==odds(ii);
            idx2=rm(:,3)==evens(ii);
            rm1=[rmx(idx1) rmy(idx1)];
            rm2=[rmx(idx2) rmy(idx2)];

            disp(['Remaking streams for divide ' num2str(ii)])
            rm1ix=coord2ind(DEM,rm1(:,1),rm1(:,2));
            S1{ii}=STREAMobj(FD,'unit','mapunits','minarea',ra,'outlets',rm1ix); 
            rm2ix=coord2ind(DEM,rm2(:,1),rm2(:,2));
            S2{ii}=STREAMobj(FD,'unit','mapunits','minarea',ra,'outlets',rm2ix);

            disp(['Finding channelheads for divide ' num2str(ii)])
            chIX1=streampoi(S1{ii},'channelheads','ix');
            chIDX1=ismember(chIX,chIX1);
            chXY1=chXY(chIDX1,:); chX1{ii}=chXY1(:,1); chY1{ii}=chXY1(:,2);
            aE1{ii}=E(chIDX1); aG1{ii}=G(chIDX1); aR1{ii}=R(chIDX1); aC1{ii}=C(chIDX1);        

            chIX2=streampoi(S2{ii},'channelheads','ix');
            chIDX2=ismember(chIX,chIX2);
            chXY2=chXY(chIDX2,:);  chX2{ii}=chXY2(:,1); chY2{ii}=chXY2(:,2);
            aE2{ii}=E(chIDX2); aG2{ii}=G(chIDX2); aR2{ii}=R(chIDX2); aC2{ii}=C(chIDX2);     
        end   

        % warning('picked outlets method is currently unfinished')

    end


    % Switch for post processing between picked outlets and everything else

    if strcmp(outlet_method,'picked_outlets')

        for jj=1:num_divides
            % Determine Channel Heads of Interest
            num1=numel(chX1{jj});
            num2=numel(chX2{jj});

            switch dbd
            case 'conservative'
                cut_dist=2.0*sqrt(ra);
            case 'moderate'
                cut_dist=2.5*sqrt(ra);
            case 'lax'
                cut_dist=3.0*sqrt(ra);
            end

            chix1=zeros(num1,1);
            for ii=1:num1
                chXoi=chX1{jj}(ii);
                chYoi=chY1{jj}(ii);

                d=EucDist(chXoi,chYoi,chX2{jj},chY2{jj});

                ix=d<cut_dist;

                if isempty(nonzeros(ix))
                    chix1(ii)=0;
                else
                    chix1(ii)=1;
                end
            end

            chix2=zeros(num2,1);
            for ii=1:num2
                chXoi=chX2{jj}(ii);
                chYoi=chY2{jj}(ii);

                d=EucDist(chXoi,chYoi,chX1{jj},chY1{jj});

                ix=d<cut_dist;

                if isempty(nonzeros(ix))
                    chix2(ii)=0;
                else
                    chix2(ii)=1;
                end
            end

            chix1=logical(chix1);
            chix2=logical(chix2);

            if dm
                f1=figure(jj);
                clf
                set(f1,'Units','normalized','Position',[0.1 0.1 0.75 0.75]);
                hold on
                imageschs(DEM,DEM,'colormap','gray');
                scatter(chX1{jj},chY1{jj},20,'y','filled');
                scatter(chX2{jj},chY2{jj},20,'g','filled');
                scatter(chX1{jj}(chix1),chY1{jj}(chix1),40,'k','filled','MarkerEdgeColor','y');
                scatter(chX2{jj}(chix2),chY2{jj}(chix2),40,'k','filled','MarkerEdgeColor','g');
                xlim([min(vertcat(chX1{jj},chX2{jj}))-100 max(vertcat(chX1{jj},chX2{jj}))+100]);
                ylim([min(vertcat(chY1{jj},chY2{jj}))-100 max(vertcat(chY1{jj},chY2{jj}))+100]);
                hold off
            end

            if sv_db
                disp('Generating drainage basins...')
                SO=union(S1{jj},S2{jj},FD);
                DB_OUT=drainagebasins(FD,SO);
                [db_ms,~,~]=GRIDobj2polygon(DB_OUT,'geometry','Polygon','waitbar',true);
                shapewrite(db_ms,[db_nm '_' num2str(jj) '.shp']);
            end

            [do,p1,p2,Ornt]=DivideOrient(chX1{jj}(chix1),chY1{jj}(chix1),chX2{jj}(chix2),chY2{jj}(chix2));

            E1=aE1{jj}(chix1);
            G1=aG1{jj}(chix1);
            R1=aR1{jj}(chix1);
            C1=aC1{jj}(chix1);

            E2=aE2{jj}(chix2);
            G2=aG2{jj}(chix2);
            R2=aR2{jj}(chix2);
            C2=aC2{jj}(chix2);

            ME1=mean(E1); MG1=mean(G1); MR1=mean(R1); MC1=mean(C1);
            ME2=mean(E2); MG2=mean(G2); MR2=mean(R2); MC2=mean(C2);

            switch wl_method
            case 'std_dev'
                stdE1=std(E1); stdG1=std(G1); stdR1=std(R1); stdC1=std(C1);
                stdE2=std(E2); stdG2=std(G2); stdR2=std(R2); stdC2=std(C2);
            case 'std_err'
                stdE1=std(E1)/sqrt(numel(E1)); stdG1=std(G1)/sqrt(numel(G1)); stdR1=std(R1)/sqrt(numel(R1)); stdC1=std(C1)/sqrt(numel(C1));
                stdE2=std(E2)/sqrt(numel(E2)); stdG2=std(G2)/sqrt(numel(G2)); stdR2=std(R2)/sqrt(numel(R2)); stdC2=std(C2)/sqrt(numel(C2));
            end

            [WL]=WinnersLosers(E1,E2,G1,G2,R1,R2,C1,C2,wl_method);

            [out_str]=OutputParser(WL,Ornt);

            dummy1=ones(numel(E1),1);
            dummy2=ones(numel(E2),1).*2;

            f2=figure(jj+num_divides);
            clf
            set(f2,'Units','normalized','Position',[0.1 0.1 0.60 0.60],'renderer','painters','PaperPositionMode','auto');

            switch plot_style
            case 'histograms'
                %%% HISTOGRAMS %%%
                subplot(2,2,3);
                hold on
                xr=linspace(min(vertcat(E1,E2)),max(vertcat(E1,E2)),50);
                [nE1]=hist(E1,xr); [nE2]=hist(E2,xr);
                b1=bar(xr,nE1,'hist');
                b2=bar(xr,nE2,'hist');
                set(b1,'FaceColor','k');
                set(b2,'FaceColor','none','EdgeColor','r');
                % d1=max(nE1)+1; d2=max(nE2)+1;
                dm1=max(nE1); dm2=max(nE2); am=max([dm1 dm2]);
                if ME1>=ME2
                    d1=am+1;
                    d2=d1+1;
                else
                    d2=am+1;
                    d1=d2+1;
                end
                plot([ME1-stdE1;stdE1+ME1],[d1;d1],'-k','LineWidth',2);
                plot([ME2-stdE2;stdE2+ME2],[d2;d2],'-r','LineWidth',2);
                h1(1)=scatter(ME1,d1,40,'k','filled');
                h1(2)=scatter(ME2,d2,40,'r','filled');
                xlabel('Elevation at Reference Area (m)');
                legend(h1,['Mean of Channels ' p1 ' of Divide'],['Mean of Channels ' p2 ' of Divide'],'location','NorthOutside')
                set(gca,'YTickLabel',[]);    
                title(['Elevation: ' out_str{1}])
                hold off

                subplot(2,2,2);
                hold on
                xr=linspace(min(vertcat(R1,R2)),max(vertcat(R1,R2)),50);
                [nR1]=hist(R1,xr); [nR2]=hist(R2,xr);
                b1=bar(xr,nR1,'hist');
                b2=bar(xr,nR2,'hist');
                set(b1,'FaceColor','k');
                set(b2,'FaceColor','none','EdgeColor','r');
                % d1=max(nR1)+1; d2=max(nR2)+1;
                dm1=max(nR1); dm2=max(nR2); am=max([dm1 dm2]);
                if MR1>=MR2
                    d2=am+1;
                    d1=d2+1;
                else
                    d1=am+1;
                    d2=d1+1;
                end                
                plot([MR1-stdR1;stdR1+MR1],[d1;d1],'-k','LineWidth',2);
                plot([MR2-stdR2;stdR2+MR2],[d2;d2],'-r','LineWidth',2);
                h1(1)=scatter(MR1,d1,40,'k','filled');
                h1(2)=scatter(MR2,d2,40,'r','filled');
                xlabel('Mean Upstream Relief above Reference Area (m)');
                legend(h1,['Mean of Channels ' p1 ' of Divide'],['Mean of Channels ' p2 ' of Divide'],'location','NorthOutside')
                set(gca,'YTickLabel',[]);     
                title(['Relief: ' out_str{2}])
                hold off

                subplot(2,2,4);
                hold on
                xr=linspace(min(vertcat(G1,G2)),max(vertcat(G1,G2)),50);
                [nG1]=hist(G1,xr); [nG2]=hist(G2,xr);
                b1=bar(xr,nG1,'hist');
                b2=bar(xr,nG2,'hist');
                set(b1,'FaceColor','k');
                set(b2,'FaceColor','none','EdgeColor','r');
                % d1=max(nG1)+1; d2=max(nG2)+1;
                dm1=max(nG1); dm2=max(nG2); am=max([dm1 dm2]);
                if MG1>=MG2
                    d2=am+1;
                    d1=d2+1;
                else
                    d1=am+1;
                    d2=d1+1;
                end                
                plot([MG1-stdG1;stdG1+MG1],[d1;d1],'-k','LineWidth',2);
                plot([MG2-stdG2;stdG2+MG2],[d2;d2],'-r','LineWidth',2);
                h1(1)=scatter(MG1,d1,40,'k','filled');
                h1(2)=scatter(MG2,d2,40,'r','filled');
                xlabel('Mean Upstream Gradient above Reference Area');
                legend(h1,['Mean of Channels ' p1 ' of Divide'],['Mean of Channels ' p2 ' of Divide'],'location','NorthOutside')
                set(gca,'YTickLabel',[]);    
                title(['Gradient: ' out_str{3}])
                hold off

                subplot(2,2,1);
                hold on
                xr=linspace(min(vertcat(C1,C2)),max(vertcat(C1,C2)),50);
                [nC1]=hist(C1,xr); [nC2]=hist(C2,xr);
                b1=bar(xr,nC1,'hist');
                b2=bar(xr,nC2,'hist');
                set(b1,'FaceColor','k');
                set(b2,'FaceColor','none','EdgeColor','r');
                % d1=max(nC1)+1; d2=max(nC2)+1;
                dm1=max(nC1); dm2=max(nC2); am=max([dm1 dm2]);
                if MC1>=MC2
                    d1=am+1;
                    d2=d1+1;
                else
                    d2=am+1;
                    d1=d2+1;
                end                
                plot([MC1-stdC1;stdC1+MC1],[d1;d1],'-k','LineWidth',2);
                plot([MC2-stdC2;stdC2+MC2],[d2;d2],'-r','LineWidth',2);
                h1(1)=scatter(MC1,d1,40,'k','filled');
                h1(2)=scatter(MC2,d2,40,'r','filled');
                xlabel('Chi at Reference Area');
                legend(h1,['Mean of Channels ' p1 ' of Divide'],['Mean of Channels ' p2 ' of Divide'],'location','NorthOutside')
                set(gca,'YTickLabel',[]);   
                title(['Chi: ' out_str{4}])
                hold off
            case 'points'
            %%% DOTS %%%%
                subplot(2,2,3);
                hold on
                plot([ME1-stdE1;stdE1+ME1],[1;1],'-k','LineWidth',2);
                scatter(E1,dummy1,20,'k');
                plot([ME2-stdE2;stdE2+ME2],[2;2],'-r','LineWidth',2);
                scatter(E2,dummy2,20,'r');
                h1(1)=scatter(ME1,1,40,'k','filled');
                h1(2)=scatter(ME2,2,40,'r','filled');
                xlabel('Mean Upstream Elevation above Reference Area (m)');
                legend(h1,['Mean of Channels ' p1 ' of Divide'],['Mean of Channels ' p2 ' of Divide'],'location','north')
                set(gca,'YTickLabel',[]);     ylim([0 3]);
                title(['Elevation: ' out_str{1}])
                hold off

                subplot(2,2,2);
                hold on
                plot([MR1-stdR1;stdR1+MR1],[1;1],'-k','LineWidth',2);
                scatter(R1,dummy1,20,'k');
                plot([MR2-stdR2;stdR2+MR2],[2;2],'-r','LineWidth',2);
                scatter(R2,dummy2,20,'r');
                h2(1)=scatter(MR1,1,40,'k','filled');
                h2(2)=scatter(MR2,2,40,'r','filled');
                xlabel('Mean Upstream Relief above Reference Area (m)');
                legend(h2,['Mean of Channels ' p1 ' of Divide'],['Mean of Channels ' p2 ' of Divide'],'location','north')
                set(gca,'YTickLabel',[]);     ylim([0 3]);
                title(['Relief: ' out_str{2}])
                hold off

                subplot(2,2,4);
                hold on
                plot([MG1-stdG1;stdG1+MG1],[1;1],'-k','LineWidth',2);    
                scatter(G1,dummy1,20,'k');
                plot([MG2-stdG2;stdG2+MG2],[2;2],'-r','LineWidth',2);   
                scatter(G2,dummy2,20,'r');
                h3(1)=scatter(MG1,1,40,'k','filled');
                h3(2)=scatter(MG2,2,40,'r','filled');
                xlabel('Mean Upstream Gradient above Reference Area');
                legend(h3,['Mean of Channels ' p1 ' of Divide'],['Mean of Channels ' p2 ' of Divide'],'location','north')
                set(gca,'YTickLabel',[]);     ylim([0 3]);
                title(['Gradient: ' out_str{3}])
                hold off

                subplot(2,2,1);
                hold on
                plot([MC1-stdC1;stdC1+MC1],[1;1],'-k','LineWidth',2); 
                scatter(C1,dummy1,20,'k');
                plot([MC2-stdC2;stdC2+MC2],[2;2],'-r','LineWidth',2); 
                scatter(C2,dummy2,20,'r');
                h3(1)=scatter(MC1,1,40,'k','filled');
                h3(2)=scatter(MC2,2,40,'r','filled');
                xlabel('Chi at Reference Area');
                legend(h3,['Mean of Channels ' p1 ' of Divide'],['Mean of Channels ' p2 ' of Divide'],'location','north')
                set(gca,'YTickLabel',[]);     ylim([0 3]);
                title(['Chi: ' out_str{4}])
                hold off
            end

            % Package Output
            X1=chX1{jj}(chix1); Y1=chY1{jj}(chix1);
            X2=chX2{jj}(chix2); Y2=chY2{jj}(chix2);

            div1=ones(size(X1)).*odds(jj);
            div2=ones(size(X2)).*evens(jj);

            val1=[X1 Y1 E1 G1 R1 C1 div1];
            val2=[X2 Y2 E2 G2 R2 C2 div2];

            head_vals{jj,1}=vertcat(val1,val2);

            if save_plot
                print(f2,'-depsc',[db_nm '_' num2str(jj) '.eps']);
            end
        % Loop end
        end
        head_vals=vertcat(head_vals{:})
    else
        % Determine Channel Heads of Interest
        num1=numel(chX1);
        num2=numel(chX2);

        switch dbd
        case 'conservative'
            cut_dist=2.0*sqrt(ra);
        case 'moderate'
            cut_dist=2.5*sqrt(ra);
        case 'lax'
            cut_dist=3.0*sqrt(ra);
        end

        for ii=1:num1
        	chXoi=chX1(ii);
        	chYoi=chY1(ii);

        	d=EucDist(chXoi,chYoi,chX2,chY2);

        	ix=d<cut_dist;

        	if isempty(nonzeros(ix))
        		chix1(ii)=0;
        	else
        		chix1(ii)=1;
        	end
        end

        for ii=1:num2
        	chXoi=chX2(ii);
        	chYoi=chY2(ii);

        	d=EucDist(chXoi,chYoi,chX1,chY1);

        	ix=d<cut_dist;

        	if isempty(nonzeros(ix))
        		chix2(ii)=0;
        	else
        		chix2(ii)=1;
        	end
        end

        chix1=logical(chix1);
        chix2=logical(chix2);

        if strcmp(outlet_method,'picked_outlets')~=1
            hold on
            scatter(chX1(chix1),chY1(chix1),40,'k','filled','MarkerEdgeColor','y');
            scatter(chX2(chix2),chY2(chix2),40,'k','filled','MarkerEdgeColor','g');
            hold off
        end

        if sv_db
            disp('Generating drainage basins...')
            IXX=GRIDobj(DEM);
            IXX.Z(chIX1)=1; IXX.Z(chIX2)=1;
            IXX.Z=logical(IXX.Z);
            SO=modify(S,'downstreamto',IXX);
            DB_OUT=drainagebasins(FD,SO);
            [db_ms,~,~]=GRIDobj2polygon(DB_OUT,'geometry','Polygon','waitbar',true);
            shapewrite(db_ms,[db_nm '.shp']);
        end

        [do,p1,p2,Ornt]=DivideOrient(chX1(chix1),chY1(chix1),chX2(chix2),chY2(chix2));

        E1=aE1(chix1);
        G1=aG1(chix1);
        R1=aR1(chix1);
        C1=aC1(chix1);

        E2=aE2(chix2);
        G2=aG2(chix2);
        R2=aR2(chix2);
        C2=aC2(chix2);

        ME1=mean(E1); MG1=mean(G1); MR1=mean(R1); MC1=mean(C1);
        ME2=mean(E2); MG2=mean(G2); MR2=mean(R2); MC2=mean(C2);

        switch wl_method
        case 'std_dev'
            stdE1=std(E1); stdG1=std(G1); stdR1=std(R1); stdC1=std(C1);
            stdE2=std(E2); stdG2=std(G2); stdR2=std(R2); stdC2=std(C2);
        case 'std_err'
            stdE1=std(E1)/sqrt(numel(E1)); stdG1=std(G1)/sqrt(numel(G1)); stdR1=std(R1)/sqrt(numel(R1)); stdC1=std(C1)/sqrt(numel(C1));
            stdE2=std(E2)/sqrt(numel(E2)); stdG2=std(G2)/sqrt(numel(G2)); stdR2=std(R2)/sqrt(numel(R2)); stdC2=std(C2)/sqrt(numel(C2));
        end

        [WL]=WinnersLosers(E1,E2,G1,G2,R1,R2,C1,C2,wl_method);

        [out_str]=OutputParser(WL,Ornt);

        dummy1=ones(numel(E1),1);
        dummy2=ones(numel(E2),1).*2;

        f2=figure(2);
        clf
        set(f2,'Units','normalized','Position',[0.1 0.1 0.60 0.60],'renderer','painters','PaperPositionMode','auto');

        switch plot_style
        case 'histograms'
            %%% HISTOGRAMS %%%
            subplot(2,2,3);
            hold on
            xr=linspace(min(vertcat(E1,E2)),max(vertcat(E1,E2)),50);
            [nE1]=hist(E1,xr); [nE2]=hist(E2,xr);
            b1=bar(xr,nE1,'hist');
            b2=bar(xr,nE2,'hist');
            set(b1,'FaceColor','k');
            set(b2,'FaceColor','none','EdgeColor','r');
            % d1=max(nE1)+1; d2=max(nE2)+1;
            dm1=max(nE1); dm2=max(nE2); am=max([dm1 dm2]);
            if ME1>=ME2
                d1=am+1;
                d2=d1+1;
            else
                d2=am+1;
                d1=d2+1;
            end            
            plot([ME1-stdE1;stdE1+ME1],[d1;d1],'-k','LineWidth',2);
            plot([ME2-stdE2;stdE2+ME2],[d2;d2],'-r','LineWidth',2);
            h1(1)=scatter(ME1,d1,40,'k','filled');
            h1(2)=scatter(ME2,d2,40,'r','filled');
            xlabel('Elevation at Reference Area (m)');
            legend(h1,['Mean of Channels ' p1 ' of Divide'],['Mean of Channels ' p2 ' of Divide'],'location','NorthOutside')
            set(gca,'YTickLabel',[]);    
            title(['Elevation: ' out_str{1}])
            hold off

            subplot(2,2,2);
            hold on
            xr=linspace(min(vertcat(R1,R2)),max(vertcat(R1,R2)),50);
            [nR1]=hist(R1,xr); [nR2]=hist(R2,xr);
            b1=bar(xr,nR1,'hist');
            b2=bar(xr,nR2,'hist');
            set(b1,'FaceColor','k');
            set(b2,'FaceColor','none','EdgeColor','r');
            % d1=max(nR1)+1; d2=max(nR2)+1;
            dm1=max(nR1); dm2=max(nR2); am=max([dm1 dm2]);
            if MR1>=MR2
                d2=am+1;
                d1=d2+1;
            else
                d1=am+1;
                d2=d1+1;
            end            
            plot([MR1-stdR1;stdR1+MR1],[d1;d1],'-k','LineWidth',2);
            plot([MR2-stdR2;stdR2+MR2],[d2;d2],'-r','LineWidth',2);
            h1(1)=scatter(MR1,d1,40,'k','filled');
            h1(2)=scatter(MR2,d2,40,'r','filled');
            xlabel('Mean Upstream Relief above Reference Area (m)');
            legend(h1,['Mean of Channels ' p1 ' of Divide'],['Mean of Channels ' p2 ' of Divide'],'location','NorthOutside')
            set(gca,'YTickLabel',[]);     
            title(['Relief: ' out_str{2}])
            hold off

            subplot(2,2,4);
            hold on
            xr=linspace(min(vertcat(G1,G2)),max(vertcat(G1,G2)),50);
            [nG1]=hist(G1,xr); [nG2]=hist(G2,xr);
            b1=bar(xr,nG1,'hist');
            b2=bar(xr,nG2,'hist');
            set(b1,'FaceColor','k');
            set(b2,'FaceColor','none','EdgeColor','r');
            % d1=max(nG1)+1; d2=max(nG2)+1;
            dm1=max(nG1); dm2=max(nG2); am=max([dm1 dm2]);
            if MG1>=MG2
                d2=am+1;
                d1=d2+1;
            else
                d1=am+1;
                d2=d1+1;
            end            
            plot([MG1-stdG1;stdG1+MG1],[d1;d1],'-k','LineWidth',2);
            plot([MG2-stdG2;stdG2+MG2],[d2;d2],'-r','LineWidth',2);
            h1(1)=scatter(MG1,d1,40,'k','filled');
            h1(2)=scatter(MG2,d2,40,'r','filled');
            xlabel('Mean Upstream Gradient above Reference Area');
            legend(h1,['Mean of Channels ' p1 ' of Divide'],['Mean of Channels ' p2 ' of Divide'],'location','NorthOutside')
            set(gca,'YTickLabel',[]);    
            title(['Gradient: ' out_str{3}])
            hold off

            subplot(2,2,1);
            hold on
            xr=linspace(min(vertcat(C1,C2)),max(vertcat(C1,C2)),50);
            [nC1]=hist(C1,xr); [nC2]=hist(C2,xr);
            b1=bar(xr,nC1,'hist');
            b2=bar(xr,nC2,'hist');
            set(b1,'FaceColor','k');
            set(b2,'FaceColor','none','EdgeColor','r');
            % d1=max(nC1)+1; d2=max(nC2)+1;
            dm1=max(nC1); dm2=max(nC2); am=max([dm1 dm2]);
            if MC1>=MC2
                d1=am+1;
                d2=d1+1;
            else
                d2=am+1;
                d1=d2+1;
            end            
            plot([MC1-stdC1;stdC1+MC1],[d1;d1],'-k','LineWidth',2);
            plot([MC2-stdC2;stdC2+MC2],[d2;d2],'-r','LineWidth',2);
            h1(1)=scatter(MC1,d1,40,'k','filled');
            h1(2)=scatter(MC2,d2,40,'r','filled');
            xlabel('Chi at Reference Area');
            legend(h1,['Mean of Channels ' p1 ' of Divide'],['Mean of Channels ' p2 ' of Divide'],'location','NorthOutside')
            set(gca,'YTickLabel',[]);   
            title(['Chi: ' out_str{4}])
            hold off
        case 'points'
        %%% DOTS %%%%
            subplot(2,2,3);
            hold on
            plot([ME1-stdE1;stdE1+ME1],[1;1],'-k','LineWidth',2);
            scatter(E1,dummy1,20,'k');
            plot([ME2-stdE2;stdE2+ME2],[2;2],'-r','LineWidth',2);
            scatter(E2,dummy2,20,'r');
            h1(1)=scatter(ME1,1,40,'k','filled');
            h1(2)=scatter(ME2,2,40,'r','filled');
            xlabel('Mean Upstream Elevation above Reference Area (m)');
            legend(h1,['Mean of Channels ' p1 ' of Divide'],['Mean of Channels ' p2 ' of Divide'],'location','north')
            set(gca,'YTickLabel',[]);     ylim([0 3]);
            title(['Elevation: ' out_str{1}])
            hold off

            subplot(2,2,2);
            hold on
            plot([MR1-stdR1;stdR1+MR1],[1;1],'-k','LineWidth',2);
            scatter(R1,dummy1,20,'k');
            plot([MR2-stdR2;stdR2+MR2],[2;2],'-r','LineWidth',2);
            scatter(R2,dummy2,20,'r');
            h2(1)=scatter(MR1,1,40,'k','filled');
            h2(2)=scatter(MR2,2,40,'r','filled');
            xlabel('Mean Upstream Relief above Reference Area (m)');
            legend(h2,['Mean of Channels ' p1 ' of Divide'],['Mean of Channels ' p2 ' of Divide'],'location','north')
            set(gca,'YTickLabel',[]);     ylim([0 3]);
            title(['Relief: ' out_str{2}])
            hold off

            subplot(2,2,4);
            hold on
            plot([MG1-stdG1;stdG1+MG1],[1;1],'-k','LineWidth',2);    
            scatter(G1,dummy1,20,'k');
            plot([MG2-stdG2;stdG2+MG2],[2;2],'-r','LineWidth',2);   
            scatter(G2,dummy2,20,'r');
            h3(1)=scatter(MG1,1,40,'k','filled');
            h3(2)=scatter(MG2,2,40,'r','filled');
            xlabel('Mean Upstream Gradient above Reference Area');
            legend(h3,['Mean of Channels ' p1 ' of Divide'],['Mean of Channels ' p2 ' of Divide'],'location','north')
            set(gca,'YTickLabel',[]);     ylim([0 3]);
            title(['Gradient: ' out_str{3}])
            hold off

            subplot(2,2,1);
            hold on
            plot([MC1-stdC1;stdC1+MC1],[1;1],'-k','LineWidth',2); 
            scatter(C1,dummy1,20,'k');
            plot([MC2-stdC2;stdC2+MC2],[2;2],'-r','LineWidth',2); 
            scatter(C2,dummy2,20,'r');
            h3(1)=scatter(MC1,1,40,'k','filled');
            h3(2)=scatter(MC2,2,40,'r','filled');
            xlabel('Chi at Reference Area');
            legend(h3,['Mean of Channels ' p1 ' of Divide'],['Mean of Channels ' p2 ' of Divide'],'location','north')
            set(gca,'YTickLabel',[]);     ylim([0 3]);
            title(['Chi: ' out_str{4}])
            hold off
        end

        % Package Output
        X1=chX1(chix1); Y1=chY1(chix1);
        X2=chX2(chix2); Y2=chY2(chix2);

        div1=ones(size(X1));
        div2=ones(size(X2)).*2;

        val1=[X1 Y1 E1 G1 R1 C1 div1];
        val2=[X2 Y2 E2 G2 R2 C2 div2];

        head_vals=vertcat(val1,val2);

        if save_plot
            print(f2,'-depsc',[db_nm '.eps']);
        end
    % If Switch End
    end

    % Add channel ID numbers to channel heads
    head_groups=unique(head_vals(:,7));

    for ii=1:numel(head_groups)
        hgOI=head_groups(ii);
        idx=head_vals(:,7)==hgOI;
        nl=[1:1:nnz(idx)]';
        head_vals(idx,8)=nl;
    end

    %% Optional Outputs
    % Generate shape file of head values 
    if sv_db
        head_valsN=double(head_vals);
        X=head_valsN(:,1); Y=head_valsN(:,2); Elev=head_valsN(:,3); Grad=head_valsN(:,4); Rlf=head_valsN(:,5); Chi=head_valsN(:,6); Div=head_valsN(:,7); ChID=head_valsN(:,8);
        % Create Geometry column
        Geometry=repmat({'Point'},numel(X),1);
        T=table(Geometry,X,Y,Elev,Grad,Rlf,Chi,Div,ChID);
        MS=table2struct(T);
        shapewrite(MS,[db_nm '_head_vals.shp']);
    end

    % Save out head values in a mat file
    if sv_hv
        save([db_nm '_head_vals.mat'],'head_vals');
    end

    % Extract channel profiles
    if ep
        % Hydrologically condition dem
        zc=mincosthydrocon(S,DEM,'interp',0.1);
        DEMc=GRIDobj(DEM);
        DEMc.Z(DEMc.Z==0)=NaN;
        DEMc.Z(S.IXgrid)=zc;

        w1=waitbar(0,'Extracting stream segments...');
        StreamSegments=struct;
        for jj=1:numel(head_groups)

            idx=head_vals(:,7)==head_groups(jj);

            x=head_vals(idx,1);
            y=head_vals(idx,2);
            num_heads=numel(x);

            StreamSgmnts=cell(num_heads,1);
            ChiSgmnts=cell(num_heads,1);
            for ii=1:num_heads
                chOI=[x1(ii) y1(ii)];

                % Build logical raster
                ix=coord2ind(DEM,chOI(:,1),chOI(:,2));
                IX=GRIDobj(DEM);
                IX.Z(ix)=1;
                [ixmat,X,Y]=GRIDobj2mat(IX);
                ixmat=logical(ixmat);
                IX=GRIDobj(X,Y,ixmat);

                % Extract stream from channel head to outlet
                Sn=modify(S,'downstreamto',IX);

                C=chiplot(Sn,DEMc,A,'a0',1,'mn',theta_ref,'plot',false);

                StreamSgmnts{ii}=Sn;
                ChiSgmnts{ii}=C;
            end
            StreamSegments(jj,1).Streams=StreamSgmnts;
            StreamSegments(jj,1).Chi=ChiSgmnts;
            StreamSegments(jj,1).Side=jj;


            waitbar(jj/numel(head_groups),w1,[ num2str(jj) ' of ' num2str(numel(head_groups)) ' stream segments extracted']);
        end
        close(w1);
        save([db_nm '_segments.mat'],'StreamSegments');
    end
% Main Function End
end

function [do,p1,p2,Ornt]=DivideOrient(x1,y1,x2,y2);
    x=vertcat(x1,x2);
    y=vertcat(y1,y2);
    f1=fit(x1,y1,'poly1');
    f2=fit(x2,y2,'poly1');

    ms=(f1.p1+f2.p1)/2;

    Ornt=struct;

    ang=rad2deg(atan(1/ms));
    if ang>-22.5 & ang<22.5
        do='N-S';
        Ornt.d=1;
        if mean(x1)<mean(x2)
            p1='West';
            Ornt.p1=270;
            p2='East';
            Ornt.p2=90;
        else
            p1='East';
            Ornt.p1=90;
            p2='West';
            Ornt.p2=270;
        end
    elseif ang>=22.5 & ang<67.5
        do='NE-SW';
        Ornt.d=2;
        if mean(x1)<mean(x2) & mean(y1)>mean(y2)
            p1='Northwest';
            Ornt.p1=315;
            p2='Southeast';
            Ornt.p2=135;
        else
            p1='Southeast';
            Ornt.p1=135;
            p2='Northwest';
            Ornt.p2=315;
        end
    elseif ang>=67.5 & ang<=90
        do='E-W';
        Ornt.d=3;
        if mean(y1)<mean(y2)
            p1='South';
            Ornt.p1=180;
            p2='North';
            Ornt.p2=0;
        else
            p1='North';
            Ornt.p1=0;
            p2='South';
            Ornt.p2=180;
        end
    elseif ang>-67.5 & ang<=-22.5
        do='SE-NW';
        Ornt.d=4;
        if mean(x1)>mean(x2) & mean(y1)>mean(y2)
            p1='Northeast';
            Ornt.p1=45;
            p2='Southwest';
            Ornt.p2=225;
        else
            p1='Southwest';
            Ornt.p1=225;
            p2='Northeast';
            Ornt.p2=45;
        end
    elseif ang>=-90 & ang<=-67.5
        do='E-W';
        Ornt.d=3;
        if mean(y1)<mean(y2)
            p1='South';
            Ornt.p1=180;
            p2='North';
            Ornt.p2=0;
        else
            p1='North';
            Ornt.p1=0;
            p2='South';
            Ornt.p2=180;
        end
    end
end

function [WL]=WinnersLosers(E1,E2,G1,G2,R1,R2,C1,C2,wl_method)

    ME1=mean(E1); MG1=mean(G1); MR1=mean(R1); MC1=mean(C1);
    ME2=mean(E2); MG2=mean(G2); MR2=mean(R2); MC2=mean(C2);

    switch wl_method
    case 'std_dev'
        stdE1=std(E1); stdG1=std(G1); stdR1=std(R1); stdC1=std(C1);
        stdE2=std(E2); stdG2=std(G2); stdR2=std(R2); stdC2=std(C2);
        WL=struct;

        if ME1<ME2 && (ME1)<(ME2-stdE2) && (ME1+stdE1)<(ME2);
            WL.E_AV1=1; % Aggressor
            WL.E_AV2=2; % Victim
        elseif ME1>ME2 && (ME1)>(ME2+stdE2) && (ME1-stdE1)>(ME2);
            WL.E_AV1=2;
            WL.E_AV2=1;
        else
            WL.E_AV1=3; % Stable
            WL.E_AV2=3;
        end
            
        if MG1>MG2 && (MG1)>(MG2+stdG2) && (MG1-stdG1)>(MG2);
            WL.G_AV1=1;
            WL.G_AV2=2;
        elseif MG1<MG2 && (MG1)<(MG2-stdG2) && (MG1+stdG1)<(MG2);
            WL.G_AV1=2;
            WL.G_AV2=1;
        else
            WL.G_AV1=3;
            WL.G_AV2=3;
        end

        if MR1>MR2 && (MR1)>(MR2+stdR2) && (MR1-stdR1)>(MR2);
            WL.R_AV1=1;
            WL.R_AV2=2;
        elseif MR1<MR2 && (MR1)<(MR2-stdR2) && (MR1+stdR1)<(MR2);
            WL.R_AV1=2;
            WL.R_AV2=1;
        else
            WL.R_AV1=3;
            WL.R_AV2=3;
        end

        if MC1<MC2 && (MC1)<(MC2-stdC2) && (MC1+stdC1)<(MC2);
            WL.C_AV1=1; % Aggressor
            WL.C_AV2=2; % Victim
        elseif MC1>MC2 && (MC1)>(MC2+stdC2) && (MC1-stdC1)>(MC2);
            WL.C_AV1=2;
            WL.C_AV2=1;
        else
            WL.C_AV1=3; % Stable
            WL.C_AV2=3;
        end

    case 'std_err'
        stdE1=std(E1)/sqrt(numel(E1)); stdG1=std(G1)/sqrt(numel(G1)); stdR1=std(R1)/sqrt(numel(R1)); stdC1=std(C1)/sqrt(numel(C1));
        stdE2=std(E2)/sqrt(numel(E2)); stdG2=std(G2)/sqrt(numel(G2)); stdR2=std(R2)/sqrt(numel(R2)); stdC2=std(C2)/sqrt(numel(C2));

        WL=struct;

        if ME1<ME2 && (ME1+stdE1)<(ME2-stdE2);
            WL.E_AV1=1; % Aggressor
            WL.E_AV2=2; % Victim
        elseif ME1>ME2 && (ME1-stdE1)>(ME2+stdE2);
            WL.E_AV1=2;
            WL.E_AV2=1;
        else
            WL.E_AV1=3; % Stable
            WL.E_AV2=3;
        end
            
        if MG1>MG2 && (MG1-stdG1)>(MG2+stdG2);
            WL.G_AV1=1;
            WL.G_AV2=2;
        elseif MG1<MG2 && (MG1+stdG1)<(MG2-stdG2);
            WL.G_AV1=2;
            WL.G_AV2=1;
        else
            WL.G_AV1=3;
            WL.G_AV2=3;
        end

        if MR1>MR2 && (MR1-stdR1)>(MR2+stdR2);
            WL.R_AV1=1;
            WL.R_AV2=2;
        elseif MR1<MR2 && (MR1+stdR1)<(MR2-stdR2);
            WL.R_AV1=2;
            WL.R_AV2=1;
        else
            WL.R_AV1=3;
            WL.R_AV2=3;
        end

        if MC1<MC2 && (MC1+stdC1)<(MC2-stdC2);
            WL.C_AV1=1; % Aggressor
            WL.C_AV2=2; % Victim
        elseif MC1>MC2 && (MC1-stdC1)>(MC2+stdC2);
            WL.C_AV1=2;
            WL.C_AV2=1;
        else
            WL.C_AV1=3; % Stable
            WL.C_AV2=3;
        end
    end
end

function [str]=OutputParser(WL,Ornt);

    if WL.E_AV1==1
        if Ornt.p1==0
            E_str='Divide is predicted to move south';
        elseif Ornt.p1==45
            E_str='Divide is predicted to move southwest';
        elseif Ornt.p1==90
            E_str='Divide is predicted to move west';
        elseif Ornt.p1==135
            E_str='Divide is predicted to move northwest';
        elseif Ornt.p1==180
            E_str='Divide is predicted to move north';
        elseif Ornt.p1==225
            E_str='Divide is predicted to move northeast';
        elseif Ornt.p1==270
            E_str='Divide is predicted to move east';
        elseif Ornt.p1==315
            E_str='Divide is predicted to move southeast';
        end
    elseif WL.E_AV2==1
        if Ornt.p2==0
            E_str='Divide is predicted to move south';
        elseif Ornt.p2==45
            E_str='Divide is predicted to move southwest';
        elseif Ornt.p2==90
            E_str='Divide is predicted to move west';
        elseif Ornt.p2==135
            E_str='Divide is predicted to move northwest';
        elseif Ornt.p2==180
            E_str='Divide is predicted to move north';
        elseif Ornt.p2==225
            E_str='Divide is predicted to move northeast';
        elseif Ornt.p2==270
            E_str='Divide is predicted to move east';
        elseif Ornt.p2==315
            E_str='Divide is predicted to move southeast';
        end
    elseif WL.E_AV1==3
        E_str='Divide is stable';
    end

    if WL.G_AV1==1
        if Ornt.p1==0
            G_str='Divide is predicted to move south';
        elseif Ornt.p1==45
            G_str='Divide is predicted to move southwest';
        elseif Ornt.p1==90
            G_str='Divide is predicted to move west';
        elseif Ornt.p1==135
            G_str='Divide is predicted to move northwest';
        elseif Ornt.p1==180
            G_str='Divide is predicted to move north';
        elseif Ornt.p1==225
            G_str='Divide is predicted to move northeast';
        elseif Ornt.p1==270
            G_str='Divide is predicted to move east';
        elseif Ornt.p1==315
            G_str='Divide is predicted to move southeast';
        end
    elseif WL.G_AV2==1
        if Ornt.p2==0
            G_str='Divide is predicted to move south';
        elseif Ornt.p2==45
            G_str='Divide is predicted to move southwest';
        elseif Ornt.p2==90
            G_str='Divide is predicted to move west';
        elseif Ornt.p2==135
            G_str='Divide is predicted to move northwest';
        elseif Ornt.p2==180
            G_str='Divide is predicted to move north';
        elseif Ornt.p2==225
            G_str='Divide is predicted to move northeast';
        elseif Ornt.p2==270
            G_str='Divide is predicted to move east';
        elseif Ornt.p2==315
            G_str='Divide is predicted to move southeast';
        end
    elseif WL.G_AV1==3
        G_str='Divide is stable';
    end

    if WL.R_AV1==1
        if Ornt.p1==0
            R_str='Divide is predicted to move south';
        elseif Ornt.p1==45
            R_str='Divide is predicted to move southwest';
        elseif Ornt.p1==90
            R_str='Divide is predicted to move west';
        elseif Ornt.p1==135
            R_str='Divide is predicted to move northwest';
        elseif Ornt.p1==180
            R_str='Divide is predicted to move north';
        elseif Ornt.p1==225
            R_str='Divide is predicted to move northeast';
        elseif Ornt.p1==270
            R_str='Divide is predicted to move east';
        elseif Ornt.p1==315
            R_str='Divide is predicted to move southeast';
        end
    elseif WL.R_AV2==1
        if Ornt.p2==0
            R_str='Divide is predicted to move south';
        elseif Ornt.p2==45
            R_str='Divide is predicted to move southwest';
        elseif Ornt.p2==90
            R_str='Divide is predicted to move west';
        elseif Ornt.p2==135
            R_str='Divide is predicted to move northwest';
        elseif Ornt.p2==180
            R_str='Divide is predicted to move north';
        elseif Ornt.p2==225
            R_str='Divide is predicted to move northeast';
        elseif Ornt.p2==270
            R_str='Divide is predicted to move east';
        elseif Ornt.p2==315
            R_str='Divide is predicted to move southeast';
        end
    elseif WL.R_AV1==3
        R_str='Divide is stable';
    end

    if WL.C_AV1==1
        if Ornt.p1==0
            C_str='Divide is predicted to move south';
        elseif Ornt.p1==45
            C_str='Divide is predicted to move southwest';
        elseif Ornt.p1==90
            C_str='Divide is predicted to move west';
        elseif Ornt.p1==135
            C_str='Divide is predicted to move northwest';
        elseif Ornt.p1==180
            C_str='Divide is predicted to move north';
        elseif Ornt.p1==225
            C_str='Divide is predicted to move northeast';
        elseif Ornt.p1==270
            C_str='Divide is predicted to move east';
        elseif Ornt.p1==315
            C_str='Divide is predicted to move southeast';
        end
    elseif WL.C_AV2==1
        if Ornt.p2==0
            C_str='Divide is predicted to move south';
        elseif Ornt.p2==45
            C_str='Divide is predicted to move southwest';
        elseif Ornt.p2==90
            C_str='Divide is predicted to move west';
        elseif Ornt.p2==135
            C_str='Divide is predicted to move northwest';
        elseif Ornt.p2==180
            C_str='Divide is predicted to move north';
        elseif Ornt.p2==225
            C_str='Divide is predicted to move northeast';
        elseif Ornt.p2==270
            C_str='Divide is predicted to move east';
        elseif Ornt.p2==315
            C_str='Divide is predicted to move southeast';
        end
    elseif WL.C_AV1==3
        C_str='Divide is stable';
    end

    str={E_str;R_str;G_str;C_str};
end


function [d]=EucDist(x1,y1,x2,y2);
	d=sqrt(((x1-x2).^2)+((y1-y2).^2));
end

function [RidgeLines,DB]=FindRidgeLines(FD,S,method,varargin)
	% Helper function to find ridge line locations, defined as top of
	% a drainage basin (i.e. ridge lines are two grids thick, with true divide between them)

	p=inputParser;
	p.FunctionName = 'FindRidgeLines';
	addRequired(p,'FD',@(x) isa(x,'FLOWobj'));
	addRequired(p,'S',@(x) isa(x,'STREAMobj'));
	addRequired(p,'method',@(x) ischar(validatestring(x,{'auto_outlets','streamorder','outlets'})));

	addParamValue(p,'minimum_order',3,@(x) isscalar(x) && isnumeric(x));
    addParamValue(p,'river_mouths',[1 1],@(x) isnumeric(x) && size(x,2)==2);

	parse(p,FD,S,method,varargin{:});
	FD=p.Results.FD;
	S=p.Results.S;
	method=p.Results.method;
    rm=p.Results.river_mouths;

	so=p.Results.minimum_order;

	switch method
	case 'auto_outlets'
		COix=streampoi(S,'outlets','ix');
        DB=drainagebasins(FD,COix);
	case 'streamorder'
		ms=max(streamorder(S));
		if ms<so;
			so=ms;
			warning('Input stream order is greater than maximum, using maximum')
		end
		crit=['>=' num2str(so)];
		SS=modify(S,'streamorder',crit);
		Oix=streampoi(SS,'outlets','ix');
		Oxy=streampoi(SS,'outlets','xy');
		Cix=streampoi(SS,'bconfluences','ix');
		Cxy=streampoi(SS,'bconfluences','xy');
		COix=vertcat(Oix,Cix);
		DB=drainagebasins(FD,COix);
    case 'outlets'
        COxy=rm;
        DB=drainagebasins(FD,COxy(:,1),COxy(:,2));
	end

	[db,X,Y]=GRIDobj2mat(DB);
	db=double(db);

	rl_l=diff(db,1,2);
	rl_d=diff(db,1,1);

	rl_l(rl_l~=0)=1;
	rl_d(rl_d~=0)=1;

	size_left=size(rl_l);
	size_down=size(rl_d);

	rl_u=vertcat(zeros(1,size_down(:,2)),rl_d);
	rl_r=horzcat(zeros(size_left(:,1),1),rl_l);

	rl_d=vertcat(rl_d,zeros(1,size_down(:,2)));
	rl_l=horzcat(rl_l,zeros(size_left(:,1),1));

	rdgs=rl_d+rl_u+rl_r+rl_l;
	rdgs(rdgs~=0)=1;
	RidgeLines=GRIDobj(X,Y,rdgs);
end
