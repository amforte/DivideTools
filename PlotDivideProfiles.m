function PlotDivideProfiles(head_vals,StreamSegments,varargin)
    %%BETA%%BETA%%BETA%%BETA%%BETA%%BETA%%BETA%%BETA%%BETA%%BETA%%BETA%%
    %
    % Function plots chi-elevation and distance-elevation profiles for the streams used to evaluate a divide section with the 'AcrossDivide' code.
    %
    % Required Inputs:
    %       head_vals - m x 8 matrix output from 'AcrossDivide' or saved in the '*_head_vals.mat' file if the optional 'save_heads' parameter was
    %			set to true in 'AcrossDivide'
    %		StreamSegments - structure saved out from 'AcrossDivide' within the '*_segments.mat' file. Will only be produced if 'extract_profiles' is 
    %			set to true when running 'AcrossDivide'
    %
    % Optional Inputs:
 	%		label [false] - option to either label the stream segments by their unique channel ID (true) or not (false, default)
 	%		color_by ['side'] - option to control how streams are colored in the plots. Recognized inputs are:
 	%			'side' - all segments on the same side of a divide will be one color and the segments on the opposite side another color
 	%			'gradient' - segments colored by the mean upstream gradient value at their heads
 	%			'relief' - segmetns colored by the mean upstream local relief value at their heads
 	%		channel_list [] - this parameter allows the user to restrict the plotting to specific channels based on their channel ID. Input is expected
 	%			to be a cell array and must have the same number of cells as the number of groups of channels (i.e. if you defined two divides, this means
 	%			you have four groups of stream segments so the cell array you provide must have four cells). If you leave a cell blank in this input, the
 	%			plot will include all channels from this side of the divide.
 	%		name_prefix [] - string to include as a filename prefix for output eps files of plots. If left empty (default), no files will be saved.
    %           
    % Examples:
    %		PlotDivideProfiles(head_vals,StreamSegments);
    %		PlotDivideProfiles(head_vals,StreamSegments,'color_by','gradient');
    %		PlodDivideProfiles(head_vals,StreamSegments,'name_prefix','exmpl');
    %
    %		ch_lst=cell(4,1);
    %		ch_lst{1}=[4,5]; ch_lst{2}=[3,4];
    %		PlotDivideProfiles(head_vals,StreamSegments,'channel_list',ch_lst); %If you had to two divides and wanted to plot only channels 4 & 5 on one side 
    %			of Divide 1 with channels 3 & 4 on the other side of Divide 1, but all the channels that define Divide 2. 
    %			
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Function Written by Adam M. Forte - Last Revised Spring 2017 %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	% Parse Inputs
	p = inputParser;
	p.FunctionName = 'PlotDivideProfiles';
	addRequired(p,'head_vals',@(x) isnumeric(x));
	addRequired(p,'StreamSegments',@(x) isstruct(x));

	addParamValue(p,'label',false,@(x) islogical(x));
	addParamValue(p,'color_by','side',@(x) ischar(validatestring(x,{'side','gradient','relief'})));
	addParamValue(p,'channel_list',[],@(x) iscell(x));
	addParamValue(p,'name_prefix',[], @(x) ischar(x));

	parse(p,head_vals,StreamSegments,varargin{:});
	hv=p.Results.head_vals;
	S=p.Results.StreamSegments;

	lab=p.Results.label;
	color_by=p.Results.color_by;
	chL=p.Results.channel_list;
	name_prefix=p.Results.name_prefix;

	divide_nums=head_vals(:,7);
	chID=head_vals(:,8);
	chX=head_vals(:,1);
	chY=head_vals(:,2);

	chG=head_vals(:,4);
	chR=head_vals(:,5);

	num_divides=round(max(divide_nums)/2);
    odds=[1:2:(num_divides*2)-1];
    evens=[2:2:(num_divides*2)];
	cMap=jet(20);

	% Error check
	if ~isempty(chL) & numel(chL)~=numel(unique(divide_nums));
		error('List of channel IDs provided with "channel_list" must have the same number of cells as there are divide numbers');
	end

    fig_num=1;

    if isempty(chL)
		for ii=1:num_divides
			C1=StreamSegments(odds(ii),1).Chi;
			C2=StreamSegments(evens(ii),1).Chi;		

			idx1=divide_nums==odds(ii);
			idx2=divide_nums==evens(ii);

			x1=chX(idx1); x2=chX(idx2);
			y1=chX(idx1); y2=chX(idx2);
			id1=chID(idx1); id2=chID(idx2);
			g1=chG(idx1); g2=chG(idx2);
			r1=chR(idx1); r2=chR(idx2);

			rL=vertcat(r1,r2); rcL=linspace(min(rL),max(rL),20);
			gL=vertcat(g1,g2); gcL=linspace(min(gL),max(gL),20);


			num_seg1=numel(x1);
			num_seg2=numel(x2);

			[pos1,pos2,ab1,ab2]=DivideOrient(x1,y1,x2,y2);

			f1=figure(fig_num);
			set(f1,'Units','normalized','Position',[0.1 0.1 0.75 0.75],'renderer','painters','PaperSize',[10 10],'PaperPositionMode','auto');
			clf
			colormap(jet);

			for jj=1:num_seg1
				C=C1{jj};

				subplot(2,1,1);
				hold on
				switch color_by
				case 'side'
					p1(1)=plot(C.chi,C.elev,'Color','k');
					if lab
						[mc,ix]=max(C.chi);
						me=C.elev(ix);
						text(mc+mc*0.01,me,[ab1 ' ' num2str(id1(jj))],'Color','k');
					end
				case 'relief'
					rVal=r1(jj);
					[~,vix]=min(abs(rVal-rcL));
					p1(1)=plot(C.chi,C.elev,'Color',cMap(vix,:));
					if lab
						[mc,ix]=max(C.chi);
						me=C.elev(ix);
						text(mc+mc*0.01,me,[ab1 ' ' num2str(id1(jj))],'Color',cMap(vix,:));
					end
				case 'gradient'
					gVal=g1(jj);
					[~,vix]=min(abs(gVal-gcL));				
					p1(1)=plot(C.chi,C.elev,'Color',cMap(vix,:));
					if lab
						[mc,ix]=max(C.chi);
						me=C.elev(ix);
						text(mc+mc*0.01,me,[ab1 ' ' num2str(id1(jj))],'Color',cMap(vix,:));
					end
				end	
				hold off

				subplot(2,1,2);
				hold on
				switch color_by
				case 'side'
					p2(1)=plot(C.distance,C.elev,'Color','k');
					if lab
						[md,ix]=max(C.distance);
						me=C.elev(ix);
						text(md+md*0.01,me,[ab1 ' ' num2str(id1(jj))],'Color','k');
					end
				case 'relief'
					rVal=r1(jj);
					[~,vix]=min(abs(rVal-rcL));
					p2(1)=plot(C.distance,C.elev,'Color',cMap(vix,:));
					if lab
						[mc,ix]=max(C.chi);
						me=C.elev(ix);
						text(mc+mc*0.01,me,[ab1 ' ' num2str(id1(jj))],'Color',cMap(vix,:));
					end
				case 'gradient'
					gVal=g1(jj);
					[~,vix]=min(abs(gVal-gcL));
					p2(1)=plot(C.distance,C.elev,'Color',cMap(vix,:));	
					if lab
						[mc,ix]=max(C.chi);
						me=C.elev(ix);
						text(mc+mc*0.01,me,[ab1 ' ' num2str(id1(jj))],'Color',cMap(vix,:));
					end
				end
				hold off
			end

			for jj=1:num_seg2
				C=C2{jj};

				subplot(2,1,1);
				hold on
				switch color_by
				case 'side'
					p1(2)=plot(C.chi,C.elev,'Color','r');
					if lab
						[mc,ix]=max(C.chi);
						me=C.elev(ix);
						text(mc+mc*0.01,me,[ab1 ' ' num2str(id2(jj))],'Color','r');
					end
				case 'relief'
					rVal=r2(jj);
					[~,vix]=min(abs(rVal-rcL));
					p1(2)=plot(C.chi,C.elev,'-.','Color',cMap(vix,:));
					if lab
						[mc,ix]=max(C.chi);
						me=C.elev(ix);
						text(mc+mc*0.01,me,[ab1 ' ' num2str(id2(jj))],'Color',cMap(vix,:));
					end
				case 'gradient'
					gVal=g2(jj);
					[~,vix]=min(abs(gVal-gcL));				
					p1(2)=plot(C.chi,C.elev,'-.','Color',cMap(vix,:));
					if lab
						[mc,ix]=max(C.chi);
						me=C.elev(ix);
						text(mc+mc*0.01,me,[ab1 ' ' num2str(id2(jj))],'Color',cMap(vix,:));
					end
				end	
				hold off

				subplot(2,1,2);
				hold on
				switch color_by
				case 'side'
					p2(2)=plot(C.distance,C.elev,'Color','r');
					if lab
						[md,ix]=max(C.distance);
						me=C.elev(ix);
						text(md+md*0.01,me,[ab1 ' ' num2str(id2(jj))],'Color','r');
					end
				case 'relief'
					rVal=r2(jj);
					[~,vix]=min(abs(rVal-rcL));
					p2(2)=plot(C.distance,C.elev,'-.','Color',cMap(vix,:));
					if lab
						[mc,ix]=max(C.chi);
						me=C.elev(ix);
						text(mc+mc*0.01,me,[ab1 ' ' num2str(id2(jj))],'Color',cMap(vix,:));
					end
				case 'gradient'
					gVal=g2(jj);
					[~,vix]=min(abs(gVal-gcL));
					p2(2)=plot(C.distance,C.elev,'-.','Color',cMap(vix,:));	
					if lab
						[mc,ix]=max(C.chi);
						me=C.elev(ix);
						text(mc+mc*0.01,me,[ab1 ' ' num2str(id2(jj))],'Color',cMap(vix,:));
					end
				end
				hold off
			end

			LegendEnt={['Streams ' pos1 ' of Divide'],['Streams ' pos2 ' of Divide']};

			subplot(2,1,1);
			hold on
				xlabel('Chi')
				ylabel('Elevation (m)')
				title('Chi - Z')
				legend(p1,LegendEnt,'location','best');
				switch color_by
				case 'relief'
					caxis([min(rL) max(rL)]);
					c1=colorbar;
					ylabel(c1,'Relief');
				case 'gradient'
					caxis([min(gL) max(gL)]);
					c1=colorbar;
					ylabel(c1,'Gradient');
				end
			hold off

			subplot(2,1,2);
			hold on
				xlabel('Distance from Mouth (m)')
				ylabel('Elevation (m)')
				title('Long Profile')
				legend(p2,LegendEnt,'location','best');
				switch color_by
				case 'relief'
					caxis([min(rL) max(rL)]);
					c2=colorbar;
					ylabel(c2,'Relief');
				case 'gradient'
					caxis([min(gL) max(gL)]);
					c2=colorbar;
					ylabel(c2,'Gradient');
				end
			hold off

			drawnow

			if ~isempty(name_prefix)
				print(f1,[name_prefix '_Divide_' num2str(ii) '.eps'],'-depsc');
			end

			fig_num=fig_num+1;
		end
	else
		for ii=1:num_divides
			C1=StreamSegments(odds(ii),1).Chi;
			C2=StreamSegments(evens(ii),1).Chi;		

			idx1=divide_nums==odds(ii);
			idx2=divide_nums==evens(ii);

			x1=chX(idx1); x2=chX(idx2);
			y1=chX(idx1); y2=chX(idx2);
			id1=chID(idx1); id2=chID(idx2);
			g1=chG(idx1); g2=chG(idx2);
			r1=chR(idx1); r2=chR(idx2);

			[pos1,pos2,ab1,ab2]=DivideOrient(x1,y1,x2,y2);			

			chL1=chL{odds(ii)};
			chL2=chL{evens(ii)};

			if isempty(chL1)
				C1=C1;
				x1=x1; y1=y1; id1=id1; g1=g1; r1=r1;
			else
				C1=C1(chL1);
				x1=x1(chL1); y1=y1(chL1); id1=id1(chL1); g1=g1(chL1); r1=r1(chL1);
			end

			if isempty(chL2)
				C2=C2;
				x2=x2; y2=y2; id2=id2; g2=g2; r2=r2;
			else
				C2=C2(chL2);
				x2=x2(chL2); y2=y2(chL2); id2=id2(chL2); g2=g2(chL2); r2=r2(chL2);
			end

			rL=vertcat(r1,r2); rcL=linspace(min(rL),max(rL),20);
			gL=vertcat(g1,g2); gcL=linspace(min(gL),max(gL),20);

			num_seg1=numel(x1);
			num_seg2=numel(x2);


			f1=figure(fig_num);
			set(f1,'Units','normalized','Position',[0.1 0.1 0.75 0.75],'renderer','painters','PaperSize',[10 10],'PaperPositionMode','auto');
			clf
			colormap(jet);

			for jj=1:num_seg1
				C=C1{jj};

				subplot(2,1,1);
				hold on
				switch color_by
				case 'side'
					p1(1)=plot(C.chi,C.elev,'Color','k');
					if lab
						[mc,ix]=max(C.chi);
						me=C.elev(ix);
						text(mc+mc*0.01,me,[ab1 ' ' num2str(id1(jj))],'Color','k');
					end
				case 'relief'
					rVal=r1(jj);
					[~,vix]=min(abs(rVal-rcL));
					p1(1)=plot(C.chi,C.elev,'Color',cMap(vix,:));
					if lab
						[mc,ix]=max(C.chi);
						me=C.elev(ix);
						text(mc+mc*0.01,me,[ab1 ' ' num2str(id1(jj))],'Color',cMap(vix,:));
					end
				case 'gradient'
					gVal=g1(jj);
					[~,vix]=min(abs(gVal-gcL));				
					p1(1)=plot(C.chi,C.elev,'Color',cMap(vix,:));
					if lab
						[mc,ix]=max(C.chi);
						me=C.elev(ix);
						text(mc+mc*0.01,me,[ab1 ' ' num2str(id1(jj))],'Color',cMap(vix,:));
					end
				end	
				hold off

				subplot(2,1,2);
				hold on
				switch color_by
				case 'side'
					p2(1)=plot(C.distance,C.elev,'Color','k');
					if lab
						[md,ix]=max(C.distance);
						me=C.elev(ix);
						text(md+md*0.01,me,[ab1 ' ' num2str(id1(jj))],'Color','k');
					end
				case 'relief'
					rVal=r1(jj);
					[~,vix]=min(abs(rVal-rcL));
					p2(1)=plot(C.distance,C.elev,'Color',cMap(vix,:));
					if lab
						[mc,ix]=max(C.chi);
						me=C.elev(ix);
						text(mc+mc*0.01,me,[ab1 ' ' num2str(id1(jj))],'Color',cMap(vix,:));
					end
				case 'gradient'
					gVal=g1(jj);
					[~,vix]=min(abs(gVal-gcL));
					p2(1)=plot(C.distance,C.elev,'Color',cMap(vix,:));	
					if lab
						[mc,ix]=max(C.chi);
						me=C.elev(ix);
						text(mc+mc*0.01,me,[ab1 ' ' num2str(id1(jj))],'Color',cMap(vix,:));
					end
				end
				hold off
			end

			for jj=1:num_seg2
				C=C2{jj};

				subplot(2,1,1);
				hold on
				switch color_by
				case 'side'
					p1(2)=plot(C.chi,C.elev,'Color','r');
					if lab
						[mc,ix]=max(C.chi);
						me=C.elev(ix);
						text(mc+mc*0.01,me,[ab1 ' ' num2str(id2(jj))],'Color','r');
					end
				case 'relief'
					rVal=r2(jj);
					[~,vix]=min(abs(rVal-rcL));
					p1(2)=plot(C.chi,C.elev,'-.','Color',cMap(vix,:));
					if lab
						[mc,ix]=max(C.chi);
						me=C.elev(ix);
						text(mc+mc*0.01,me,[ab1 ' ' num2str(id2(jj))],'Color',cMap(vix,:));
					end
				case 'gradient'
					gVal=g2(jj);
					[~,vix]=min(abs(gVal-gcL));				
					p1(2)=plot(C.chi,C.elev,'-.','Color',cMap(vix,:));
					if lab
						[mc,ix]=max(C.chi);
						me=C.elev(ix);
						text(mc+mc*0.01,me,[ab1 ' ' num2str(id2(jj))],'Color',cMap(vix,:));
					end
				end	
				hold off

				subplot(2,1,2);
				hold on
				switch color_by
				case 'side'
					p2(2)=plot(C.distance,C.elev,'Color','r');
					if lab
						[md,ix]=max(C.distance);
						me=C.elev(ix);
						text(md+md*0.01,me,[ab1 ' ' num2str(id2(jj))],'Color','r');
					end
				case 'relief'
					rVal=r2(jj);
					[~,vix]=min(abs(rVal-rcL));
					p2(2)=plot(C.distance,C.elev,'-.','Color',cMap(vix,:));
					if lab
						[mc,ix]=max(C.chi);
						me=C.elev(ix);
						text(mc+mc*0.01,me,[ab1 ' ' num2str(id2(jj))],'Color',cMap(vix,:));
					end
				case 'gradient'
					gVal=g2(jj);
					[~,vix]=min(abs(gVal-gcL));
					p2(2)=plot(C.distance,C.elev,'-.','Color',cMap(vix,:));	
					if lab
						[mc,ix]=max(C.chi);
						me=C.elev(ix);
						text(mc+mc*0.01,me,[ab1 ' ' num2str(id2(jj))],'Color',cMap(vix,:));
					end
				end
				hold off
			end

			LegendEnt={['Streams ' pos1 ' of Divide'],['Streams ' pos2 ' of Divide']};

			subplot(2,1,1);
			hold on
				xlabel('Chi')
				ylabel('Elevation (m)')
				title('Chi - Z')
				legend(p1,LegendEnt,'location','best');
				switch color_by
				case 'relief'
					caxis([min(rL) max(rL)]);
					c1=colorbar;
					ylabel(c1,'Relief');
				case 'gradient'
					caxis([min(gL) max(gL)]);
					c1=colorbar;
					ylabel(c1,'Gradient');
				end
			hold off

			subplot(2,1,2);
			hold on
				xlabel('Distance from Mouth (m)')
				ylabel('Elevation (m)')
				title('Long Profile')
				legend(p2,LegendEnt,'location','best');
				switch color_by
				case 'relief'
					caxis([min(rL) max(rL)]);
					c2=colorbar;
					ylabel(c2,'Relief');
				case 'gradient'
					caxis([min(gL) max(gL)]);
					c2=colorbar;
					ylabel(c2,'Gradient');
				end
			hold off

			drawnow

			if ~isempty(name_prefix)
				print(f1,[name_prefix '_Divide_' num2str(ii) '.eps'],'-depsc');
			end

			fig_num=fig_num+1;
		end
	end
end

function [p1,p2,ab1,ab2]=DivideOrient(x1,y1,x2,y2);
	x1=double(x1); y1=double(y1);
	x2=double(x2); y2=double(y2);
    f1=fit(x1,y1,'poly1');
    f2=fit(x2,y2,'poly1');

    ms=(f1.p1+f2.p1)/2;

    ang=rad2deg(atan(1/ms));
    if ang>-22.5 & ang<22.5
        if mean(x1)<mean(x2)
            p1='West'; ab1='W';
            p2='East'; ab2='E';
        else
            p1='East'; ab1='E';
            p2='West'; ab2='W';
        end
    elseif ang>=22.5 & ang<67.5
        if mean(x1)<mean(x2) & mean(y1)>mean(y2)
            p1='Northwest'; ab1='NW';
            p2='Southeast'; ab2='SE';
        else
            p1='Southeast'; ab1='SE';
            p2='Northwest'; ab2='NW';
        end
    elseif ang>=67.5 & ang<=90
        if mean(y1)<mean(y2)
            p1='South'; ab1='S';
            p2='North'; ab2='N';
        else
            p1='North'; ab1='N';
            p2='South'; ab2='S';
        end
    elseif ang>-67.5 & ang<=-22.5
        if mean(x1)>mean(x2) & mean(y1)>mean(y2)
            p1='Northeast'; ab1='NE';
            p2='Southwest'; ab2='SW';
        else
            p1='Southwest'; ab1='SW';
            p2='Northeast'; ab2='NE';
        end
    elseif ang>=-90 & ang<=-67.5
        if mean(y1)<mean(y2)
            p1='South'; ab1='S';
            p2='North'; ab2='N';
        else
            p1='North'; ab1='N';
            p2='South'; ab2='S';
        end
    end
end
