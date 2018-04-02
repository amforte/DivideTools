function AlongDividePlot(head_vals,varargin)
    % Function produces plots of values of the four divide metrics along a series of divide segments. Will function if only one divide segment
    % is defined, but is really only useful when multiple divide segments have been defined. Will produce three plots (1) a plot of mean and
    % standard deviation values for both sides of a divide, (2) delta values for each of the four metrics along with propagated uncertainty, and (3)
    % a pseuod normalized delta plot (e.g. positive delta values for all metrics imply the same direction of divide motion).
    %
    % If you use the result of this code in a publication, please cite Forte, A.M. & Whipple, K.X., In Review, Criteria and Tools for Determining
    % Drainage Divide Stability, submitted to EPSL.
    %
    % Required Inputs:
    %       head_vals - output of 'AcrossDivide' function
    % Optional Inputs:
    %       wl_method ['std_dev'] - switch to determine restrictiveness of criteria for defining whether a divide is stable, 'std_dev' uses the standard deviation to determine if the means of
    %               a give metric overlap, alternatively 'std_err' uses the standard error, 'bootstrap' uses the 95% confidence interval from a normal bootstrap statistic, and 'ttest' uses a
    %               paired t-test to assess whether the means overlap.
    %		prefix ['DIV'] - text prefix for naming divide segments in the resulting plots
    %		side_1_name ['Side 1'] - text to name the side of the divide with odd number pour points, e.g. 'North' or 'Interior of Plateau'
    %		side_2_name ['Side 2'] - text to name the side of the divide with even number pour points, e.g. 'South' or 'Exterior of Plateau'
    %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Function Written by Adam M. Forte - Last Revised Spring 2018 %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	p=inputParser;
    p.FunctionName = 'AlongDividePlot';
    addRequired(p,'head_vals',@(x) isnumeric(x) & size(x,2)==8 );

    addParamValue(p,'wl_method','std_dev',@(x) ischar(validatestring(x,{'std_dev','std_err','bootstrap','ttest'})));
    addParamValue(p,'prefix','DIV',@(x) ischar(x));
    addParamValue(p,'side_1_name','Side 1', @(x) ischar(x));
    addParamValue(p,'side_2_name','Side 2', @(x) ischar(x));   

    parse(p,head_vals,varargin{:});
    head_vals=p.Results.head_vals;

    wl_method=p.Results.wl_method;
    pf=p.Results.prefix;
    s1=p.Results.side_1_name;
    s2=p.Results.side_2_name;

	% Break into sepearte variables
	chX=head_vals(:,1); chY=head_vals(:,2); 
	Elev=head_vals(:,3); Grad=head_vals(:,4); Rlf=head_vals(:,5); Chi=head_vals(:,6); 
	Div=head_vals(:,7); ChID=head_vals(:,8);

	num_sides=unique(Div);
	num_divides=max(num_sides)/2;

    odds=[1:2:(num_divides*2)-1];
    evens=[2:2:(num_divides*2)];

	mo=zeros(num_divides,17);
	tick_names=cell(1,num_divides);
	side1=cell(1,num_divides);
	side2=cell(1,num_divides);

	dlt=zeros(num_divides,9);

	tick_E_names=cell(1,num_divides);
	tick_C_names=cell(1,num_divides);
	tick_R_names=cell(1,num_divides);
	tick_G_names=cell(1,num_divides);

    for ii=1:num_divides
		idx1=Div==odds(ii);
		idx2=Div==evens(ii);

		tick_names{1,ii}=[pf num2str(ii)];

		warning off
		[do,p1,p2,Ornt]=DivideOrient(chX(idx1),chY(idx1),chX(idx2),chY(idx2));
		warning on

		side1{1,ii}=p1;
		side2{1,ii}=p2;

		E1=Elev(idx1);
		G1=Grad(idx1);
		R1=Rlf(idx1);
		C1=Chi(idx1);

		E2=Elev(idx2);
		G2=Grad(idx2);
		R2=Rlf(idx2);
		C2=Chi(idx2);

		ME1=mean(E1); MG1=mean(G1); MR1=mean(R1); MC1=mean(C1);
        ME2=mean(E2); MG2=mean(G2); MR2=mean(R2); MC2=mean(C2);

        switch wl_method
        case 'std_dev'
            stdE1=std(E1); stdG1=std(G1); stdR1=std(R1); stdC1=std(C1);
            stdE2=std(E2); stdG2=std(G2); stdR2=std(R2); stdC2=std(C2);
        case 'ttest'
            stdE1=std(E1); stdG1=std(G1); stdR1=std(R1); stdC1=std(C1);
            stdE2=std(E2); stdG2=std(G2); stdR2=std(R2); stdC2=std(C2);
        case 'std_err'
            stdE1=std(E1)/sqrt(numel(E1)); stdG1=std(G1)/sqrt(numel(G1)); stdR1=std(R1)/sqrt(numel(R1)); stdC1=std(C1)/sqrt(numel(C1));
            stdE2=std(E2)/sqrt(numel(E2)); stdG2=std(G2)/sqrt(numel(G2)); stdR2=std(R2)/sqrt(numel(R2)); stdC2=std(C2)/sqrt(numel(C2));
        case 'bootstrap'
            stdE1=bootCI(E1); stdG1=bootCI(G1); stdR1=bootCI(R1); stdC1=bootCI(C1);
            stdE2=bootCI(E2); stdG2=bootCI(G2); stdR2=bootCI(R2); stdC2=bootCI(C2);            
        end

        mo(ii,1)=ii;

        mo(ii,2)=ME1; mo(ii,3)=stdE1; mo(ii,4)=ME2; mo(ii,5)=stdE2;
        mo(ii,6)=MG1; mo(ii,7)=stdG1; mo(ii,8)=MG2; mo(ii,9)=stdG2;
        mo(ii,10)=MR1; mo(ii,11)=stdR1; mo(ii,12)=MR2; mo(ii,13)=stdR2;
        mo(ii,14)=MC1; mo(ii,15)=stdC1; mo(ii,16)=MC2; mo(ii,17)=stdC2;

        dlt(ii,1)=ii;
		[WL]=WinnersLosers(E1,E2,G1,G2,R1,R2,C1,C2,wl_method);
		dlt(ii,2)=WL.E_DLT;
		dlt(ii,3)=WL.G_DLT;
		dlt(ii,4)=WL.R_DLT;
		dlt(ii,5)=WL.C_DLT;
		dlt(ii,6)=WL.E_UNC;
		dlt(ii,7)=WL.G_UNC;
		dlt(ii,8)=WL.R_UNC;
		dlt(ii,9)=WL.C_UNC;

		[str]=OutputParser(WL,Ornt);
		ESTR=str{1};
		RSTR=str{2};
		GSTR=str{3};
		CSTR=str{4};

		tick_E_names{1,ii}=[pf num2str(ii) '\newline' ESTR];
		tick_R_names{1,ii}=[pf num2str(ii) '\newline' RSTR];
		tick_G_names{1,ii}=[pf num2str(ii) '\newline' GSTR];
		tick_C_names{1,ii}=[pf num2str(ii) '\newline' CSTR];

    end

    f1=figure(1);
    clf
    set(f1,'Units','normalized','Position',[0.1 0.1 0.50 0.75],'renderer','painters','PaperPositionMode','auto');

    subplot(4,1,1)
    hold on 
    errorbar(mo(:,1),mo(:,2),mo(:,3),'-k','CapSize',4); 
	errorbar(mo(:,1),mo(:,4),mo(:,5),'-r','CapSize',4);  
    scatter(mo(:,1),mo(:,2),20,'k','filled'); 
	scatter(mo(:,1),mo(:,4),20,'r','filled');

	for ii=1:num_divides

		if mo(ii,2)>=mo(ii,4)
			pos1=mo(ii,2)+(2*mo(ii,3));
			pos2=mo(ii,4)-(2*mo(ii,5));
		else
			pos1=mo(ii,2)-(2*mo(ii,3));
			pos2=mo(ii,4)+(2*mo(ii,5));
		end

		text(mo(ii,1),pos1,side1{1,ii},'Color','k');
		text(mo(ii,1),pos2,side2{1,ii},'Color','r');

		pos1L(ii,1)=pos1;
		pos2L(ii,1)=pos2;
	end

	posL=vertcat(pos1L,pos2L);
	ylim([min(posL)-0.1*min(posL) max(posL)+0.1*max(posL)]);

	title('Elevation');
	xlim([0,num_divides+1]);
	xticks([1:1:num_divides]);
	xticklabels(tick_names);
	hold off

    subplot(4,1,3)
    hold on
    errorbar(mo(:,1),mo(:,6),mo(:,7),'-k','CapSize',4); 
	errorbar(mo(:,1),mo(:,8),mo(:,9),'-r','CapSize',4);  
    scatter(mo(:,1),mo(:,6),20,'k','filled'); 
	scatter(mo(:,1),mo(:,8),20,'r','filled');

	for ii=1:num_divides

		if mo(ii,6)>=mo(ii,8)
			pos1=mo(ii,6)+(2*mo(ii,7));
			pos2=mo(ii,8)-(2*mo(ii,9));
		else
			pos1=mo(ii,6)-(2*mo(ii,7));
			pos2=mo(ii,8)+(2*mo(ii,9));
		end

		text(mo(ii,1),pos1,side1{1,ii},'Color','k');
		text(mo(ii,1),pos2,side2{1,ii},'Color','r');

		pos1L(ii,1)=pos1;
		pos2L(ii,1)=pos2;
	end

	posL=vertcat(pos1L,pos2L);
	ylim([min(posL)-0.1*min(posL) max(posL)+0.1*max(posL)]);

	title('Gradient');
	xlim([0,num_divides+1]);
	xticks([1:1:num_divides]);
	xticklabels(tick_names);
	hold off	

    subplot(4,1,4)
    hold on
    errorbar(mo(:,1),mo(:,10),mo(:,11),'-k','CapSize',4); 
	errorbar(mo(:,1),mo(:,12),mo(:,13),'-r','CapSize',4);  
    scatter(mo(:,1),mo(:,10),15,'k','filled'); 
	scatter(mo(:,1),mo(:,12),15,'r','filled');

	for ii=1:num_divides

		if mo(ii,10)>=mo(ii,12)
			pos1=mo(ii,10)+(2*mo(ii,11));
			pos2=mo(ii,12)-(2*mo(ii,13));
		else
			pos1=mo(ii,10)-(2*mo(ii,11));
			pos2=mo(ii,12)+(2*mo(ii,13));
		end

		text(mo(ii,1),pos1,side1{1,ii},'Color','k');
		text(mo(ii,1),pos2,side2{1,ii},'Color','r');

		pos1L(ii,1)=pos1;
		pos2L(ii,1)=pos2;
	end

	posL=vertcat(pos1L,pos2L);
	ylim([min(posL)-0.1*min(posL) max(posL)+0.1*max(posL)]);

	title('Relief');
	xlim([0,num_divides+1]);
	xticks([1:1:num_divides]);
	xticklabels(tick_names);
	hold off

    subplot(4,1,2)
    hold on
    errorbar(mo(:,1),mo(:,14),mo(:,15),'-k','CapSize',4); 
	errorbar(mo(:,1),mo(:,16),mo(:,17),'-r','CapSize',4);  
    scatter(mo(:,1),mo(:,14),15,'k','filled'); 
	scatter(mo(:,1),mo(:,16),15,'r','filled');

	for ii=1:num_divides

		if mo(ii,14)>=mo(ii,16)
			pos1=mo(ii,14)+(2*mo(ii,15));
			pos2=mo(ii,16)-(2*mo(ii,17));
		else
			pos1=mo(ii,14)-(2*mo(ii,15));
			pos2=mo(ii,16)+(2*mo(ii,17));
		end

		text(mo(ii,1),pos1,side1{1,ii},'Color','k');
		text(mo(ii,1),pos2,side2{1,ii},'Color','r');

		pos1L(ii,1)=pos1;
		pos2L(ii,1)=pos2;
	end

	posL=vertcat(pos1L,pos2L);
	ylim([min(posL)-0.1*min(posL) max(posL)+0.1*max(posL)]);

	title('\chi');
	xlim([0,num_divides+1]);
	xticks([1:1:num_divides]);
	xticklabels(tick_names);
	hold off

    f2=figure(2);
    clf
    set(f2,'Units','normalized','Position',[0.1 0.1 0.50 0.75],'renderer','painters','PaperPositionMode','auto');

    subplot(4,1,1)
    hold on
    plot([0,num_divides+1],[0,0],'--k');
    errorbar(dlt(:,1),dlt(:,2),dlt(:,6),'-r');
    scatter(dlt(:,1),dlt(:,2),20,'r','filled');
	xlim([0,num_divides+1]);
	xticks([1:1:num_divides]);
	xticklabels(tick_E_names);
	ylabel('\Delta Elevation')
	hold off

    subplot(4,1,2)
    hold on
    plot([0,num_divides+1],[0,0],'--k');
    errorbar(dlt(:,1),dlt(:,5),dlt(:,9),'-r');
    scatter(dlt(:,1),dlt(:,5),20,'r','filled');
	xlim([0,num_divides+1]);
	xticks([1:1:num_divides]);
	xticklabels(tick_C_names);
	ylabel('\Delta \chi')
	hold off

    subplot(4,1,3)
    hold on
	plot([0,num_divides+1],[0,0],'--k');
    errorbar(dlt(:,1),dlt(:,3),dlt(:,7),'-r');
    scatter(dlt(:,1),dlt(:,3),20,'r','filled');
	xlim([0,num_divides+1]);
	xticks([1:1:num_divides]);
	xticklabels(tick_G_names);
	ylabel('\Delta Gradient')
	hold off

    subplot(4,1,4)
    hold on
    plot([0,num_divides+1],[0,0],'--k');
    errorbar(dlt(:,1),dlt(:,4),dlt(:,8),'-r');
    scatter(dlt(:,1),dlt(:,4),20,'r','filled');
	xlim([0,num_divides+1]);
	xticks([1:1:num_divides]);
	xticklabels(tick_R_names);
	ylabel('\Delta Relief')
	hold off

    f3=figure(3);
    clf
    set(f3,'Units','normalized','Position',[0.1 0.1 0.50 0.75],'renderer','painters','PaperPositionMode','auto');

    subplot(4,1,1)
    hold on
    plot([0,num_divides+1],[0,0],'--k');
    errorbar(dlt(:,1),-1*dlt(:,2),dlt(:,6),'-r');
    scatter(dlt(:,1),-1*dlt(:,2),20,'r','filled');
	xlim([0,num_divides+1]);
	xticks([1:1:num_divides]);
	xticklabels(tick_names);
	ax1=gca;
	yticks([ax1.YLim(1),0,ax1.YLim(2)]);
	yticklabels({['Divide Moves Toward\newline' s1],'Divide is Stable',['Divide Moves Toward\newline' s2]})
	title('\Delta Elevation')
	hold off

    subplot(4,1,2)
    hold on
    plot([0,num_divides+1],[0,0],'--k');
    errorbar(dlt(:,1),-1*dlt(:,5),dlt(:,9),'-r');
    scatter(dlt(:,1),-1*dlt(:,5),20,'r','filled');
	xlim([0,num_divides+1]);
	xticks([1:1:num_divides]);
	xticklabels(tick_names);
	ax1=gca;
	yticks([ax1.YLim(1),0,ax1.YLim(2)]);
	yticklabels({['Divide Moves Toward\newline' s1],'Divide is Stable',['Divide Moves Toward\newline' s2]})
	title('\Delta \chi')
	hold off

    subplot(4,1,3)
    hold on
	plot([0,num_divides+1],[0,0],'--k');
    errorbar(dlt(:,1),dlt(:,3),dlt(:,7),'-r');
    scatter(dlt(:,1),dlt(:,3),20,'r','filled');
	xlim([0,num_divides+1]);
	xticks([1:1:num_divides]);
	xticklabels(tick_names);
	ax1=gca;
	yticks([ax1.YLim(1),0,ax1.YLim(2)]);
	yticklabels({['Divide Moves Toward\newline' s1],'Divide is Stable',['Divide Moves Toward\newline' s2]})
	title('\Delta Gradient')
	hold off

    subplot(4,1,4)
    hold on
    plot([0,num_divides+1],[0,0],'--k');
    errorbar(dlt(:,1),dlt(:,4),dlt(:,8),'-r');
    scatter(dlt(:,1),dlt(:,4),20,'r','filled');
	xlim([0,num_divides+1]);
	xticks([1:1:num_divides]);
	xticklabels(tick_names);
	ax1=gca;
	yticks([ax1.YLim(1),0,ax1.YLim(2)]);
	yticklabels({['Divide Moves Toward\newline' s1],'Divide is Stable',['Divide Moves Toward\newline' s2]})
	title('\Delta Relief')
	hold off

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
            p1='W';
            Ornt.p1=270;
            p2='E';
            Ornt.p2=90;
        else
            p1='E';
            Ornt.p1=90;
            p2='W';
            Ornt.p2=270;
        end
    elseif ang>=22.5 & ang<67.5
        do='NE-SW';
        Ornt.d=2;
        if mean(x1)<mean(x2) & mean(y1)>mean(y2)
            p1='NW';
            Ornt.p1=315;
            p2='SE';
            Ornt.p2=135;
        else
            p1='SE';
            Ornt.p1=135;
            p2='NW';
            Ornt.p2=315;
        end
    elseif ang>=67.5 & ang<=90
        do='E-W';
        Ornt.d=3;
        if mean(y1)<mean(y2)
            p1='S';
            Ornt.p1=180;
            p2='N';
            Ornt.p2=0;
        else
            p1='N';
            Ornt.p1=0;
            p2='S';
            Ornt.p2=180;
        end
    elseif ang>-67.5 & ang<=-22.5
        do='SE-NW';
        Ornt.d=4;
        if mean(x1)>mean(x2) & mean(y1)>mean(y2)
            p1='NE';
            Ornt.p1=45;
            p2='SW';
            Ornt.p2=225;
        else
            p1='SW';
            Ornt.p1=225;
            p2='NE';
            Ornt.p2=45;
        end
    elseif ang>=-90 & ang<=-67.5
        do='E-W';
        Ornt.d=3;
        if mean(y1)<mean(y2)
            p1='S';
            Ornt.p1=180;
            p2='N';
            Ornt.p2=0;
        else
            p1='N';
            Ornt.p1=0;
            p2='S';
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

case 'ttest'
        WL=struct;

        [Et,~]=ttest2(E1,E2,'Vartype','unequal');
        [Gt,~]=ttest2(G1,G2,'Vartype','unequal');  
        [Rt,~]=ttest2(R1,R2,'Vartype','unequal');  
        [Ct,~]=ttest2(C1,C2,'Vartype','unequal');    

        if ME1<ME2 && Et~=0;
            WL.E_AV1=1; % Aggressor
            WL.E_AV2=2; % Victim
        elseif ME1>ME2 && Et~=0;
            WL.E_AV1=2;
            WL.E_AV2=1;
        else
            WL.E_AV1=3; % Stable
            WL.E_AV2=3;
        end
            
        if MG1>MG2 && Gt~=0;
            WL.G_AV1=1;
            WL.G_AV2=2;
        elseif MG1<MG2 && Gt~=0;
            WL.G_AV1=2;
            WL.G_AV2=1;
        else
            WL.G_AV1=3;
            WL.G_AV2=3;
        end

        if MR1>MR2 && Rt~=0;
            WL.R_AV1=1;
            WL.R_AV2=2;
        elseif MR1<MR2 && Rt~=0;
            WL.R_AV1=2;
            WL.R_AV2=1;
        else
            WL.R_AV1=3;
            WL.R_AV2=3;
        end

        if MC1<MC2 && Ct~=0;
            WL.C_AV1=1; % Aggressor
            WL.C_AV2=2; % Victim
        elseif MC1>MC2 && Ct~=0;
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

    case 'bootstrap'
        stdE1=bootCI(E1); stdG1=bootCI(G1); stdR1=bootCI(R1); stdC1=bootCI(C1);
        stdE2=bootCI(E2); stdG2=bootCI(G2); stdR2=bootCI(R2); stdC2=bootCI(C2);  

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
            E_str='S';
        elseif Ornt.p1==45
            E_str='SW';
        elseif Ornt.p1==90
            E_str='W';
        elseif Ornt.p1==135
            E_str='NW';
        elseif Ornt.p1==180
            E_str='N';
        elseif Ornt.p1==225
            E_str='NE';
        elseif Ornt.p1==270
            E_str='E';
        elseif Ornt.p1==315
            E_str='SE';
        end
    elseif WL.E_AV2==1
        if Ornt.p2==0
            E_str='S';
        elseif Ornt.p2==45
            E_str='SW';
        elseif Ornt.p2==90
            E_str='W';
        elseif Ornt.p2==135
            E_str='NW';
        elseif Ornt.p2==180
            E_str='N';
        elseif Ornt.p2==225
            E_str='NE';
        elseif Ornt.p2==270
            E_str='E';
        elseif Ornt.p2==315
            E_str='SE';
        end
    elseif WL.E_AV1==3
        E_str='STBL';
    end

    if WL.G_AV1==1
        if Ornt.p1==0
            G_str='S';
        elseif Ornt.p1==45
            G_str='SW';
        elseif Ornt.p1==90
            G_str='W';
        elseif Ornt.p1==135
            G_str='NW';
        elseif Ornt.p1==180
            G_str='N';
        elseif Ornt.p1==225
            G_str='NE';
        elseif Ornt.p1==270
            G_str='E';
        elseif Ornt.p1==315
            G_str='SE';
        end
    elseif WL.G_AV2==1
        if Ornt.p2==0
            G_str='S';
        elseif Ornt.p2==45
            G_str='SW';
        elseif Ornt.p2==90
            G_str='W';
        elseif Ornt.p2==135
            G_str='NW';
        elseif Ornt.p2==180
            G_str='N';
        elseif Ornt.p2==225
            G_str='NE';
        elseif Ornt.p2==270
            G_str='E';
        elseif Ornt.p2==315
            G_str='SE';
        end
    elseif WL.G_AV1==3
        G_str='STBL';
    end

    if WL.R_AV1==1
        if Ornt.p1==0
            R_str='S';
        elseif Ornt.p1==45
            R_str='SW';
        elseif Ornt.p1==90
            R_str='W';
        elseif Ornt.p1==135
            R_str='NW';
        elseif Ornt.p1==180
            R_str='N';
        elseif Ornt.p1==225
            R_str='NE';
        elseif Ornt.p1==270
            R_str='E';
        elseif Ornt.p1==315
            R_str='SE';
        end
    elseif WL.R_AV2==1
        if Ornt.p2==0
            R_str='S';
        elseif Ornt.p2==45
            R_str='SW';
        elseif Ornt.p2==90
            R_str='W';
        elseif Ornt.p2==135
            R_str='NW';
        elseif Ornt.p2==180
            R_str='N';
        elseif Ornt.p2==225
            R_str='NE';
        elseif Ornt.p2==270
            R_str='E';
        elseif Ornt.p2==315
            R_str='SE';
        end
    elseif WL.R_AV1==3
        R_str='STBL';
    end

    if WL.C_AV1==1
        if Ornt.p1==0
            C_str='S';
        elseif Ornt.p1==45
            C_str='SW';
        elseif Ornt.p1==90
            C_str='W';
        elseif Ornt.p1==135
            C_str='NW';
        elseif Ornt.p1==180
            C_str='N';
        elseif Ornt.p1==225
            C_str='NE';
        elseif Ornt.p1==270
            C_str='E';
        elseif Ornt.p1==315
            C_str='SE';
        end
    elseif WL.C_AV2==1
        if Ornt.p2==0
            C_str='S';
        elseif Ornt.p2==45
            C_str='SW';
        elseif Ornt.p2==90
            C_str='W';
        elseif Ornt.p2==135
            C_str='NW';
        elseif Ornt.p2==180
            C_str='N';
        elseif Ornt.p2==225
            C_str='NE';
        elseif Ornt.p2==270
            C_str='E';
        elseif Ornt.p2==315
            C_str='SE';
        end
    elseif WL.C_AV1==3
        C_str='STBL';
    end

    str={E_str;R_str;G_str;C_str};
end

function [CI]=bootCI(data)
    m=zeros(1000,1);
    for ii=1:1000
        d=datasample(data,numel(data));
        m(ii)=mean(d);
    end

    CI=std(m);

    m=sort(m,'ascend');

    [f,x]=ecdf(m);

    cil_ix=find(f>0.025,1,'first');
    cih_ix=find(f<0.975,1,'first');

    cil=x(cil_ix);
    cih=x(cih_ix);

    CI=mean([abs(cil-mean(m)),abs(cih-mean(m))]);
end



