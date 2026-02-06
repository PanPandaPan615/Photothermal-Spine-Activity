clc, clear, close all

% Load Data
load("Phys Temp.mat"),    PhysTemp = ExpData; 
load("Room Temp.mat"),    RoomTemp = ExpData.PulseTrain; 
RoomTemp = rmfield(RoomTemp,'Post');
load("Pharmacology.mat"), Pharma = ExpData.PulseTrain;
Pharma = rmfield(Pharma,'Post');
clear ExpData

DoseType = {'Single', 'Multi'};
artifactpts = 31:5:51;

Colors = [0 0 0; 
         .5 .5 0; 
          0 .5 0;];
yellow_lut = [linspace(0,1,256)', linspace(0,1,256)',zeros(256,1)];

% Separate Room Temp experiments into experimental categories
mask = contains({RoomTemp.Single.filename},'30A');
RoomTemp.Subthreshold.Single = RoomTemp.Single(mask);
RoomTemp.Subthreshold.Multi = RoomTemp.Multi(mask);
RoomTemp.Threshold.Single = RoomTemp.Single(~mask);
RoomTemp.Threshold.Multi = RoomTemp.Multi(~mask);

% Good ol Peak counting. This is the way. [labels,locs height,prom,width]
% Soma = [loction, height, prominance, width ]
% Dendrites = [DendriteLabel, loction, height, prominance, width, bAP]
% Spines = [SpineLabel, ParentDendrite, loction, height, prominance, width, bAP, Reactive, NeuronParent]

PhysTemp.PulseTrain.Single  = ExtractPks(PhysTemp.PulseTrain.Single,752,artifactpts(1));
PhysTemp.PulseTrain.Multi  = ExtractPks(PhysTemp.PulseTrain.Multi,900,artifactpts);
PhysTemp.SinglePulse.Single  = ExtractPks(PhysTemp.SinglePulse.Single,752,artifactpts(1));
PhysTemp.SinglePulse.Multi  = ExtractPks(PhysTemp.SinglePulse.Multi,900,artifactpts);

RoomTemp.Subthreshold.Single  = ExtractPks(RoomTemp.Subthreshold.Single,752,artifactpts(1));
RoomTemp.Subthreshold.Multi  = ExtractPks(RoomTemp.Subthreshold.Multi,900,artifactpts);
RoomTemp.Threshold.Single  = ExtractPks(RoomTemp.Threshold.Single,752,artifactpts(1));
RoomTemp.Threshold.Multi  = ExtractPks(RoomTemp.Threshold.Multi,900,artifactpts);

Pharma.Single = ExtractPks(Pharma.Single,752,artifactpts(1));
Pharma.Multi = ExtractPks(Pharma.Multi,899,artifactpts);

% Append Data
PhysTemp.PulseTrain = AppendData(PhysTemp.PulseTrain);clc
PhysTemp.SinglePulse = AppendData(PhysTemp.SinglePulse);
RoomTemp.Subthreshold = AppendData(RoomTemp.Subthreshold);
RoomTemp.Threshold = AppendData(RoomTemp.Threshold);

% Peak Stats
PhysTemp.PulseTrain.SpinePkStats  = SpinePeakStats(PhysTemp.PulseTrain);
PhysTemp.SinglePulse.SpinePkStats = SpinePeakStats(PhysTemp.SinglePulse);
RoomTemp.Subthreshold.SpinePkStats = SpinePeakStats(RoomTemp.Subthreshold);
RoomTemp.Threshold.SpinePkStats = SpinePeakStats(RoomTemp.Threshold);
Pharma.SpinepkStats = SpinePeakStats(Pharma);


%% Panel 3A
close all
    ExpTyp = PhysTemp.SinglePulse.Single;    
    neuron = 4;
figure('Theme','light'); 
    subplot(8,1,1:5), 
    spines = RemoveDC(ExpTyp(neuron).Spines,200);
    t = (1:length(spines))/7.5;
    imagesc(spines'), colormap("hot"), 
    axpos =  get(gca,"Position");
    c =colorbar;  c.Label.String = '\DeltaF/F';  c.Location = 'manual'; 
    c.Position = [axpos(1)+axpos(3)+.01 axpos(2) .02 axpos(4)]; 
    xline(30*7.5,'LineWidth',2,LineStyle="--",Alpha=1,Color='w')
    set(gca, 'XTick', [], 'YTick', [],'FontSize',24)

    ylabel('Denditic Spines')

    subplot(8,1,6)
    mask = ~ExpTyp(neuron).SpinePks(:,end-2);
    xline(ExpTyp(neuron).SpinePks(mask,3),'Color','r','LineWidth',2,Alpha=1 ); 
    xline(ExpTyp(neuron).SomaPks(:,1),'LineWidth',2,Alpha=1 ); 
    xlim([t(1) t(end)]), xline(30,'LineWidth',2,LineStyle="--",Alpha=1 )
    set(gca, 'XTick', [], 'YTick', [])

    subplot(8,1,7:8)
    soma = RemoveDC(ExpTyp(neuron).Soma,200);
    plot(t,soma,'Color','k','LineWidth',2),
    xlim([t(1) t(end)]), xline(30,'LineWidth',2,LineStyle="--",Alpha=1 )
    set(gca, 'FontSize',24)    
    ylabel('\DeltaF/F'), xlabel('time (sec)')

%% Figure 3 comparing Soma, and Spine Responses
close all

SomaLabels = []; SomaMean =[]; SomaPhenotype= []; SpineLabels =[]; 
SpineMean =[]; SpinePhenotype= []; SpikePercentage = [];
duration = 75;
 % Mean Plot responses
  for i = 1:4
     if i ==1,     PulseType = PhysTemp.PulseTrain;      ExposureType = PulseType.SingleAppened;
     stimpt = round(30*7.5); thresh = stimpt+duration;
     ReactiveFraction = PulseType.SpinePkStats.SingleRespondingFraction(:,2);
     elseif i ==2, PulseType = PhysTemp.SinglePulse;     ExposureType = PulseType.SingleAppened; 
     ReactiveFraction = PulseType.SpinePkStats.SingleRespondingFraction(:,2);
     elseif i ==3, PulseType = PhysTemp.PulseTrain;      ExposureType = PulseType.MultiAppened;
     stimpt = round(50*7.5); thresh = stimpt+duration;
     ReactiveFraction = PulseType.SpinePkStats.MultiRespondingFraction(:,2);
     else,         PulseType = PhysTemp.SinglePulse;     ExposureType = PulseType.MultiAppened;
     ReactiveFraction = PulseType.SpinePkStats.MultiRespondingFraction(:,2);
     end
     
     soma = RemoveDC(ExposureType.Soma,225);
     spines = RemoveDC(ExposureType.Spines,225);
     SomaMean = [SomaMean; mean(soma(stimpt:thresh,:))'];
     SpineMean = [SpineMean; mean(spines(stimpt:thresh,:))']; 
       
         
     mask = median(sum(soma(stimpt:thresh,:))')<sum(soma(stimpt:thresh,:))'; % Reactive neurons are above the median of the AUC 
     SomaLabels = [SomaLabels; ones(size(mask))*i];
     SomaPhenotype = logical([SomaPhenotype; mask]);

     SpikePercentage = [SpikePercentage; ReactiveFraction];
     pv = ranksum(ReactiveFraction(mask),ReactiveFraction(~mask));
        
 % Panel C
    mask = (PulseType.Labels.Structure=='Spine') & ismember(PulseType.Labels.NeuronID,find(mask));
    temp = (PulseType.Labels.Structure=='Soma') | (PulseType.Labels.Structure=='Dendrite');
    mask(temp) = [];
    SpineLabels = [SpineLabels; ones(size(mask))*i];
    SpinePhenotype = logical([SpinePhenotype; mask]);
  end

% 3B
    Labels = SomaLabels+SomaPhenotype*4;
    TickLabels = {'PTSE_(_-_)' 'SPSE_(_-_)' 'PTME_(_-_)' 'SPME_(_-_)'...
                  'PTSE_(_+_)' 'SPSE_(_+_)' 'PTME_(_+_)'  'SPME_(_+_)'};
     
    Colors = [repmat([.5 .5 .5],[4 1]); repmat([1 .5 .5],[4 1])];
    f1 = figure(Theme='light');
    boxchart(Labels, SomaMean,'GroupByColor',Labels,'ColorGroupLayout',...
        'overlaid', 'BoxEdgeColor','k','MarkerSize',10, 'MarkerStyle',...
         '.','MarkerColor','k','BoxFaceAlpha', 1);     
    ax = gca; ax.FontSize = 24; ax.TickDir = 'none'; 
       
    colororder(Colors),         ylabel('\DeltaF/F'),   
    set(gca, 'XTick', 1:max(Labels), 'XTickLabels', TickLabels)
    ax = gca; ax.FontSize = 24; ax.TickDir = 'none'; theme('light')

     % Anderson-Darling normality test
    for i = 1:max(Labels)
        data = SomaMean(Labels==i);
       isnormal(i) = adtest(data);
    end, isnormal
    
    % P Vals and heat map
    dFSoma = PairwiseTest(SomaMean,Labels,TickLabels,'rank');
    mask = dFSoma.P_Value<0.05;
    figure(f1),    hold on
    sigstar(dFSoma.Pairs(mask), dFSoma.P_Value(mask)); hold off; 
    legend('Reactive Soma', 'Unreactive Soma')
    
    for i = 5:8
        mask = Labels==i;
        M = mean(SomaMean(mask))*100;
        sem = std(SomaMean(mask))/sqrt(sum(mask))*100;
        sprintf('%s: %.2f+/-%.2f', TickLabels{i}, M, sem)
    end

% Panel 3C 
    Labels = SpineLabels+SpinePhenotype*4;
   
figure('Theme','light');
   swarmchart(SpineLabels(~SpinePhenotype),SpineMean(~SpinePhenotype),5,'black','filled'), hold on
   swarmchart(SpineLabels(SpinePhenotype),SpineMean(SpinePhenotype),5,'red','filled'), hold off
   set(gca, 'XTick', 1:4, 'XTickLabels', {'PTSE' 'SPSE' 'PTME'  'SPME'})
   ylabel("\DeltaF/F")
   ax = gca; ax.FontSize = 30; ax.TickDir = 'none';
   legend('Unreactive (-)','Reactive (+)',"Location", "northwest",...
           "Color",'none','EdgeColor','none'); 
    
    % Kolmgrov-smirnov normality test & effect size
    for i = 1:max(SpineLabels)
        spines = SpineLabels==i;
        isnormal(i) = kstest(SpineMean(spines));
        x = spines & SpinePhenotype;
        x = SpineMean(x);
        y = spines & ~SpinePhenotype;
        y = SpineMean(y);
        meanEffectSize(x,y,"Effect","cohen")
        d(i,1) = ans.Effect;
        CI = abs(ans.ConfidenceIntervals-ans.Effect);
        d(i,2) = mean(CI);
    end, isnormal
    
    [~,~,stats] =anova1(SpineMean,Labels,"off");
    c = multcompare(stats,'Display','off','CriticalValueType','bonferroni');
    pvVal = orderpvalues(c); 
    plotheatmap(TickLabels, pvVal);
    title('ANOVA P-Values')
    
    dFSpines = PairwiseTest(SpineMean,Labels,TickLabels,'ttest');


   for i = 5:8
        mask = Labels==i;
        M = mean(SpineMean(mask))*100;
        sem = std(SpineMean(mask))/sqrt(sum(mask))*100;
        sprintf('%s: %.2f+/-%.2f', TickLabels{i}, M, sem)
    end

% 3D Spiking activity striated by soma phenotype
close all
f1 = figure(Theme='light');
    Labels = SomaLabels+SomaPhenotype*4;
    boxchart(Labels, SpikePercentage*100,'GroupByColor',Labels,'ColorGroupLayout',...
        'overlaid', 'BoxEdgeColor','k','MarkerSize',10, 'MarkerStyle',...
         '.','MarkerColor','k','BoxFaceAlpha', .5),        
    colororder(Colors),         ylabel('% DS Responding/Neuron'),   ylim([0 100])
    set(gca, 'XTick', 1:max(Labels), 'XTickLabels', TickLabels,'FontSize',24)


    PhenoSpike = PairwiseTest(SpikePercentage,Labels,TickLabels,'rank');

    % Anderson-Darling normality test
    for i = 1:max(Labels)
        data = SpikePercentage(Labels==i);
        isnormal(i) = adtest(data);
    end, isnormal

    mask = PhenoSpike.P_Value<0.05;
    figure(f1),    hold on
    sigstar(PhenoSpike.Pairs(mask), PhenoSpike.P_Value(mask)); hold off; % Add significance markers


%% Figure 4
close all

% 4B Distances from Soma for each experimental protocol
Distances = horzcat([PhysTemp.PulseTrain.Single.SpineDistance],...
                    [PhysTemp.PulseTrain.Multi.SpineDistance],...
                    [PhysTemp.SinglePulse.Single.SpineDistance],...
                    [PhysTemp.SinglePulse.Multi.SpineDistance]);

% Labels for all Spines
Labels = horzcat(ones(1,length([PhysTemp.PulseTrain.Single.SpineEventMask])),...
                 ones(1,length([PhysTemp.PulseTrain.Multi.SpineEventMask]))*2,...
                 ones(1,length([PhysTemp.SinglePulse.Single.SpineEventMask]))*3,...
                 ones(1,length([PhysTemp.SinglePulse.Multi.SpineEventMask]))*4);

% Mask of activated spines
mask = logical(horzcat([PhysTemp.PulseTrain.Single.SpineEventMask],...
                    [PhysTemp.PulseTrain.Multi.SpineEventMask],...
                    [PhysTemp.SinglePulse.Single.SpineEventMask],...
                    [PhysTemp.SinglePulse.Multi.SpineEventMask]));


figure('Theme','light');      
set(gcf, 'Units', 'inches', 'Position', [5, 5, 5.5, 4]);
   swarmchart(Labels(~mask),Distances(~mask),10,'black','filled'), hold on
   swarmchart(Labels(mask),Distances(mask)',10,'red','filled'), hold off
   set(gca, 'XTick', 1:4, 'XTickLabels', {'PTSE' 'PTME' 'SPSE' 'SPME'})
   ylabel("Distance (\mum)")
   legend('Unreactive','Reactive',"Location", "northwest",...
           "Color",'none','EdgeColor','none'); 
   ax = gca; ax.FontSize = 24; ax.TickDir = 'none';

for i = 1:4
    imask = Labels == i;
    [r(i)] = corr(Distances(imask)', mask(imask)', 'Type', 'Pearson');  
end
   Labels = Labels+mask*4;

   TickLabels = {'PTSE (-)' 'PTME (-)' 'SPSE (-)' 'SPME (-)'...
                  'PTSE (+)' 'PTME  (+)' 'SPSE  (+)' 'SPME (+)'};

    [~,~,stats] = anova1(Distances, Labels,"off");
    c = multcompare(stats,'Display','off');
    p = orderpvalues(c); 
    plotheatmap(TickLabels, p)
    title('ANOVA Bonferroni Post-Hoc P-Values')

% Figure 4 C & D
    PTSS = PhysTemp.PulseTrain.SpinePkStats.Post_Single;
    PTMS = PhysTemp.PulseTrain.SpinePkStats.Post_Multi;
    SPSS = PhysTemp.SinglePulse.SpinePkStats.Post_Single;
    SPMS = PhysTemp.SinglePulse.SpinePkStats.Post_Multi;

    PkStats = vertcat(PTSS, PTMS, SPSS, SPMS);
    Labels = [ones(length(PTSS),1); ones(length(PTMS),1)*2; ...
              ones(length(SPSS),1)*3; ones(length(SPMS),1)*4];
    IRPks = horzcat(PkStats, Labels);
    
    PTSS = PhysTemp.PulseTrain.SpinePkStats.Single_bAP;
    PTMS = PhysTemp.PulseTrain.SpinePkStats.Multi_bAP;
    SPSS = PhysTemp.SinglePulse.SpinePkStats.Single_bAP;
    SPMS = PhysTemp.SinglePulse.SpinePkStats.Multi_bAP;

    PkStats = vertcat(PTSS, PTMS, SPSS, SPMS);
    Labels = [ones(length(PTSS),1); ones(length(PTMS),1)*2; ...
              ones(length(SPSS),1)*3; ones(length(SPMS),1)*4];
    bAPPks = horzcat(PkStats, Labels);
    mask = bAPPks(:,3)>30 & bAPPks(:,7)==1;
    post_bAPPks = bAPPks(mask,:);
    pre_bAPPks = bAPPks(~mask,:);
    

figure('Theme','light');    set(gcf, 'Units', 'inches', 'Position', [5, 5, 5.5, 4]);
   swarmchart(IRPks(:,end),IRPks(:,4),5,'red','filled'), hold on
   swarmchart(post_bAPPks(:,end),post_bAPPks(:,4),5,'black','filled'), hold off
   set(gca, 'XTick', 1:4, 'XTickLabels', {'PTSE' 'PTME' 'SPSE' 'SPME'})
   ylabel("\DeltaF/F (A.U.)")
   ax = gca; ax.FontSize = 24; ax.TickDir = 'none';
   legend('ICS','bAP',"Location", "best",...
           "Color",'none','EdgeColor','none'); 

   TickLabels = {'PTSE bAP' 'PTME bAP' 'SPSE bAP' 'SPME bAP'...
                 'PTSE ICS' 'PTME ICS' 'SPSE ICS' 'SPME ICS'};

    Data = [post_bAPPks(:,4); IRPks(:,4)];
    Labels = [post_bAPPks(:,end); IRPks(:,end)+4];
    
    for i = 1:max(Labels)
        neurons = Labels==i;
        tempstruct(i) = datastatsJake(Data(neurons));
        isnormal(i) = jbtest(Data(neurons));
    end,    dFstats = struct2table(tempstruct); clear tempstruct
    dFstats.Labels = TickLabels'; isnormal
    
    [~,~,stats] = anova1(Data, Labels,"off");
    c = multcompare(stats,'CriticalValueType','bonferroni','Display','off');
    p = orderpvalues(c); 
    plotheatmap(TickLabels, p);
    title('Kruskal-Wallis P-Values')
    disp(dFstats)


figure('Theme','light');   set(gcf, 'Units', 'inches', 'Position', [5, 5, 5.5, 4]);  
   swarmchart(IRPks(:,end),IRPks(:,5),5,'red','filled'), hold on
   swarmchart(post_bAPPks(:,end),post_bAPPks(:,5),5,'black','filled'), hold off
   set(gca, 'XTick', 1:4, 'XTickLabels', {'PTSE' 'PTME' 'SPSE' 'SPME'})
   ylabel("FWHM (Sec)")
   legend('ICS','bAP',"Location", "best",...
           "Color",'none','EdgeColor','none'); 
   ax = gca; ax.FontSize = 24; ax.TickDir = 'none';

   Data = [post_bAPPks(:,5); IRPks(:,5)];
   Labels = [post_bAPPks(:,end); IRPks(:,end)+4];


   for i = 1:max(Labels)
        neurons = Labels==i;
        tempstruct(i) = datastatsJake(Data(neurons));
        isnormal(i) = jbtest(Data(neurons));
    end,    FWHMstats = struct2table(tempstruct); clear tempstruct
    FWHMstats.Labels = TickLabels'; isnormal

    [~,~,stats] = anova1(Data, Labels,"off");
    c = multcompare(stats,'CriticalValueType','bonferroni','Display','off');
    p = orderpvalues(c);   
    plotheatmap(TickLabels, p);
    disp(FWHMstats)

% Figure 4E
PkCount = horzcat(PhysTemp.PulseTrain.SpinePkStats.SingleRespondingFraction,...
                  PhysTemp.PulseTrain.SpinePkStats.MultiRespondingFraction,...
                  PhysTemp.SinglePulse.SpinePkStats.SingleRespondingFraction,...
                  PhysTemp.SinglePulse.SpinePkStats.MultiRespondingFraction);
PkCount(:,[3 7]) = [];


    [j,k] = size(PkCount);
    PTPkCount = reshape(PkCount,  [j*k 1])*100;
    PTLabels = reshape(repmat(1:k,[j 1]),  size(PTPkCount));
    PTLabels(PTLabels==4) = 1; PTLabels(PTLabels==5) = 4; PTLabels(PTLabels==6) = 5;
    Colors = [0.5, 0.5, 0.5;  % Gray
                1, 0, 0;        % Bright Red
                0.8, 0, 0;      % Darker Red
                1, 0.6, 0;      % Light Orange
            0.8, 0.4, 0;];      % Darker Orange

figure('Theme','light'); set(gcf, 'Units', 'inches', 'Position', [1, 1, 7, 6]);
     boxchart(PTLabels, PTPkCount,'GroupByColor',PTLabels,'ColorGroupLayout',...
            'overlaid', 'BoxEdgeColor','k','MarkerSize',10, 'MarkerStyle',...
             '.','MarkerColor','k','BoxFaceAlpha', 1)
  %violinplot(PTLabels, PTPkCount,'GroupByColor',PTLabels,'ColorGroupLayout',...
   % 'overlaid','FaceAlpha', 1,'EdgeColor','k','DensityScale','count')
    colororder(Colors)
    %b = gbar(PkCount, Labels)
    ylabel("% DS Responding/Neuron"), ylim([0 100]), xlim([0.5 max(PTLabels)+0.5])
    set(gca, 'XTick', 1:5, 'XTickLabels', {'Pre Stim' 'PTSE_3_0' 'PTME_3_0' 'SPSE_3_0' 'SPME_3_0'})
    ax = gca; ax.FontSize = 24; ax.TickDir = 'none';
    
    [~,~,stats] = kruskalwallis(PTPkCount, PTLabels,'off');
    c = multcompare(stats,'Display','off',CriticalValueType='dunn-sidak');
    pairs = num2cell(c(:, 1:2),2);    mask = c(:,6)<0.05;
    sigstar(pairs(mask), c(mask,6)); hold off; % Add significance markers

    p = orderpvalues(c);   
    plotheatmap({'Pre Stim' 'PTSE' 'PTME' 'SPSE' 'SPME'}, p);
    title('Kruskal-Wallis  P-Values')

   for i = 1:max(PTLabels)
        neurons = PTLabels==i;
        tempstruct(i) = datastatsJake(PTPkCount(neurons));
        isnormal(i) = adtest(PTPkCount(neurons));
    end,    RespFraction = struct2table(tempstruct); clear tempstruct
    isnormal
    RespFraction.Labels = {'Pre Stim' 'PTSE' 'PTME' 'SPSE' 'SPME'}';
    disp(RespFraction)

% Figure 4F Room Temperature Data

PkCount = horzcat(RoomTemp.Subthreshold.SpinePkStats.SingleRespondingFraction,...
                    RoomTemp.Subthreshold.SpinePkStats.MultiRespondingFraction,...
                    RoomTemp.Threshold.SpinePkStats.SingleRespondingFraction,...
                    RoomTemp.Threshold.SpinePkStats.MultiRespondingFraction);

PkCount(:,[3, 5, 7]) = []; % pre exposure indexes for later time points

[j,k] = size(PkCount);
RTPkCount = reshape(PkCount,  [j*k 1])*100;
RTLabels = reshape(repmat(1:k,[j 1]),  size(RTPkCount));


titles =  {'Pre Stim' 'PTSE_S_u_b' 'PTME_S_u_b' 'SPSE_T_h_r_e_s_h' 'SPME_T_h_r_e_s_h'}; 
Colors = [0.5, 0.5, 0.5;  % Gray
            1, 0, 0;        % Bright Red
            0.8, 0, 0;      % Darker Red
            0.5, 0.8, 1;
            .1, .3, 1;];      

figure('Theme','light'); set(gcf, 'Units', 'inches', 'Position', [1, 1, 7, 6]);
     boxchart(RTLabels, RTPkCount,'GroupByColor',RTLabels,'ColorGroupLayout',...
            'overlaid', 'BoxEdgeColor','k','MarkerSize',10, 'MarkerStyle',...
             '.','MarkerColor','k','BoxFaceAlpha', 1)
    colororder(Colors)
    ylabel("% DS Responding/Neuron"), ylim([0 100])
    set(gca, 'XTick', 1:5, 'XTickLabels',titles)
    ax = gca; ax.FontSize = 24; ax.TickDir = 'none';

    [~,~,stats] = kruskalwallis(RTPkCount, RTLabels,'off');
    c = multcompare(stats,'Display','off',CriticalValueType='dunn-sidak');
    pairs = num2cell(c(:, 1:2),2);    mask = c(:,6)<0.05;
    sigstar(pairs(mask), c(mask,6)); hold off; % Add significance markers

    p = orderpvalues(c);
    plotheatmap(titles, p)
    title('Kruskal-Wallis P-Values')

   for i = 1:max(PTLabels)
        neurons = PTLabels==i;
        tempstruct(i) = datastatsJake(RTPkCount(neurons));
        isnormal(i) = adtest(RTPkCount(neurons));
    end,    RTRespFraction = struct2table(tempstruct); clear tempstruct
    RTRespFraction.Labels = titles';
    disp(RTRespFraction), isnormal
    
%% Figure 6 Physiological Mechanism
close all
PharmaGroups ={'Control','CaFree','AP5', 'CNQX', 'RR','Gd','TRP'};

CountLabels = []; PkCount = [];    spines = Pharma.Multi;
Locs = []; LocLabels = [];
Colors = orderedcolors("gem");
Colors([2, 5], :) = Colors([5, 2], :);

ex = [3 1 1 1 1 1 1 1 1];
for i = 1:length(PharmaGroups)

    mask = contains({Pharma.Multi.filename},PharmaGroups{i});
    temp = spines(mask);
    TempLocs = vertcat(temp.SpinePks);
    TempLocs = TempLocs(~TempLocs(:,6),3); % Time of firing that is not an AP
    Locs = [Locs; TempLocs(TempLocs>30)];
    LocLabels = [LocLabels; ones(numel(TempLocs(TempLocs>30)),1)*i];
    PkCount = [PkCount; Pharma.SpinepkStats.MultiRespondingFraction(mask,2)*100];
    CountLabels = [CountLabels; ones(sum(mask),1)*i];
    

figure('Theme','light'); set(gcf, 'Units', 'inches', 'Position', [1, 1, 6, 2]);
    histogram(TempLocs,24,'Normalization','countdensity',...
             'FaceColor', Colors(i,:), 'FaceAlpha', 1), %hold on
    ylim([0 15]), xlim([0 120]), 
    xlabel('Time (sec)'), ylabel('Counts/5 Sec')
    xline(artifactpts-1)
    set(gca, "FontSize", 20)
    legend(PharmaGroups{i},"Location", "best",...
           "Color",'none','EdgeColor','none'); 
end

for i = 1:max(CountLabels)
    neurons = CountLabels==i;
    tempstruct(i) = datastatsJake(PkCount(neurons));
    isnormal(i) = adtest(PkCount(neurons));
end,    PRespFraction = struct2table(tempstruct); clear tempstruct
PRespFraction.Labels = PharmaGroups';
disp(PRespFraction),

figure('Theme','light');
    set(gcf, 'Units', 'inches', 'Position', [1, 1, 7, 6]);
    boxchart(CountLabels, PkCount,'GroupByColor',CountLabels,'ColorGroupLayout','overlaid', ...
             'BoxEdgeColor','k','MarkerSize',10, 'MarkerStyle',...
             '.','MarkerColor','k','BoxFaceAlpha', 1)
    colororder(Colors)
    ylabel("% DS Responding/Neuron"), ylim([0 119])
    set(gca, 'XTick', 1:i, 'XTickLabels',PharmaGroups)
    set(gca, 'YTick', 0:20:100)%, 'XTickLabels',PharmaGroups)
    ax = gca; ax.FontSize = 24; ax.TickDir = 'none';

    [~,~,stats] = kruskalwallis(PkCount, CountLabels,'off');
    c = multcompare(stats,'Display','off',CriticalValueType='bonferroni');
    pairs = num2cell(c(:, 1:2),2);    mask = c(:,end)<0.05;
    sigstar(pairs(mask), c(mask,6)); hold off; % Add significance markers

    p = orderpvalues(c);
    plotheatmap(PharmaGroups, p);

%% Energy Calculations
[Area,~] = SpotSizeCalc(0.22,1.33,1.33,400,200,0);
PTAmps = [ 28 28 30 28 27 28 28 28 28 27.5 30 30 30 30 30];
z = 0.02; %cm
lambda = [1860 1880];  alpha = [14.19 31.08]; % Hale & Querry 1973 (1/cm)
alpha = linterp(lambda, alpha, 1875); 

% Pulse Train
%m = (121.9-88.7)/(37-28); b = 121.9-m*37;
X = [28 37];
Y = [88.7 121.9];
Period = 0.01;
%PT_Energy = mean(linterp(X,Y,PTAmps)*Period*100);
PT_Energy = 3.52*.25*100 % 7/10 CW measurement I=28.7A
PT_H = PT_Energy*exp(-alpha*z)*1e5/Area;

% Single Pulse
SPAmps = [ 25 25 25 25 24 24 24 23.5 25 25 23 24 24 23 23.5];
X = [23 35];
Y = [233.9 256.7];
Period = 0.1; %seconds
SP_Energy = mean(linterp(X,Y,SPAmps)*Period)
SP_Energy = 2.92*8; % 7/10 CW measurement I=24.2A
SP_H = SP_Energy*exp(-alpha*z)*1e5/Area;


rho = .9932; % density of water at 37C
cp = 4.18; % isobaric specific heat kJ/(kg*K) from engineering toolbox
mean((alpha*SP_H*2)./(rho*cp))+30 % Temperature conversions


% Pulse train energy at 37C: 91.28 +/- 4.17 mJ
% Single 8 ms pulse energy at 37C: 27.02 +/- 0.83 mJ
% Pulse train energy at ~18C: 121.9 mJ

%% Functions
function [ExperimentType] = AppendData(ExperimentType)

    ExperimentType.SingleAppened.Soma = []; 
    ExperimentType.SingleAppened.Dendrites = []; 
    ExperimentType.SingleAppened.Spines = [];

    ExperimentType.MultiAppened.Soma = []; 
    ExperimentType.MultiAppened.Dendrites = []; 
    ExperimentType.MultiAppened.Spines = [];

    temp = 0; DendriteParent = []; SpineNeuronParent = [];SpineParent =[];
    for N = 1:length(ExperimentType.Single)
        time = 752;
        ExperimentType.SingleAppened.Soma = [ExperimentType.SingleAppened.Soma ExperimentType.Single(N).Soma(1:time,1)];
        ExperimentType.SingleAppened.Dendrites = [ExperimentType.SingleAppened.Dendrites ExperimentType.Single(N).Dendrites(1:time,:)];
        ExperimentType.SingleAppened.Spines = [ExperimentType.SingleAppened.Spines ExperimentType.Single(N).Spines(1:time,:)];


        time = 900;
        ExperimentType.MultiAppened.Soma = [ExperimentType.MultiAppened.Soma ExperimentType.Multi(N).Soma(1:time,1)];
        ExperimentType.MultiAppened.Dendrites = [ExperimentType.MultiAppened.Dendrites ExperimentType.Multi(N).Dendrites(1:time,:)];
        ExperimentType.MultiAppened.Spines = [ExperimentType.MultiAppened.Spines ExperimentType.Multi(N).Spines(1:time,:) ];


        sz = size(ExperimentType.Single(N).Dendrites,2);
        DendriteParent = [DendriteParent; ones(sz,1)*N];
        sz = size(ExperimentType.Single(N).SpineParent,2);
        SpineNeuronParent = [SpineNeuronParent; ones(sz,1)*N];
        SpineParent =  [SpineParent ExperimentType.Single(N).SpineParent+temp];
        temp = max(SpineParent);
    end
    
    NeuronID = [(1:N)'; DendriteParent; SpineNeuronParent];
    DendriteID = [nan(N,1); (1:length(DendriteParent))'; SpineParent'];
    Structure = categorical([repelem("Soma", N),...
                repelem("Dendrite", length(DendriteParent)),...
                repelem("Spine", length(SpineParent))])';

    ExperimentType.Labels = table(NeuronID, DendriteID, Structure);
end

function [Traces, trough] = RemoveDC(Traces,Frame)
    smoothed = smoothdata(Traces,1,"movmean",30);
    trough = min(smoothed(Frame-100:Frame,:),[],1);
    Traces = (Traces-trough)./trough;
end

function ExpTyp = ExtractPks(ExpTyp,Frames,artifactpts)

    for N = 1:length(ExpTyp)
        
    % Soma
        Soma =  ExpTyp(N).Soma(1:Frames);              
        Fs = length(Soma)/ExpTyp(N).ndinfo.duration; % sampling Fq
        time = (1:length(Soma))/Fs;
        frame = round(30*Fs);
        
        Soma  = RemoveDC(Soma,frame); 
        [Somalocs, ~, height, width] = FindDerivativePks(Soma,artifactpts, Fs, 1);
    
        ExpTyp(N).SomaPks = [Somalocs, height, width];       
          

    % Spines 
        Spines  =  ExpTyp(N).Spines(1:Frames,:);  
        Spines  = RemoveDC(Spines,frame);
        [spinelocs, spinelabel, pkheight, FWHM, bAPlabel] ...
            = FindDerivativePks(Spines,artifactpts, Fs, 3, Somalocs);
        SomaLabel  = ones(numel(spinelabel),1)*N;
        SpineParent =  ExpTyp(N).SpineParent(spinelabel)';
    
        UniqueSpines = unique(spinelabel.*~bAPlabel);
        ActiveSpine = ismember(spinelabel, UniqueSpines);
        SpinePks = [spinelabel, SpineParent, spinelocs, pkheight, FWHM, bAPlabel, ActiveSpine, SomaLabel]; 
    
        SpineEventMask = zeros(size(Spines,2),1); 
        SpineEventMask(UniqueSpines(2:end)) = 1; 

        ExpTyp(N).SpinePks = SpinePks; 
        ExpTyp(N).SpineEventMask = SpineEventMask';
        ExpTyp(N).NumberSpineEvent = sum(~bAPlabel);% number of independent spine events
   

    end
end

function [SPkS] = SpinePeakStats(ExpTyp)

    timethresh = 30; % Seconds

% Single
if ~isempty(ExpTyp.Single)
    PksSpine = vertcat(ExpTyp.Single.SpinePks);
    ICS = ~logical(PksSpine(:,end-2)); % not bAP events
    ActiveSpines = logical(PksSpine(:,end-1)); % Active Spines only
    timemask = PksSpine(:,3)>timethresh; % After Stimulation

    for N = 1:max(PksSpine(:,end))
        NeuronMask = PksSpine(:,end)==N;
        TotalSpines = size(ExpTyp.Single(N).Spines,2);

        postmask = ICS & ActiveSpines & NeuronMask & timemask;
        premask = ICS & ActiveSpines & NeuronMask & ~timemask;
    
        SpikeFreq(N, 1) = numel(PksSpine(premask,1))/30;
        SpikeFreq(N, 2) = numel(PksSpine(postmask,1))/70;

        RespondingFraction(N, 1) = numel(unique(PksSpine(premask,1)))/TotalSpines;
        RespondingFraction(N, 2) = numel(unique(PksSpine(postmask,1)))/TotalSpines;
    end
    SPkS.SingleSpikeFreq = SpikeFreq;
    SPkS.SingleRespondingFraction = RespondingFraction;
    clear  SpikeFreq RespondingFraction,

    % Peaks before the stimulus that are not AP
    mask = ICS & ~timemask;
    SPkS.Pre_Spine = PksSpine(mask,:);
    
    % Peaks after the stimulus that are not AP
    mask = ICS & timemask;
    SPkS.Post_Single = PksSpine(mask,:);
    SPkS.Single_bAP = PksSpine(~ICS,:);
end

% Multi
if ~isempty(ExpTyp.Multi)
    PksSpine = vertcat(ExpTyp.Multi.SpinePks);
    ICS = ~logical(PksSpine(:,end-2)); % not bAP events
    ActiveSpines = logical(PksSpine(:,end-1)); % Active Spines only
    timemask = PksSpine(:,3)>timethresh; % After Stimulation

    for N = 1:max(PksSpine(:,end))
        NeuronMask = PksSpine(:,end)==N;
        TotalSpines = size(ExpTyp.Multi(N).Spines,2);

        postmask = ICS & ActiveSpines & NeuronMask & timemask;
        premask = ICS & ActiveSpines & NeuronMask & ~timemask;
    
        SpikeFreq(N, 1) = numel(PksSpine(premask,1))/30;
        SpikeFreq(N, 2) = numel(PksSpine(postmask,1))/90;

        RespondingFraction(N, 1) = numel(unique(PksSpine(premask,1)))/TotalSpines;
        RespondingFraction(N, 2) = numel(unique(PksSpine(postmask,1)))/TotalSpines;
    end
    SPkS.MultiSpikeFreq = SpikeFreq;
    SPkS.MultiRespondingFraction = RespondingFraction;

    
    % Peaks after the stimulus that are not AP
    mask = ICS & timemask;
    SPkS.Post_Multi = PksSpine(mask,:);
    SPkS.Multi_bAP = PksSpine(~ICS,:);
end

end
