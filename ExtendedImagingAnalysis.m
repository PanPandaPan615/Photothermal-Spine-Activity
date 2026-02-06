clc, clear, close all

load("ExtendImaging.mat"), Extended = ExpData; 
ActinData = load("Actin_GCaMP.mat");
ActinData = ActinData.ExpData; clear ExpData

Extended = ExtractPks(Extended);
PkStats = vertcat(Extended.SpinePks);
bAP = logical(PkStats(:,end));
SpineEvent = PkStats(~bAP,1:end-1);

artifactpts = (30:5:50)+1;
TimeLabels = [0,3,10:10:80];

for tp = 1:max(SpineEvent(:,2))
    tmask = SpineEvent(:,2)==tp;
    Locs = SpineEvent(tmask,5);
    for N = 1:max(SpineEvent(:,1))
        mask = SpineEvent(:,1)==N & tmask;
        spines = SpineEvent(mask,3);
        count(tp,N) = numel(unique(spines))/numel(Extended(N).Spines(1,:,1))*100;
    end
    isnormal(tp) = adtest(count(tp,:));
    tempstruct(tp) = datastatsJake(count(tp,:)');
    colors(tp,:) = [ 1 tp/12 0];
    % % Time Histograms
    % figure('Theme','light'); set(gcf, 'Units', 'inches', 'Position', [1, 1, 6, 2]);
    %     histogram(Locs,24,'Normalization','countdensity',...
    %              'FaceColor', colors(tp,:), 'FaceAlpha', 1), %hold on
    %     ylim([0 15]), xlim([0 120]), 
    %     xlabel('Time (sec)'), ylabel('ICS/5 Sec')
    %     set(gca, "FontSize", 20)
    %     time = num2str(TimeLabels(tp));
    %     legend([time,' min'],"Location", "northwest",...
    %           "Color",'none','EdgeColor','none'); 
end
RespFraction = struct2table(tempstruct);% clear tempstruct
RespFraction.Labels = num2cell(TimeLabels)';
disp(RespFraction)
abs(RespFraction.mean-RespFraction.median)

TimeLabels = repmat(TimeLabels,[size(count,2) 1])';
[pkw,~,stats] = kruskalwallis(count(:), TimeLabels(:),'off');
c = multcompare(stats,'Display','off',CriticalValueType='dunnett');

figure(Theme='light'), set(gcf, 'Units', 'inches', 'Position', [1, 1, 13, 5]);
    boxchart(TimeLabels(:),count(:),'GroupByColor',TimeLabels(:),'ColorGroupLayout','overlaid', ...
                 'BoxEdgeColor','k','MarkerSize',10, 'MarkerStyle',...
                 '.','MarkerColor','k','BoxFaceAlpha', 1, BoxWidth=1.5); 
    colororder(colors)    
    xlim([-1 85]), ylim([0 100]), 
    ylabel("% Spines Responding"), xlabel("Time (min)")
    

    pairs = [TimeLabels(c(:, 1),1) TimeLabels(c(:, 2),1)];
    pairs = num2cell(pairs,2);    mask = c(:,6)<0.05;
    sigstar(pairs(mask), c(mask,6)); hold off
    set(gca, "FontSize", 24)

%% 
close all
spines = vertcat(ActinData.Spines);
mask = ~sum(spines==0,2)>0;
spines = spines(mask,:);

normalized_spines = (spines(:,2:end)-spines(:,1))./spines(:,1);
TimeLabels = repmat(1:9,[length(spines),1]);
SpineVect = normalized_spines(:);

ExtIdxs = [4:11]; 
count = []; SpineCount = []; NeuronLabels = [];
for i = 1:length(ActinData)
    spines = ActinData(i).Spines;
    spines = (spines(:,2:end)-spines(:,1))./spines(:,1);
    PkStats = Extended(ExtIdxs(i)).SpinePks;
    bAP = logical(PkStats(:,end));
    SpineEvent = PkStats(~bAP,2:end-5);
    for j = 1:length(spines)
        SpineCount(j) = sum(SpineEvent(:,2) == j);
    end
    count = [count, SpineCount]; clear SpineCount

    NeuronLabels = [NeuronLabels; ones(size(spines))*i];

end
count = count(mask);
NeuronLabels = NeuronLabels(mask,:);
copy = count;
count(count>prctile(count,95)) = floor(prctile(count,95));
count = repmat(count,[9 1])';


[pkw,~,stats] =kruskalwallis(SpineVect,TimeLabels(:),"off");
c = multcompare(stats,'CriticalValueType','bonferroni','Display','off');
p = orderpvalues(c);
Labels = [1,10:10:80];

plotheatmap(num2cell(Labels), p),

count = copy';

% Activity Masks
ActivitlLab = {'No Activity','Low Activity','High Activity'};
thresh = median(count)+mad(count);
HA = count>thresh;%count>prctile(count,75);
LA = ~HA & count>0;%count>prctile(count,50) & count<=prctile(count,75);
NA = count==0;%count>prctile(count,25) & count<=prctile(count,50);
ActivityMask = cat(2,NA,LA,HA);

colors = colororder('sail');

% Box chart separating mCardinal Fluorescence by spiking level
close all
ActivityLabels = ActivityMask*(1:3)';
figure(Theme='light');   set(gcf, 'Units', 'inches', 'Position', [1, 1, 13, 5]);
    boxchart(TimeLabels(:),SpineVect,'GroupByColor',ActivityLabels(:), ...
             'BoxEdgeColor','k','MarkerSize',10, 'MarkerStyle',...
             '.','MarkerColor','k','BoxFaceAlpha', 1, BoxWidth=0.5);  hold on


    xlim([0.5 9.5]); colororder('sail')
    ylabel('mCardinal \DeltaF/F'),  xlabel('Time (min)'), 
    set(gca, 'FontSize', 24, 'XTickLabel', [1,10:10:80]);   
    legend([ActivitlLab{:}, 'Sham'], "Color",'none','EdgeColor','none',...
      'Orientation','horizontal','Location','northwest')  

for tp = 1:max(TimeLabels(:))
    mask = TimeLabels==tp & ActivityLabels==1;
    NAstruct(tp) = datastatsJake(normalized_spines(mask));
    mask = TimeLabels==tp & ActivityLabels==2;
    LAstruct(tp) = datastatsJake(normalized_spines(mask));
    mask = TimeLabels==tp & ActivityLabels==3;
    HAstruct(tp) = datastatsJake(normalized_spines(mask));
end

NAFraction = struct2table(NAstruct); NAFraction.Labels = num2cell(Labels)';
LAFraction = struct2table(LAstruct); LAFraction.Labels = num2cell(Labels)';
HAFraction = struct2table(HAstruct); HAFraction.Labels = num2cell(Labels)';

disp(NAFraction)
%% P-Vals per timepoint
close all

for i = 1:size(normalized_spines,2)
    [p_time(i),tbl,stats] = kruskalwallis(normalized_spines(:,i), ActivityLabels(:,1), 'off');
    c = multcompare(stats,'CriticalValueType','bonferroni','Display','off');
    p = orderpvalues(c);   
    plotheatmap({'NA','LA','HA'}, p);

    [rho(i),pval(i)] = corr(normalized_spines(:,i),count,'Type','Spearman');

    H = tbl{2,5};       % Kruskal-Wallis H (chi-square)
    k = numel(stats.gnames);    % number of groups
    n = numel(normalized_spines(:,i));  % total sample size
    
    % Epsilon-squared
    epsilon2(i) = (H - k + 1) / (n - k);
    
    % Eta-squared (H-based)
    eta2_H(i) = H / (n - 1);
end

% P-vals per phenotype
for i = 1:3
    mask = ActivityLabels==i;
    [~,~,stats] = kruskalwallis(normalized_spines(mask), TimeLabels(mask), 'off');
    c = multcompare(stats,'CriticalValueType','bonferroni','Display','off');
    p = orderpvalues(c);   
    plotheatmap(num2cell(Labels), p);
end

%% Distribution per neurons
close all
figure(Theme='light');
    histogram(count,BinMethod="fd",FaceAlpha=1)
    xline(1,LineWidth=2,Alpha=1)
    xline(thresh,LineWidth=2,LineStyle="--",Alpha=1)
    legend('','1','Median+MAD', "Color",'none','EdgeColor','none')
     xlabel('# ICS Observed'), ylabel('Dendritic Spines'),
     ax = gca; ax.FontSize = 24;


figure(Theme='light');
    for i = 1:max(NeuronLabels(:))
        temp = NeuronLabels(:,1)==i;
        NDist(i,:) = ([sum(temp&NA) sum(temp&LA) sum(temp&HA)]./sum(temp))*100;
        
        for j = 1:9
            [indyrho(i,j), indypval(i,j)] =corr(normalized_spines(temp,j), count(temp),'Type','Spearman');
        end
    end
    bar(NDist, 'stacked'); colororder('sail')
    set(gca, 'FontSize', 24)  
    ylabel('Fraction Spines'), xlabel('Neuron ID')
    legend(ActivitlLab, "Color",'none','EdgeColor','none',...
          'Orientation','horizontal','Location','northoutside')  
    
    indyrho = indyrho.* (indypval<0.1)
%% Individual Neurons separated by reaction type 
close all, clear p
NeuronID = concatCellVectors({'Neuron'}, {1:max(NeuronLabels(:))}');

figNA = figure(Theme='light');
ylabel('mCardinal \DeltaF/F'),  xlabel('Time (min)'), 
title(ActivitlLab{1}), set(gca, 'FontSize', 24), ylim([-1, .1]), xlim([0 85])
figLA = figure(Theme='light');
ylabel('mCardinal \DeltaF/F'),  xlabel('Time (min)'), 
title(ActivitlLab{2}), set(gca, 'FontSize', 24) , ylim([-1, .1]), xlim([0 85])
figHA = figure(Theme='light');
ylabel('mCardinal \DeltaF/F'),  xlabel('Time (min)'), 
title(ActivitlLab{3}), set(gca, 'FontSize', 24) , ylim([-1, .1]), xlim([0 85])
space = .5;
for i = 1:max(NeuronLabels(:))
    nmask = NeuronLabels == i;
    for j =1:max(ActivityLabels(:))
        mask = nmask & ActivityLabels == j;
        DS = normalized_spines(mask);
        DS = reshape(DS, [length(DS)/9 9])';
        
        MED(:,j) = mean(DS,2); %median(DS,2)
        n = sum(mask);
        QTLS(:,:,j) = [-std(DS,[],2), std(DS,[],2)]./sqrt(n(1)); %[prctile(DS,25,2), prctile(DS,75,2)]
        %QTLS(:,:,j) = QTLS(:,:,j) + MED(:,j);
    end
    t = Labels;%+space*i;
    figure(figNA), hold on
    errorbar(t,MED(:,1),QTLS(:,1,1),QTLS(:,2,1) ,'LineWidth',1.5),
    figure(figLA), hold on
    errorbar(t,MED(:,2),QTLS(:,1,2),QTLS(:,2,2) ,'LineWidth',1.5),
    figure(figHA), hold on
    errorbar(t,MED(:,3),QTLS(:,1,3),QTLS(:,2,3) ,'LineWidth',1.5),
end
figure(figNA),
legend( NeuronID{:}, "Color",'none',...
        'EdgeColor','none','Location','eastoutside') 
figure(figLA),
legend( NeuronID{:}, "Color",'none',...
        'EdgeColor','none','Location','eastoutside') 
figure(figHA),
legend( NeuronID{:}, "Color",'none',...
        'EdgeColor','none','Location','eastoutside') 


figure(Theme='light');
alpha = 0.05/max(ActivityLabels(:)); 
for i = 1:max(NeuronLabels(:))
    nmask = NeuronLabels(:,1) == i;
    for j = 1:length(Labels)
    p(j) = kruskalwallis(normalized_spines(nmask,j),ActivityLabels(nmask,1),"off");
    %c = multcompare(stats,Display='off',CriticalValueType='dunn-sidak');
    end
    plot(Labels,p,"LineWidth",1), hold on
    pass(:,i) = p<alpha;
end
yline(alpha,LineWidth=2,LineStyle="--",Alpha=1)
set(gca, 'FontSize', 24) 
ylabel('P-value'),  xlabel('Time (min)'), 
legend( [NeuronID{:},"\alpha = 0.017"], "Color",'none',...
        'EdgeColor','none','Location','eastoutside') 

figure(Theme='light');
    bar(sum(pass,1)./size(pass,2)*100), ylim([0 100]) 
    set(gca, 'XTickLabels', Labels);
    xlabel('Time (min)'),  ylabel('Fraction of Neurons Passing Test (%)')
    set(gca, 'FontSize', 24) 
    
Rho = corr(normalized_spines,count,"Type","Spearman");

%% Functions
function ExpTyp = ExtractPks(ExpTyp)
artifactpts = (30:5:50)+1;

    for N = 1:length(ExpTyp) % for each neuron
% Init  
        Soma  =  ExpTyp(N).Soma;
        Spines  =  ExpTyp(N).Spines;  
        Fs = length(Soma)/ExpTyp(N).ndinfo.duration; % sampling Fq
        frame = round(30*Fs);
        time = (1:length(Soma))/Fs;
        SomaPks = []; SpinePks = [];
        

        %plot(diffs), hold on 
        %yline(SomaThresh), hold off
    
        for tp = 1:size(Soma,2) % for each timepoint
% Soma
            SomaIter  = Soma(:,tp); 
            SomaIter  = RemoveDC(SomaIter,frame);
            if tp == 2
              [locs, ~, height, width] = FindDerivativePks(SomaIter,artifactpts, Fs, 1);
            else
              [locs, ~, height, width] = FindDerivativePks(SomaIter,0, Fs, 1);
            end

            SomaPks = [SomaPks; repmat(tp, size(locs)), locs, height, width];       
            SomaPkLocs =  locs;

            %plot(time, SomaIter), hold on, scatter(locs, height,"filled"), hold off
            %set(gcf, 'Units', 'inches', 'Position', [1, 1, 10, 5]);            tp
          
% Spines - for each neuron load spines into  temporary variables  
            SpinesIter = Spines(:,:,tp); 
            SpinesIter  = RemoveDC(SpinesIter,frame);

            if tp == 2
                [spinelocs, spinelabel, pkheight, FWHM, bAPlabel] ...
                = FindDerivativePks(SpinesIter,artifactpts, Fs, 3, SomaPkLocs);            
            else
                [spinelocs, spinelabel, pkheight, FWHM, bAPlabel] ...
                 = FindDerivativePks(SpinesIter,0, Fs, 3, SomaPkLocs);
            end

            SomaLabel  = ones(numel(spinelabel),1)*N;
            SpineParent =  ExpTyp(N).SpineParent(spinelabel)';
        
            UniqueSpines = unique(spinelabel.*~bAPlabel);
            spinepkstats = [SomaLabel, repmat(tp, size(spinelocs)),...
            spinelabel, SpineParent, spinelocs, pkheight, FWHM, bAPlabel]; 

            SpinePks = [SpinePks; spinepkstats]; 
        end
        mask = SpinePks(:,end-1) ~= 0;
        ExpTyp(N).SpinePks = SpinePks(mask,:);

        ExpTyp(N).SomaPks = SomaPks;
    end
end

function [Traces, trough] = RemoveDC(Traces,Frame)
    smoothed = smoothdata(Traces,1,"movmean",30);
    trough = min(smoothed(Frame-100:Frame,:),[],1);
    Traces = (Traces-trough)./trough;
end

function result = concatCellVectors(vec1, vec2)
    % Convert all elements to strings (if they are numeric)
    strVec1 = cellfun(@string, vec1, 'UniformOutput', false);
    strVec2 = cellfun(@string, vec2, 'UniformOutput', false);

    % Generate all combinations using ndgrid
    [A, B] = ndgrid(1:numel(strVec1), 1:numel(strVec2));

    % Concatenate strings with a space between
    result = arrayfun(@(i, j) strcat(strVec1{i}, " ", strVec2{j}), A(:), B(:), 'UniformOutput', false);
end
