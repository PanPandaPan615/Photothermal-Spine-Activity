function [pklocs, spinelabel, pkheight, FWHM, bAPlabel]...
          = FindDerivativePks(trace, artifacts, Fs, threshold, SomaLocs)

    pklocs = []; spinelabel = []; bAPlabel =[]; 
    artifacts = round(artifacts*Fs);
    
    if nargin<4,  threshold = 3; end

    [pklocs, spinelabel, pkheight, FWHM] = locatepks(trace,artifacts,Fs,threshold);

    if ~isvector(trace)
        
        mtrace = mean(trace,2);
        bAPlocs = locatepks(mtrace,artifacts,Fs,1);
        if isempty(bAPlocs), bAPlocs=0; end
        if isempty(SomaLocs), SomaLocs=0; end
        bAPlocs = [bAPlocs; SomaLocs];

        if ~isempty(pklocs)
            bAPlabel = abs(repmat(pklocs',[length(bAPlocs) 1]) - bAPlocs)<1;
            bAPlabel = logical(sum(bAPlabel,1)');
        end

        mask = FWHM==0 | pkheight==0; % less than temporal resolution is noise
        FWHM(mask) = [];   pklocs(mask) = [];
        spinelabel(mask) = []; pkheight(mask) = [];
        bAPlabel(mask) = [];
    end

    % for i = 1:size(trace,2)
    %     plot(strace(:,i)),hold on
    %     plot(outlier(:,i)), hold off
    %     pause(1)
    %     i
    % end
   
end

function [locs, labels, pkheight, FWHM] = locatepks(trace,artifacts,Fs,threshold)
    L = length(trace);
    dtrace = [zeros(size(trace,2),1)'; diff(trace,1,1)];

    %simple detection
    threshold = std(dtrace,[],1)*threshold;
    outlier = isoutlier(dtrace,'grubbs')  & dtrace>threshold & islocalmax(dtrace); 
    outlier = movmax(outlier,[5 2],1);
    ftrace = outlier.*trace;
    outlier = islocalmax(ftrace,MinProminence=mean(threshold),ProminenceWindow=Fs);
    [locs, labels] = find(outlier);
    
    % %advanced detection
    % dtrace = abs(dtrace);
    % strace = movmedian(dtrace,7);
    % 
    % threshold = median(dtrace) + threshold.*mad(dtrace);
    % ftrace = strace>threshold;
    % ftrace = trace.*ftrace;
    % outlier = islocalmax(ftrace,MinProminence=3*mad(ftrace(:))); 
    % [pks, labels] = find(outlier);
    
    % Throw out potential half responses
    mask = ~isbetween(locs,Fs*2,L-Fs*2);
    locs(mask) = [];    labels(mask) = [];

    % Throw out artifacts
    if sum(artifacts ~=0)>1
        [locs, labels] = removeartifacts(locs, labels, artifacts);
    end
    % Get metrics
    [locs, labels, pkheight, FWHM] = getpkmetrics(trace, locs, labels);
    
    FWHM = FWHM/Fs;    locs = locs/Fs;
end

function [locs, labels, pkheight, FWHM] = getpkmetrics(trace, locs, labels)
    pkheight = [];     FWHM = [];

    tscale = 0.01;
    trace = trace - smoothdata(trace,1,"movmedian",35);
    trace(trace<0) = 0;

    n = unique(labels);
    for j = 1:numel(n)
        
        mask = labels == n(j);
        TempLocs = locs(mask);
        DS = trace(:,n(j));
        tempheight = DS(TempLocs);

        for i = 1:sum(mask)
            loc = TempLocs(i);
            window = DS( loc-14 : loc+14 );
            window = interp1(window,1:tscale:numel(window),'linear');
            if tempheight(i) < 0.01
                widths(i) = 0;
                continue
            end
            pk = find(window==tempheight(i),1,"first");
            HM = tempheight(i)/2;

            LH = abs(fliplr(window(1:pk))-HM);
            RH = abs(window(pk:end)-HM);
            
            t1 = find(islocalmin(RH) & RH <0.05,1,"first"); 
            t2 = find(islocalmin(LH) & LH <0.05,1,"first");
            widths(i) = (t1 + t2)*tscale;

             [~, minima] = mink(abs(window-HM),10);
            % minima = sort(minima);
            % widths(i) = (minima(end)-minima(1))*tscale;

        end

        pkheight = [pkheight; tempheight];
        FWHM = [FWHM; widths']; 
        clear widths
    end

end

function [locs, labels] = removeartifacts(locs, labels, artifacts)

    for i = 1:length(artifacts)
        mask = abs(locs-artifacts(i))<7;
        labels(mask) = [];   locs(mask) = [];
    end

end