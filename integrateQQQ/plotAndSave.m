function dataFound = plotAndSave(name, spectraData, fileID, spectraColor, lineStyle, plotScale, files) 
    
    global colors;
    global maxNoTimePoints;
    global peakSearchWinLeft;
    global peakSearchWinRight;
    
    noSpectra = size(spectraData, 2);
    minPeakHalfWidth = 15;       % at least 2*15 data points (1.5 s) peak half width 
    minProm = 5;                 % min peak prominence
    minHeight = 10;              % min peak height

    combinedIntData = zeros(maxNoTimePoints, 1);
    noChrom = 0;

    chromNr = zeros(noSpectra, 1);        % chromatogram number corresponding to compound "name" in each file
    
 
    % sum up selected compound spectra to establish common integration limits
    for i = 1:noSpectra

        comp = find(spectraData(i).name == name);

        if isempty(comp)
            continue;
        end

        if size(comp, 1) > 1        % if there is a primary transition, then only combine primary transitions 
            comp = comp(spectraData(i).primary(comp) == 1);
        end
        
        chromNr(i) = comp;

        combinedIntData = combinedIntData + plotScale(i) * transpose(spectraData(i).intArr(comp, :));

        noChrom = noChrom + 1;
    end
    
    if noChrom == 0
        dataFound = false;
        return;
    end
    
    combinedIntData = combinedIntData/noChrom;      % normalize combined against # chromatograms
    
    winStart = find(combinedIntData ~= 0, 1, 'first');
    winEnd = find(combinedIntData ~= 0, 1, 'last');
    
    if (peakSearchWinLeft ~= 0)
        winMid = floor((winEnd - winStart) / 2) + winStart;
        winStart = winMid - peakSearchWinLeft;
    end

    if (peakSearchWinRight ~= 0)
        winMid = floor((winEnd - winStart) / 2) + winStart;
        winEnd = winMid + peakSearchWinRight; 
    end
    
    % Find peaks in combined data
    [~, locs] = findpeaks(combinedIntData(winStart:winEnd),"MinPeakProminence", minProm, "MinPeakHeight", minHeight, "MinPeakDistance", minPeakHalfWidth);
    
    if ~locs
        return;
    end

    locs = locs + winStart;
    peakStartArr = zeros(size(locs, 1), 1);
    peakEndArr = zeros(size(locs, 1), 1);
    peakAreaArr = zeros(size(locs, 1), 1);
    
    noPeaks = 0;

    % Locate and save start and end points for peaks in combined data 
    for p = 1:size(locs, 1)
        peakStartArr(p) = size(combinedIntData, 1); 
        peakEndArr(p) = 1;
   
        for backw = locs(p):-1:2 
            
            if (combinedIntData(backw) - combinedIntData(backw - 1)) / combinedIntData(locs(p)) < 0.0001 ...
                && locs(p) - backw >= minPeakHalfWidth
      
                peakStartArr(p) = backw; 
                break;

            elseif backw == 2
                peakStartArr(p) = 1;
                break;
            end
                        
        end
        
        for forw = locs(p):maxNoTimePoints-1 

             if (combinedIntData(forw) - combinedIntData(forw + 1)) / combinedIntData(locs(p)) < 0.0001 ...
                && forw - locs(p) >= minPeakHalfWidth
                                                
                peakEndArr(p) = forw;
                break;

            elseif forw == maxNoTimePoints-1
             
                peakEndArr(p) = maxNoTimePoints;
                break;
             end

        end
                       
        % Save areas of combined data
        
        if peakStartArr(p) < peakEndArr(p)
            noPeaks = noPeaks + 1;
            
            peakAreaArr(p) = trapz(combinedIntData(peakStartArr(p):peakEndArr(p)));
        end
    end

    % Sort combined peaks from largest to smallest to select the largest peaks for individual integration

    [peakAreaArr, index] = sort(peakAreaArr, 'descend');
    peakStartArr = peakStartArr(index);
    peakEndArr = peakEndArr(index);
    locs = locs(index);
    
    if noPeaks > 1 && peakAreaArr(2) / peakAreaArr(1) < 0.01        % do not consider second peak if it is has an area less than 1% of the main peak
        noPeaks = 1;
    end

    noPeaks = min(noPeaks, 3);  % maximum 3 peaks are integrated
    
    % Integrate peaks of selected compound in each file
    
    for f = 1:noSpectra
        if chromNr(f) == 0      % selected compound, "name", is not present in file f
            continue;
        end

        for p = 1:noPeaks
                       
            ar = trapz(spectraData(f).intArr(chromNr(f), peakStartArr(p):peakEndArr(p))); 
            spectraData(f).area(chromNr(f), p) = ar * 0.05;     % convert area to counts * s
                        
        end
    end

    % Set up plot layout
    plotDescr = strings(noSpectra*2, 1);
    plArr = zeros(noChrom, 1);  
    clf(gcf);
    set(gcf,'position',[100,100,1000,300])
    c = 0;
    
    % Plot and save data
   
    for f = 1:noSpectra
        if chromNr(f) == 0      % selected compound, "name", is not present in file f
            continue;
        end

        c = c + 1;
        
        spectraData(f).intArr(chromNr(f), (spectraData(f).intArr(chromNr(f), :) == 0)) = NaN;
        spectraData(f).timeArr(chromNr(f), (spectraData(f).timeArr(chromNr(f), :) == 0)) = NaN;

        % Plot aligned, but not baseline corrected data for main product in overview plot
        subplot(1, 5, [1 2]);   
        hold('on');
        plArr(c) = plot(spectraData(f).timeArr(chromNr(f), :), spectraData(f).intArr(chromNr(f), :) + spectraData(f).baseLine(chromNr(f)));
        
        set(plArr(c), 'Color', colors(spectraColor(f), 1:3), 'lineStyle', lineStyle(f)); 
        
        plotDescr(c) = spectraData(f).shortDescr;
        
        if spectraData(f).primary(chromNr(f)) == 0
            plotDescr(c) = strcat(plotDescr(c), " 2:nd");
        end

        % Write compound information to file fileID     
        fprintf(fileID, "%s;%s;%s; ", files(f), spectraData(f).name(chromNr(f)), plotDescr(c));
           
        for p = 1:noPeaks

            % Save main product rt and area to file
            fprintf(fileID, "%.2f;%.0f;", locs(p)/1200, spectraData(f).area(chromNr(f), p));
          
            subplot(1, 5, 2 + p);   % integration window subplot
            hold('on');
              
            plot(spectraData(f).timeArr(chromNr(f), peakStartArr(p):peakEndArr(p)), plotScale(f) * spectraData(f).intArr(chromNr(f), peakStartArr(p):peakEndArr(p)), 'Color', colors(spectraColor(f), 1:3), 'LineStyle', lineStyle(f));
            
        end
          
        fprintf(fileID, "\n");
    end

    % Add legend to plot
    
    subplot(1, 5, [1 2]);
    
    % remove duplicated lables
    [plotDescr, ia, ~] = unique(plotDescr, 'stable');
    
    if plotDescr(end) == ""
        ia = ia(1:end-1);
        plotDescr = plotDescr(1:end-1);
    end
    
    plArr = plArr(ia);

    legend(plArr, plotDescr, 'location', 'eastoutside');
    
    n = find(chromNr ~= 0);
    n = n(1);

    t = strcat(name, " ", string(round(spectraData(n).precIon(chromNr(n)), 1)), "->", string(round(spectraData(n).prodIon(chromNr(n)), 1)), " @ ", string(spectraData(n).CE(chromNr(n))), "V");
    
    title(t);
  
    for p = 1:noPeaks
        subplot(1, 5, p + 2);
        t = strcat("Peak ", num2str(p));
        title(t); 
    end

    dataFound = true;

    return;
end