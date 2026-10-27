    clearvars;
    global maxNoTimePoints;
    global peakSearchWinLeft;
    global peakSearchWinRight;
    global colors;

    maxNoTimePoints = 15000;                        % maximum number of data points in a single chromatogram (1200 data points per min)
    peakSearchWinLeft = 300;                        % search for peaks only up to 300 datapoints (15 sec in our dataset) to the left of the center of the recorded data for each peak (set to 0 for full window search)  
    peakSearchWinRight = 300;                       % search for peaks only up to 300 datapoints (15 sec in our dataset) to the right of the center of the recorded data for each peak (set to 0 for full window search)
    maxDelay = 100;                                 % max +-100 datapoints delay (+-5 sec ) between samples - used for data alignment of TIC:s
    
    massTolerance = 0.2;                            % tolerance in m/z between recorded and observed mass
    plotGroupFile = "Plotgroups.xlsx";              % file that describes what and how to plot data
    compTransFile = "CompoundTransitions.xlsx";     % file that describes MS-transitions to search for in the data
    outfile = "Chromatograms.pdf";                  % pdf-file where the integrated MSMS-chromatograms are saved

    colors = [0 0.45 0.75
              0.85 0.35 0.1
              0.9 0.7 0.1
              0.5 0.2 0.5
              0.5 0.7 0.2
              0.35 0.8 0.9
              0.65 0.1 0.2
              0.2 0.9 0.4
              0.5 0.5 0.5
              1.0 0.6 0.0];
    
   
    % Read plot and transition information from .xlsx-files
    plotTable = readtable(plotGroupFile, 'Sheet','Plot');
    groupTable = readtable(plotGroupFile, 'Sheet','Groups');
    compTransTable = readtable(compTransFile);
 
   % peakTable = readtable("chosenPeaksSupr.xlsx", 'Sheet','Peaks');
 %   peakName = table2array(peakTable(1:size(peakTable,1),"Compound"));
  %  peakNumber = table2array(peakTable(1:size(peakTable,1),"Peak"));
    
    files = table2array(plotTable(1:size(plotTable,1),"File"));
    compGroup = table2array(plotTable(1:size(plotTable,1),"CompoundGroup"));
    spectraColor = table2array(plotTable(1:size(plotTable,1),"Color"));
    shortDescr = table2array(plotTable(1:size(plotTable,1),"ShortDescription"));
    isPrimary = table2array(plotTable(1:size(plotTable,1),"PrimaryTransition"));
    lineStyle = table2array(plotTable(1:size(plotTable,1),"LineStyle"));
    plotScale = table2array(plotTable(1:size(plotTable,1),"PlotScale"));
    nameGroup = table2array(groupTable);

    noFiles = size(files, 1);
    dataFile = string(files);
    fig = figure();
    
    for i = 1:noFiles
        tmp = convertStringsToChars(dataFile(i));
        st = length(tmp) - 5;
        if st < 1
            st = 1;
        end
        
        if ~contains(tmp(st:end), '.')
            dataFile(i) = strcat(dataFile(i), ".mzML");
        end
        
        if ~strcmp(lineStyle(i), '--') && ~strcmp(lineStyle(i), ':') && ~strcmp(lineStyle(i), '-.') 
            lineStyle(i) = {'-'};
        end
        
        if class(plotScale) == "cell"
            if strcmp(plotScale(i), '')
                plotScale(i) = {'1'};
            end
        elseif isnan(plotScale(i))
            plotScale(i) = 1;
        end
    end
    
    if class(plotScale) == "cell"
        plotScale = str2double(plotScale);
    end

    if class(shortDescr) ~= "cell"
        tmp = isnan(shortDescr); 
        shortDescr = num2cell(shortDescr);
        shortDescr(tmp) = {' '};
    end

    compNr = table2array(compTransTable(1:size(compTransTable,1),"CompoundNr"));
    names = table2array(compTransTable(1:size(compTransTable,1),"Compound"));
    precursorMass = table2array(compTransTable(1:size(compTransTable,1),"Precursor"));
    fragmentMass = table2array(compTransTable(1:size(compTransTable,1),"Product"));
    primary = table2array(compTransTable(1:size(compTransTable,1),"Primary"));
    rt = table2array(compTransTable(1:size(compTransTable,1),"RT"));
    cEn = table2array(compTransTable(1:size(compTransTable,1),"CE"));
   
    % import data from the mzML-files specified in "Plotgroups.xlsx"
    
    clearvars spectraData;
    spectraData(noFiles) = spectraFile;
    noCompGroups = max(compGroup);
    averageTIC = zeros(maxNoTimePoints, noCompGroups);
    
    for n = 1:noFiles
        spectraData(n) = readmzML(dataFile(n));
   
        averageTIC(1:spectraData(n).TICsize, compGroup(n)) = averageTIC(1:spectraData(n).TICsize, compGroup(n)) + spectraData(n).TIC(1:spectraData(n).TICsize);
        spectraData(n).shortDescr = shortDescr(n);
    end

    compTIC = zeros(noCompGroups, 1);
    
    % Find the TIC-file in each compound group that is most similar to average TIC in that group 
    for n = 1:noCompGroups
    
        if sum(compGroup == n) == 0
            continue;
        end
        
        averageTIC(:, n) = averageTIC(:, n) / sum(compGroup == n);
        sqSum = sum(([spectraData(compGroup == n).TIC] - averageTIC(:, n)).^2, 1);
        [~, compTIC(n)] = min(sqSum);               % select TIC most similar to average TIC to align against
         
        comp = find(compGroup == n);
        
        % Find alignments for all chromatograms in respective compound group
        for f = 1:noFiles
            if compGroup(f) ~= n
                continue;
            end
    
            spectraData(f).delay = align(spectraData(f).TIC(1:min(spectraData(f).TICsize, spectraData(comp(compTIC(n))).TICsize)), spectraData(comp(compTIC(n))).TIC(1:min(spectraData(f).TICsize, spectraData(comp(compTIC(n))).TICsize)), maxDelay);
            
            % align data
            if spectraData(f).delay > 0
                for c = 1:spectraData(f).noChromatograms
                    spectraData(f).intArr(c, 1:end-spectraData(f).delay) = spectraData(f).intArr(c, spectraData(f).delay+1:end);
                    spectraData(f).intArr(c, end-spectraData(f).delay+1:end) = 0;
                end
            elseif spectraData(f).delay < 0
                for c = 1:spectraData(f).noChromatograms
                    spectraData(f).intArr(c, 1-spectraData(f).delay:end) = spectraData(f).intArr(c, 1:end+spectraData(f).delay);
                    spectraData(f).intArr(c, 1:-spectraData(f).delay) = 0;
                end
            end
        end
    end
    
    % match each transition found in the files to a compound based on the information in "CompoundTransitions.xlsx"
  
    for n = 1:noFiles
          
        tr = find(primary == isPrimary(n));

        for i = 1:spectraData(n).noChromatograms
          
            pm = find(abs(precursorMass - spectraData(n).precIon(i)) < massTolerance);
            fm = find(abs(fragmentMass - spectraData(n).prodIon(i)) < massTolerance);
            c = find(cEn == spectraData(n).CE(i));
            
            comp = pm(ismember(pm, fm) & ismember(pm, c) & ismember(pm, tr));
            
            if isempty(comp)
                continue;
            end
            
            c = 1;
            while (c <= size(comp, 1))      % Ignore compounds not part of the selected compound group 
                if sum(strcmp(names(comp(c)), nameGroup(:, compGroup(n)))) == 0   
                      comp(c) = 0;
                end
                c = c + 1;
            end
            
            comp = comp(comp ~= 0);
            if isempty(comp)
                continue;
            end

            if size(comp, 1) > 1            % If there is more than one compound with the same transition, the closest in retention time is assumed to be the right compound
                [~, closestComp] = min(abs(rt(comp) - spectraData(n).estRt(i)));
                comp = comp(closestComp);
            end
            
            % Check if the compound already have an associated recorded transition
            a = find(compNr(comp) == spectraData(n).compNr);
            b = find(names(comp) == spectraData(n).name);
            c = find(primary(comp) == spectraData(n).primary);
            
            tmp = a(ismember(a, b) & ismember(a, c));

            if isempty(tmp)
                spectraData(n).compNr(i) = compNr(comp);
                spectraData(n).name(i) = names(comp);         
                spectraData(n).primary(i) = primary(comp);
            else
                rt1 = abs(rt(comp) - spectraData(n).estRt(tmp));
                rt2 = abs(rt(comp) - spectraData(n).estRt(i));
                
                if rt2 < rt1
                    spectraData(n).compNr(i) = compNr(comp);
                    spectraData(n).name(i) = names(comp);         
                    spectraData(n).primary(i) = primary(comp);
                    
                    spectraData(n).compNr(tmp) = 0;
                    spectraData(n).name(tmp) = "";         
                end
                    
            end
            
        end
   
    end
  
    % Integrate, plot and save data according to specifications in the file specified in plotGroupFile
    
    fileID = fopen("Integrals.txt", "w");
    fprintf(fileID, "Filename;Compound;Data description;Retention time1 (min);Peak1 area;Retention time2 (min);Peak2 area;Retention time3 (min);Peak3 area\n");
  
    delete(outfile);
    
    uniqueNames = unique(names, 'stable');
    
    for i = 1:size(uniqueNames, 1)
       
        dataFound = plotAndSave(uniqueNames(i), spectraData, fileID, spectraColor, lineStyle, plotScale, string(files));
        
        if dataFound ~= 0
            exportgraphics(fig, outfile, 'Append', true);     % export chromatograms as vector images
 
            %exportgraphics(fig, strcat("im", outfile), 'Append', true, 'ContentType', 'image');   % alternatively export as bitmap to save space
        end
        
    end

    close(fig);

    fclose(fileID);

    return;