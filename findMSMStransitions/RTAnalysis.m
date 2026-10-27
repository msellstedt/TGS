% Uses MS-data from both data-sets to assign tentative retention times for all products and isomers
% The results are saved in tab-separated text format in RT.txt and as graphs in RTAnalysis.pdf

global productName;

minHeight = 4000;           % minimal peak height (in m/z counts) to be counted as a peak in dataset 1
minHeight2 = 4000;          % minimal peak height (in m/z counts) to be counted as a peak in dataset 2
minMSPeakSeparation = 0.02; % minimal peak rt separation (in min) for each dataset
rtTolerance = 0.1;          % allowed dist (in min) between peaks in the two datasets to be counted as the same
maxRTBetweenIsomers = 0.75; % maximal dist (in min) between peaks in each dataset to be counted as possible isomers
minRatioOfIsomers = 0.15;   % minimal peak height ratio between peaks in a dataset to be counted as possible isomers
outputPDF = "RTAnalysis.pdf";

fileID = fopen("RT.txt", 'w');

delete(outputPDF);

fig = figure(1);
set(fig,'position',[100,100,700,400]);
clf;

pause('on');

hold off;

peakDetIn2Sum = 0;

for i = 1:size(triazoleMSProd1, 1)
    for j = 1:size(triazoleMSProd1, 2)
  
        hold on;
 
        rt1 = [-30
               -30];
        rt2 = [-100
               -100];
 
        calcFoundString = string;
        
        if ~isempty(triazoleMSProd1(i, j).foundMass)
 
            plot(triazoleMSProd1(i,j).retentionTimes, triazoleMSProd1(i,j).intensity, 'Color', [0 0.45 0.75]);
                    
            locsMS = double.empty();
            [pks, locsMS] = findpeaks(triazoleMSProd1(i, j).intensity, triazoleMSProd1(i, j).retentionTimes, "MinPeakHeight", minHeight, "MinPeakDistance", minMSPeakSeparation); % find peaks in data set 1
            
            if ~isempty(locsMS) 
                % peaks > minPeakHeight found in data set 1

                for l = 1:size(locsMS, 1)  % print found peaks on graph
                    text(locsMS(l), pks(l), sprintf("%.2f", locsMS(l)), 'FontSize',10, 'Color', [0 0.45 0.75]);
                end                
                
                if size(locsMS, 1) < 2                                      % fill with dummy peak if no 2:nd peak detected
                    locsMS = cat(1, locsMS, -30);
                    pks = cat(1, pks, 0);
                end

                [highest, rtpos] = maxk(pks, 2);

                if size(highest, 1) > 1
                    if highest(2) / highest(1) < minRatioOfIsomers          % remove potential isomeric peaks < minRatioOfIsomers
                        locsMS(rtpos(2)) = -30;
                    end
                end

                rt1 = locsMS(rtpos);
                calcFoundString = strcat(triazoleMSProd1(i, j).name, " calcd. mass [M+H]+: ", num2str(triazoleMSProd1(i, j).theoreticalMass, '%.4f'), ", found: ", num2str(triazoleMSProd1(i, j).foundMass, '%.4f'), " in dataset1");      
            end
        end
        
        
        if ~isempty(triazoleMSProd2(i, j).foundMass) 
            % plot extracted MS-traces in dataset 2
            
            plot(triazoleMSProd2(i, j).retentionTimes, triazoleMSProd2(i, j).intensity, 'Color', [0.8 0.2 0]);
            
            locsMS = double.empty();
            [pks, locsMS] = findpeaks(triazoleMSProd2(i, j).intensity, triazoleMSProd2(i, j).retentionTimes, "MinPeakHeight", minHeight2, "MinPeakDistance", minMSPeakSeparation); % find peaks in data set 1
            
            if ~isempty(locsMS)
                % peaks found in data set 2

                for l = 1:size(locsMS, 1) % print found peaks on graph
                    text(locsMS(l), pks(l), sprintf("%.2f", locsMS(l)), 'FontSize', 10, 'Color', [0.8 0.2 0]);
                end
                
                if size(locsMS, 1) < 2
                    locsMS = cat(1, locsMS, -100);
                    pks = cat(1, pks, 0);
                end

                [highest, rtpos] = maxk(pks, 2);
                
                if size(highest, 1) > 1
                    if highest(2) / highest(1) < minRatioOfIsomers          % remove potential isomeric peaks < minRatioOfIsomers
                        locsMS(rtpos(2)) = -100;
                    end
                end
                
                rt2 = locsMS(rtpos);   
                if calcFoundString ~= ""
                    calcFoundString = strcat(calcFoundString, " & ", num2str(triazoleMSProd2(i, j).foundMass, '%.4f'), " in dataset2");      
                else
                    calcFoundString = strcat(triazoleMSProd2(i, j).name, " calcd. mass [M+H]+: ", num2str(triazoleMSProd2(i, j).theoreticalMass, '%.4f'), ", found: ", num2str(triazoleMSProd2(i, j).foundMass, '%.4f'), " in dataset2");      
                end
            end
        end  
        
        if calcFoundString ~= ""
            title(calcFoundString, ' ');
        end

        [memb, membpos] = ismembertol(rt1, rt2, rtTolerance, 'DataScale', 1);
        rt = zeros(2, 1);
        

        if memb(1)                                                  % save large peaks common in both data-sets
            rt(1) = 0.5 * (rt1(1) + rt2(membpos(1))); 
        end
                
        if memb(2)
 
            if memb(1)
                rt(2) = 0.5 * (rt1(2) + rt2(membpos(2)));
            else
                rt(1) = 0.5 * (rt1(2) + rt2(membpos(2)));
            end
        end


        if ~rt(1)
            if rt1(1) > 0
                rt(1) = rt1(1);
                if rt1(2) > 0
                    rt(2) = rt1(2);
                else
                    if rt2(1) > 0 && abs(rt(1) - rt2(1)) < maxRTBetweenIsomers
                        rt(2) = rt2(1);
                    elseif rt2(1) > 0 && abs(rt(1) - rt2(2)) < maxRTBetweenIsomers
                        rt(2) = rt2(2);
                    end
                end
            else 
                if rt1(2) > 0
                    rt(1) = rt1(2);
                    if rt2(1) > 0 && abs(rt(1) - rt2(1)) < maxRTBetweenIsomers
                        rt(2) = rt2(1);
                    elseif rt2(1) > 0 && abs(rt(1) - rt2(2)) < maxRTBetweenIsomers
                        rt(2) = rt2(2);
                    end
                else
                    if rt2(1) > 0
                        rt(1) = rt2(1);
                        if rt2(2) > 0
                            rt(2) = rt2(2);
                        end
                    end   
                end                        
            end
            
        elseif ~rt(2)
            if ~memb(1) && abs(rt(1) - rt1(1)) < maxRTBetweenIsomers
                rt(2) = rt1(1);
            elseif ~memb(2) && abs(rt(1) - rt1(2)) < maxRTBetweenIsomers 
                rt(2) = rt1(2);
            elseif membpos(1) ~= 1 && membpos(2) ~= 1 && rt2(1) > 0
                rt(2) = rt2(1);
            elseif membpos(1) ~= 2 && membpos(2) ~= 2 && rt2(2) > 0
                rt(2) = rt2(2); 
            end
        end

        if abs(rt(2) - rt(1)) > maxRTBetweenIsomers
            rt(2) = 0;
        end
        
        fprintf(fileID, "%s", productName(i, j));
        
        if rt(1) && rt(2)
            fprintf(fileID, "\t%.2f\t%.2f", rt);
            yl = ylim();
            text(7, 0.5*yl(2), sprintf("RT: %.2f, RT2: %.2f", rt), 'FontSize', 14, 'Color', 'black', 'FontWeight','bold');

        elseif rt(1)
            fprintf(fileID, "\t%.2f", rt(1));
            yl = ylim();
            text(7, 0.5*yl(2), sprintf("RT: %.2f", rt(1)), 'FontSize', 14, 'Color', 'black', 'FontWeight','bold');
        end
        
        fprintf(fileID, "\n");     

        if or(rt(1), rt(2))
            hold off;
            drawnow nocallbacks;       
           
            exportgraphics(fig, outputPDF, 'Append', true);
            
            clf;
        end
    end
end

close(fig);
fclose(fileID);
