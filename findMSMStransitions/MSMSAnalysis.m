global noCollisionEnergies;

fig = figure(1);

subplot(2, noCollisionEnergies, 1);
set(fig,'Position',[40,20,1000,600]);
clf;

M = readtable("Excluded fragment losses.xlsx");
exclusionLoss = table2array(M(1:size(M, 1), 1));    % read list of fragment mass losses to exclude

pause('on');

fileID = fopen("Transitions.txt", "w");

outputPDF = "MSMSAnalysis.pdf";
delete(outputPDF);

for i = 1:size(MSMSDataSummary, 1)
    for j = 1:size(MSMSDataSummary, 2)
        dataToSave = 0;

        highestCommonInt = zeros(12, 1);
        highestCommonMass = zeros(12, 1);
        cEn = zeros(12, 1);

        for e = 1:noCollisionEnergies
            % normalize intensity data
            norm1 = sum(MSMSDataSummary(i, j, e).fragmentIntensityPeak1(1:12));
            if norm1 ~= 0 
                MSMSDataSummary(i, j, e).fragmentIntensityPeak1(1:12) = MSMSDataSummary(i, j, e).fragmentIntensityPeak1(1:12) ./ norm1;
            end

            norm2 = sum(MSMSDataSummary(i, j, e).fragmentIntensityPeak2(1:12));
            if norm2 ~= 0 
                MSMSDataSummary(i, j, e).fragmentIntensityPeak2(1:12) = MSMSDataSummary(i, j, e).fragmentIntensityPeak2(1:12) ./ norm2;
            end

            if MSMSDataSummary(i, j, e).retentionTimePeak1 ~= 0 && MSMSDataSummary(i, j, e).retentionTimePeak2 ~= 0
                
                % MSMS-data exists in both datasets
                
                % Find peaks present in both datasets
                [commonMass, p1, p2] = intersect(round(MSMSDataSummary(i, j, e).fragmentMassPeak1, 1), round(MSMSDataSummary(i, j, e).fragmentMassPeak2, 1), 'stable');
                
                combInt = MSMSDataSummary(i, j, e).fragmentIntensityPeak1(p1) + MSMSDataSummary(i, j, e).fragmentIntensityPeak2(p2);
                [cint, cpos] = maxk(combInt, 12);
                
                for a = 1:size(cint, 1)
                    
                    if ~ismember(round(MSMSDataSummary(i, j, e).theoreticalMass - commonMass(cpos(a)), 1), exclusionLoss(:)) % exclude mass losses defined in exclusionLoss
                        
                        for b = 1:12
                            if cint(a) > highestCommonInt(b)
                                if b < 12
                                    highestCommonInt(b+1:12) = highestCommonInt(b:11);
                                    highestCommonMass(b+1:12) = highestCommonMass(b:11);
                                    cEn(b+1:12) = cEn(b:11);
                                end
                                highestCommonInt(b) = cint(a);                 
                                highestCommonMass(b) = commonMass(cpos(a));
                                cEn(b) = e;           
                                break;
                            end
                        end
                    end
           
                end    
                
                % remove fragments present several times in highestCommonMass
                
                [highestCommonMass, locs] = unique(highestCommonMass, 'stable');
                highestCommonInt = highestCommonInt(locs);
                cEn = cEn(locs);
                highestCommonMass = cat(1, highestCommonMass, zeros(12-size(highestCommonMass, 1), 1));
                highestCommonInt = cat(1, highestCommonInt, zeros(12-size(highestCommonInt, 1), 1));
                cEn = cat(1, cEn, zeros(12-size(cEn, 1), 1));
         
            end
            
            if MSMSDataSummary(i, j, e).retentionTimePeak1 ~= 0
                % MSMS-data present in dataset 1

                dataToSave = 1;
               
                % Plot the MSMS-data
         
                subplot(2, noCollisionEnergies, e);
                stem(MSMSDataSummary(i, j, e).fragmentMassPeak1, MSMSDataSummary(i, j, e).fragmentIntensityPeak1, 'Marker','none', 'LineWidth', 1);
         
                title(sprintf("%s, %.4f @ %.2f min, %.0f eV", MSMSDataSummary(i, j, e).name, MSMSDataSummary(i, j, e).foundParentMass1, MSMSDataSummary(i, j, e).retentionTimePeak1, MSMSDataSummary(i, j, e).collisionEnergy))
                axis([40, MSMSDataSummary(i, j, e).foundParentMass1+20, 0, 1.1 * max(MSMSDataSummary(i, j, e).fragmentIntensityPeak1)]);
                
                for p=1:6
                    text(MSMSDataSummary(i, j, e).fragmentMassPeak1(p), MSMSDataSummary(i, j, e).fragmentIntensityPeak1(p), sprintf("%.1f", MSMSDataSummary(i, j, e).fragmentMassPeak1(p)), "FontSize", 8);
                end

                for p=1:3        
                    if e == noCollisionEnergies && cEn(p) ~= 0     % print massdata that exists in both datasets
                        yl = ylim();
                        xl = xlim();
                        text(xl(1)-80 + (p-1)*(50 + xl(2) - xl(1))/2, -yl(2)/6, sprintf("%.1f @ %.0f eV", highestCommonMass(p), MSMSDataSummary(i, j, cEn(p)).collisionEnergy), "FontSize", 12, "FontWeight","bold");                        
                    end

                end
               
            end
            
            if MSMSDataSummary(i, j, e).retentionTimePeak2 ~= 0
                % MSMS-data present in dataset 2
                
                dataToSave = 1;
                
                % Plot the MSMS-data
                
                subplot(2, noCollisionEnergies, e + noCollisionEnergies);
                stem(MSMSDataSummary(i, j, e).fragmentMassPeak2, MSMSDataSummary(i, j, e).fragmentIntensityPeak2, 'Marker','none', 'LineWidth', 1);
                
                title(sprintf("%s, %.4f @ %.2f min, %.0f eV", MSMSDataSummary(i, j, e).name, MSMSDataSummary(i, j, e).foundParentMass2, MSMSDataSummary(i, j, e).retentionTimePeak2, MSMSDataSummary(i, j, e).collisionEnergy))
                axis([40, MSMSDataSummary(i, j, e).foundParentMass2+20, 0, 1.1 * max(MSMSDataSummary(i, j, e).fragmentIntensityPeak2)]);
                
                for p=1:6
                    text(MSMSDataSummary(i, j, e).fragmentMassPeak2(p), MSMSDataSummary(i, j, e).fragmentIntensityPeak2(p), sprintf("%.1f", MSMSDataSummary(i, j, e).fragmentMassPeak2(p)), "FontSize", 8);
                end
            end

        end
        
        if dataToSave ~= 0  % save data

            if cEn(1) ~= 0 % mass data exists in both datasets            

                fprintf(fileID, "%s, rt %.2f: %.1f ->", MSMSDataSummary(i, j, cEn(1)).name, MSMSDataSummary(i, j, 1).retentionTimePeak1, MSMSDataSummary(i, j, cEn(1)).theoreticalMass);

                for a = 1:4
                    if cEn(a) ~= 0  
                        % Save the four most intense transitions present in both datasets. (Intensity relative to the sum of the top ten fragments per collision energy.) 
                        
                        fprintf(fileID, " %.1f @%.0f eV", highestCommonMass(a), MSMSDataSummary(i, j, cEn(a)).collisionEnergy);    
                    end
                end

                fprintf(fileID, "\n");
           
            else
                if MSMSDataSummary(i, j, 1).retentionTimePeak1 ~= 0
                    
                    % MSMSdata exists in only dataset 1 for the product
                    
                    [~, locs] = maxk(cat(1, MSMSDataSummary(i, j, 1).fragmentIntensityPeak1, MSMSDataSummary(i, j, 2).fragmentIntensityPeak1, MSMSDataSummary(i, j, 3).fragmentIntensityPeak1), 12);
                    
                    m = cat(1, MSMSDataSummary(i, j, 1).fragmentMassPeak1, MSMSDataSummary(i, j, 2).fragmentMassPeak1, MSMSDataSummary(i, j, 3).fragmentMassPeak1);
                    masses = m(locs);
                    
                    b = 0;
                    usedFragments = zeros(4, 1);
                    
                    fprintf(fileID, "%s, rt %.2f: %.1f ->", MSMSDataSummary(i, j, 1).name, MSMSDataSummary(i, j, 1).retentionTimePeak1, MSMSDataSummary(i, j, 1).theoreticalMass);
                    for a = 1:4     % save the four highest MS-fragments
                        b = b + 1;
                        if b > 12
                            break;
                        end

                        while or(or(ismember(round(MSMSDataSummary(i, j, 1).theoreticalMass - masses(b), 1), exclusionLoss(:)), ismember(round(masses(b), 1), usedFragments)), masses(b) >= MSMSDataSummary(i, j, 1).theoreticalMass)
                            % exclude fragments with mass losses defined in exclusionLoss 
                            % or already selected at other collision energies

                            b = b + 1;
                            if b > 12
                                break;
                            end
                        end
                        
                        if b <= 12
                            usedFragments(a) = round(masses(b), 1);
                  
                            fprintf(fileID, " %.1f @%.0f eV", masses(b), MSMSDataSummary(i, j, idivide(int32(locs(b)), 12, 'floor') + 1).collisionEnergy);    

                        end
                    end

                    fprintf(fileID, "\n");
                end

                if MSMSDataSummary(i, j, 1).retentionTimePeak2 ~= 0
                    
                    % MSMSdata exists only in dataset 2 for the product

                    [~, locs] = maxk(cat(1, MSMSDataSummary(i, j, 1).fragmentIntensityPeak2, MSMSDataSummary(i, j, 2).fragmentIntensityPeak2, MSMSDataSummary(i, j, 3).fragmentIntensityPeak2), 12);
                    
                    m = cat(1, MSMSDataSummary(i, j, 1).fragmentMassPeak2, MSMSDataSummary(i, j, 2).fragmentMassPeak2, MSMSDataSummary(i, j, 3).fragmentMassPeak2);
                    masses = m(locs);

                    b = 0;
                    usedFragments = zeros(4, 1);
           
                    fprintf(fileID, "%s, rt %.2f: %.1f ->", MSMSDataSummary(i, j, 1).name, MSMSDataSummary(i, j, 1).retentionTimePeak2, MSMSDataSummary(i, j, 1).theoreticalMass);    
                            
                    for a = 1:4     % save the four most intense fragments
                        b = b + 1;
                        if b > 12
                            break;
                        end
                        
                        while or(or(ismember(round(MSMSDataSummary(i, j, 1).theoreticalMass - masses(b), 1), exclusionLoss(:)), ismember(round(masses(b), 1), usedFragments)), masses(b) >= MSMSDataSummary(i, j, 1).theoreticalMass)
                            % exclude fragments with mass losses defined in exclusionLoss
                            % or already selected at other collision energies
                            
                            b = b + 1;
                            if b > 12
                                break;
                            end
                        end
                        
                        if b <= 12               
                            usedFragments(a) = round(masses(b), 1);
                            
                            fprintf(fileID, " %.1f @%.0f eV", masses(b), MSMSDataSummary(i, j, idivide(int32(locs(b)), 12, 'floor') + 1).collisionEnergy);    
                            
                        end
                    end

                    fprintf(fileID, "\n");
                end
            end
 
            drawnow nocallbacks;
        
            exportgraphics(fig, outputPDF, 'Append', true);    % save graphic 
            
            clf;
        end
 
    end
end

close(fig);

fclose(fileID);