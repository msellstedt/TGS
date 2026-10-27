% Combines MSMS-data into one MSMS-spectra per compound and dataset
% Expected retention times are read from "RT.txt"
% The twelve most intense fragments per ionization energy and product is saved in MSMSData.txt. 
% Products without recorded MSMS-data in the expected retention time windows are output in the file "MissingMSMSdata.txt"
% The twelve most common fragments from subsets of products containing a specific building block is saved in "FreqBBFragments.txt"

global productName;
global noCollisionEnergies;
rtTolerance = 0.1;              % MSMS-data is searched for in the retention time window of RT +- rtTolerance. (The RT is read from "RT.txt".)
massLimit = 1200;               % upper limit of recorded m/z

[rtName, RT] = readRT("RT.txt");  % read retention times

fileID = fopen("MSMSDataC.txt", 'w');
fileID2 = fopen("MSMSDataCompr.txt", 'w');
missingDataFileID = fopen("MissingMSMSdata.txt", 'w');

bb1FragmentMass = zeros(size(triazoleMSMSProd1, 1), massLimit, noCollisionEnergies);
bb2FragmentMass = zeros(size(triazoleMSMSProd1, 2), massLimit, noCollisionEnergies); 

clearvars MSMSDataSummary;
MSMSDataSummary(size(triazoleMSMSProd1, 1), size(triazoleMSMSProd1, 2), noCollisionEnergies) = MSMSdata;

for i = 1:size(triazoleMSMSProd1, 1)
    for j = 1:size(triazoleMSMSProd1, 2)
       
        pos = find(rtName == triazoleMSMSProd1(i, j).name);
   
        if ~isempty(pos)

            % MSMS data exists for the product (i, j) in dataset 1
            
            for r = 1:size(triazoleMSMSProd1(i,j).retentionTimes, 1)        
                rcomp = idivide(int32(r-1), noCollisionEnergies, "floor")*noCollisionEnergies + 1;

                if or(triazoleMSMSProd1(i,j).retentionTimes(rcomp) - rtTolerance < RT(pos, 1) && triazoleMSMSProd1(i,j).retentionTimes(rcomp) + rtTolerance > RT(pos, 1), ...
                      triazoleMSMSProd1(i,j).retentionTimes(rcomp) - rtTolerance < RT(pos, 2) && triazoleMSMSProd1(i,j).retentionTimes(rcomp) + rtTolerance > RT(pos, 2))
                    
                    % MSMS-data exist for product (i,j) in dataset 1 at one of the selected retention time points 
                    
                    tmpCE = mod(r-1, noCollisionEnergies) + 1;
                    validPeaks = find(triazoleMSMSProd1(i,j).fragmentMass(r,:) <= triazoleMSMSProd1(i,j).foundMass(tmpCE) + 1.1);
                  
                    %[MSMSDataSummary(i, j, tmpCE).fragmentMassPeak1(:), MSMSDataSummary(i, j, tmpCE).fragmentIntensityPeak1(:)] = addFragment(MSMSDataSummary(i, j, tmpCE).fragmentMassPeak1(:), MSMSDataSummary(i, j, tmpCE).fragmentIntensityPeak1(:), triazoleMSMSProd1(i,j).fragmentMass(r,:)', triazoleMSMSProd1(i,j).fragmentIntensity(r,:)');
                    [MSMSDataSummary(i, j, tmpCE).fragmentMassPeak1(:), MSMSDataSummary(i, j, tmpCE).fragmentIntensityPeak1(:)] = addFragment(MSMSDataSummary(i, j, tmpCE).fragmentMassPeak1(:), MSMSDataSummary(i, j, tmpCE).fragmentIntensityPeak1(:), triazoleMSMSProd1(i,j).fragmentMass(r, validPeaks)', triazoleMSMSProd1(i,j).fragmentIntensity(r, validPeaks)');
                    
                    MSMSDataSummary(i, j, tmpCE).name = rtName(pos);
                    MSMSDataSummary(i, j, tmpCE).collisionEnergy = triazoleMSMSProd1(i,j).collisionEnergy(r);
                    MSMSDataSummary(i, j, tmpCE).theoreticalMass = triazoleMSMSProd1(i,j).theoreticalMass;
                    MSMSDataSummary(i, j, tmpCE).foundParentMass1 = triazoleMSMSProd1(i,j).foundMass(r);        
                    
                    if triazoleMSMSProd1(i, j).retentionTimes(rcomp) - rtTolerance < RT(pos, 1) && triazoleMSMSProd1(i,j).retentionTimes(rcomp) + rtTolerance > RT(pos, 1)
                        MSMSDataSummary(i, j, tmpCE).retentionTimePeak1 = RT(pos, 1);
                    else
                        MSMSDataSummary(i, j, tmpCE).retentionTimePeak1 = RT(pos, 2);
                    end
      
                    bb1FragmentMass(i, int32(triazoleMSMSProd1(i, j).fragmentMass(r,:)), tmpCE) = bb1FragmentMass(i, int32(triazoleMSMSProd1(i,j).fragmentMass(r,:)), tmpCE) + triazoleMSMSProd1(i,j).fragmentIntensity(r,:) / max(triazoleMSMSProd1(i,j).fragmentIntensity(r));
                    bb2FragmentMass(j, int32(triazoleMSMSProd1(i, j).fragmentMass(r,:)), tmpCE) = bb2FragmentMass(j, int32(triazoleMSMSProd1(i,j).fragmentMass(r,:)), tmpCE) + triazoleMSMSProd1(i,j).fragmentIntensity(r,:) / max(triazoleMSMSProd1(i,j).fragmentIntensity(r)); 
                
                end
            end
        end
        
        pos = find(rtName == triazoleMSMSProd2(i, j).name);
        
        if ~isempty(pos)
            
            % MSMS-data for product (i,j) exists in dataset 2
            
            for r = 1:size(triazoleMSMSProd2(i, j).retentionTimes, 1)
                rcomp = idivide(int32(r-1), noCollisionEnergies, "floor")*noCollisionEnergies + 1;

                if or(triazoleMSMSProd2(i, j).retentionTimes(rcomp) - rtTolerance < RT(pos, 1) && triazoleMSMSProd2(i,j).retentionTimes(rcomp) + rtTolerance > RT(pos, 1), ...
                      triazoleMSMSProd2(i, j).retentionTimes(rcomp) - rtTolerance < RT(pos, 2) && triazoleMSMSProd2(i,j).retentionTimes(rcomp) + rtTolerance > RT(pos, 2))
                    
                    % MSMS-data exist for product (i,j) in dataset 2 at one of the selected retention time points 
                    
                    tmpCE = mod(r-1, noCollisionEnergies) + 1;
                    validPeaks = find(triazoleMSMSProd2(i,j).fragmentMass(r,:) <= triazoleMSMSProd2(i,j).foundMass(tmpCE) + 1.1);
                    
                    %[MSMSDataSummary(i, j, tmpCE).fragmentMassPeak2(:), MSMSDataSummary(i, j, tmpCE).fragmentIntensityPeak2(:)] = addFragment(MSMSDataSummary(i, j, tmpCE).fragmentMassPeak2(:), MSMSDataSummary(i, j, tmpCE).fragmentIntensityPeak2(:), triazoleMSMSProd2(i,j).fragmentMass(r,:)', triazoleMSMSProd2(i,j).fragmentIntensity(r,:)');
                    [MSMSDataSummary(i, j, tmpCE).fragmentMassPeak2(:), MSMSDataSummary(i, j, tmpCE).fragmentIntensityPeak2(:)] = addFragment(MSMSDataSummary(i, j, tmpCE).fragmentMassPeak2(:), MSMSDataSummary(i, j, tmpCE).fragmentIntensityPeak2(:), triazoleMSMSProd2(i,j).fragmentMass(r,validPeaks)', triazoleMSMSProd2(i,j).fragmentIntensity(r,validPeaks)');
                    
                    MSMSDataSummary(i, j, tmpCE).name = rtName(pos);
                    MSMSDataSummary(i, j, tmpCE).collisionEnergy = triazoleMSMSProd2(i,j).collisionEnergy(r);
                    MSMSDataSummary(i, j, tmpCE).theoreticalMass = triazoleMSMSProd2(i,j).theoreticalMass;
                    MSMSDataSummary(i, j, tmpCE).foundParentMass2 = triazoleMSMSProd2(i,j).foundMass(r);        
                    
                    if triazoleMSMSProd2(i, j).retentionTimes(rcomp) - rtTolerance < RT(pos, 1) && triazoleMSMSProd2(i,j).retentionTimes(rcomp) + rtTolerance > RT(pos, 1)
                        MSMSDataSummary(i, j, tmpCE).retentionTimePeak2 = RT(pos, 1);
                    else
                        MSMSDataSummary(i, j, tmpCE).retentionTimePeak2 = RT(pos, 2);
                    end
     
                    bb1FragmentMass(i, int32(triazoleMSMSProd2(i, j).fragmentMass(r,:)), tmpCE) = bb1FragmentMass(i, int32(triazoleMSMSProd2(i,j).fragmentMass(r,:)), tmpCE) + triazoleMSMSProd2(i,j).fragmentIntensity(r,:) / max(triazoleMSMSProd2(i,j).fragmentIntensity(r));
                    bb2FragmentMass(j, int32(triazoleMSMSProd2(i, j).fragmentMass(r,:)), tmpCE) = bb2FragmentMass(j, int32(triazoleMSMSProd2(i,j).fragmentMass(r,:)), tmpCE) + triazoleMSMSProd2(i,j).fragmentIntensity(r,:) / max(triazoleMSMSProd2(i,j).fragmentIntensity(r));
  
                end
            end
        end

    end
end

rt1 = 0;
rt2 = 0;
noMSMSdata = 0;
for i = 1:size(MSMSDataSummary, 1)

    for j = 1:size(MSMSDataSummary, 2)
    
        for r = 1:size(MSMSDataSummary, 3)
            if MSMSDataSummary(i, j, r).retentionTimePeak1 ~= 0    
                fprintf(fileID, "%s @ %.2f min, %.0f eV, calc: %.4f, found in dataset 1: %.4f\n Fragment Mass: %.3f %.3f %.3f %.3f %.3f %.3f %.3f %.3f %.3f %.3f %.3f %.3f", MSMSDataSummary(i, j, r).name, MSMSDataSummary(i, j, r).retentionTimePeak1, MSMSDataSummary(i, j, r).collisionEnergy, MSMSDataSummary(i, j, r).theoreticalMass, MSMSDataSummary(i, j, r).foundParentMass1, MSMSDataSummary(i, j, r).fragmentMassPeak1(:));
                fprintf(fileID, "\n Fragment Intensity: %.0f %.0f %.0f %.0f %.0f %.0f %.0f %.0f %.0f %.0f %.0f %.0f\n", MSMSDataSummary(i, j, r).fragmentIntensityPeak1(:));                  
                rt1 = rt1 + 1;

                fprintf(fileID2, "%s %.3f %.3f %.3f %.3f %.3f %.3f %.3f %.3f %.3f %.3f %.3f %.3f %.0f %.0f %.0f %.0f %.0f %.0f %.0f %.0f %.0f %.0f %.0f %.0f\n", MSMSDataSummary(i, j, r).name, MSMSDataSummary(i, j, r).fragmentMassPeak1(1:12), MSMSDataSummary(i, j, r).fragmentIntensityPeak1(1:12));
            end
            
            if MSMSDataSummary(i, j, r).retentionTimePeak2 ~= 0        
                fprintf(fileID, "%s @ %.2f min, %.0f eV, calc: %.4f, found in dataset 2: %.4f\n Fragment Mass: %.3f %.3f %.3f %.3f %.3f %.3f %.3f %.3f %.3f %.3f %.3f %.3f", MSMSDataSummary(i, j, r).name, MSMSDataSummary(i, j, r).retentionTimePeak2, MSMSDataSummary(i, j, r).collisionEnergy, MSMSDataSummary(i, j, r).theoreticalMass, MSMSDataSummary(i, j, r).foundParentMass2, MSMSDataSummary(i, j, r).fragmentMassPeak2(:));
                fprintf(fileID, "\n Fragment Intensity: %.0f %.0f %.0f %.0f %.0f %.0f %.0f %.0f %.0f %.0f %.0f %.0f\n", MSMSDataSummary(i, j, r).fragmentIntensityPeak2(:));                  
                rt2 = rt2 + 1;
 
                fprintf(fileID2, "%s %.3f %.3f %.3f %.3f %.3f %.3f %.3f %.3f %.3f %.3f %.3f %.3f %.0f %.0f %.0f %.0f %.0f %.0f %.0f %.0f %.0f %.0f %.0f %.0f\n", MSMSDataSummary(i, j, r).name, MSMSDataSummary(i, j, r).fragmentMassPeak2(1:12), MSMSDataSummary(i, j, r).fragmentIntensityPeak2(1:12));
            end
       
            if r == 1 && MSMSDataSummary(i, j, r).retentionTimePeak1 == 0 && MSMSDataSummary(i, j, r).retentionTimePeak2 == 0
                fprintf(fileID, "%s not found\n", productName(i, j));

                fprintf(missingDataFileID, "%s\n", productName(i, j));
      
                fprintf(fileID2, "%s 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0\n", MSMSDataSummary(i, j, r).name);
      
                noMSMSdata = noMSMSdata + 1;
            end
           
        end
    end
end

fclose(fileID);
fclose(fileID2);
fclose(missingDataFileID);
