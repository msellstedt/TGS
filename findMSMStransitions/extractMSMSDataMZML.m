function [MSMSdata, MSdata] = extractMSMSDataMZML(filename, mainBB, dim, productName, productMass)
    
    global massTolerance;
 
    xmlData = readstruct(filename, "FileType", "xml");          % read mzML-file into xmlData
    noSpectra = xmlData.mzML.run.spectrumList.countAttribute;
    
    precursor = zeros(noSpectra, 1);
    rt = zeros(noSpectra, 1);
    cE = zeros(noSpectra, 1);
    mzArr = zeros(noSpectra, 10000);
    intArr = zeros(noSpectra, 10000);
    
    for i = 1:noSpectra         % convert each spectra in the file and save rt, cE, mzArr, and intArr 
        
        if isa(xmlData.mzML.run.spectrumList.spectrum(i).precursorList, 'struct')
            precursor(i) = xmlData.mzML.run.spectrumList.spectrum(i).precursorList.precursor.selectedIonList.selectedIon.cvParam(1).valueAttribute;
            cE(i) = xmlData.mzML.run.spectrumList.spectrum(i).precursorList.precursor.activation.cvParam(2).valueAttribute;
        end
        
        rt(i) = xmlData.mzML.run.spectrumList.spectrum(i).scanList.scan.cvParam.valueAttribute;
        
        tmpMZArr = xmlData.mzML.run.spectrumList.spectrum(i).binaryDataArrayList.binaryDataArray(1).binary;
        tmpIntArr = xmlData.mzML.run.spectrumList.spectrum(i).binaryDataArrayList.binaryDataArray(2).binary;
    
        %convert binary base64 data to double arrays
        tmpMZDecoded = typecast(matlab.net.base64decode(tmpMZArr), 'double');
        tmpIntDecoded = typecast(matlab.net.base64decode(tmpIntArr), 'double');
        mzArr(i, 1:length(tmpMZDecoded)) = tmpMZDecoded;
        intArr(i, 1:length(tmpIntDecoded)) = tmpIntDecoded;
    end
     
    MSMSdata(size(productMass, dim), 1) = triazoleMSMS;
    MSdata(size(productMass, dim), 1) = triazoleMS;
       
    for i = 1:noSpectra
        pmass = precursor(i); 
        if pmass == 0       % not MSMS-data
            continue;
        end
            
        [prodName, prodMass, pos] = searchedMass(pmass, mainBB, dim, productName, productMass);   % check if pmass is an expected product mass in current file
        
        if prodName ~= ""  % the product mass pmass corresponds to a desired product present in the current file        
            MSMSdata(pos).noOccurences = MSMSdata(pos).noOccurences + 1;
            MSMSdata(pos).name = prodName;
            MSMSdata(pos).theoreticalMass = prodMass;
            MSMSdata(pos).foundMass(MSMSdata(pos).noOccurences) = pmass;
            
            MSMSdata(pos).collisionEnergy(MSMSdata(pos).noOccurences) = cE(i);
            
            MSMSdata(pos).fragmentMass(MSMSdata(pos).noOccurences, 1:12) = zeros(12, 1);
            MSMSdata(pos).fragmentIntensity(MSMSdata(pos).noOccurences, 1:12) = zeros(12, 1);

            MSMSdata(pos).retentionTimes(MSMSdata(pos).noOccurences) = rt(i);   % save rt in minutes    

            [MSMSdata(pos).fragmentIntensity(MSMSdata(pos).noOccurences, 1:12), Ipos] = maxk(intArr(i, :), 12);
            MSMSdata(pos).fragmentMass(MSMSdata(pos).noOccurences, 1:12) = mzArr(i, Ipos);             
        
        end
    end 

    if dim == 1
        for a=1:size(productMass, 1)
            
            peak = 1;
            MSdata(a).name = productName(a, mainBB);
            MSdata(a).theoreticalMass = productMass(a, mainBB);
            
            for i = 1:noSpectra
                
                pmass = precursor(i); 
                if pmass ~= 0       % not MS-data
                    continue;
                end

                MSdata(a).retentionTimes(peak) = rt(i);

                index = find(abs(mzArr(i, :) - productMass(a, mainBB)) < massTolerance);
                
                if ~isempty(index)
                    [MSdata(a).intensity(peak), pos] = max(intArr(i, index));
                
                    if MSdata(a).intensity(peak) == max(MSdata(a).intensity)
                        MSdata(a).foundMass = mzArr(i, index(pos));               % extract mass at highest intensity
                    end
                else
                    MSdata(a).intensity(peak) = 0;          
                end
                
                peak = peak + 1;
            end
        end
             
    else
        for a=1:size(productMass, 2)
            
            peak = 1;
            MSdata(a).name = productName(mainBB, a);
            MSdata(a).theoreticalMass = productMass(mainBB, a);
            
            for i = 1:noSpectra
                
                pmass = precursor(i); 
                if pmass ~= 0       % not MS-data
                    continue;
                end

                MSdata(a).retentionTimes(peak) = rt(i);
                   
                index = find(abs(mzArr(i, :) - productMass(mainBB, a)) < massTolerance);
                
                if ~isempty(index)
                    
                    [MSdata(a).intensity(peak), pos] = max(intArr(i, index));
                
                    if MSdata(a).intensity(peak) == max(MSdata(a).intensity) 
                        MSdata(a).foundMass = mzArr(i, index(pos));              % extract mass at highest intensity
                    end
                else
                    MSdata(a).intensity(peak) = 0;
                end
                
                peak = peak + 1;
            end
        end
    
    end
   
end