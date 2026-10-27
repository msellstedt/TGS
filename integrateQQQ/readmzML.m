function spectraData = readmzML(mzMLFile)

    global maxNoTimePoints;
    spectraData = spectraFile;

    xmlData = readstruct(mzMLFile, "FileType", "xml");
    spectraData.noChromatograms = size(xmlData.mzML.run.chromatogramList.chromatogram, 2) - 1;
    spectraData.timeArr = zeros(spectraData.noChromatograms, maxNoTimePoints);
    spectraData.intArr = zeros(spectraData.noChromatograms, maxNoTimePoints);
    spectraData.precIon = zeros(spectraData.noChromatograms, 1);
    spectraData.prodIon = zeros(spectraData.noChromatograms, 1);
    spectraData.CE = zeros(spectraData.noChromatograms, 1);
    spectraData.TIC = zeros(maxNoTimePoints, 1);
    spectraData.TICsize = 0;
    spectraData.delay = 0;
    spectraData.baseLine = zeros(spectraData.noChromatograms, 1);
    spectraData.estRt = zeros(spectraData.noChromatograms, 1);
    spectraData.compNr = zeros(spectraData.noChromatograms, 1);
    spectraData.primary = zeros(spectraData.noChromatograms, 1);     
    spectraData.name = strings(spectraData.noChromatograms, 1);
    spectraData.shortDescr = "";
    spectraData.area = zeros(spectraData.noChromatograms, 3);     

    for i = 2:spectraData.noChromatograms+1   % skip chromatogram 1 (= TIC)
    
        if ~isa(xmlData.mzML.run.chromatogramList.chromatogram(i).precursor, 'struct')  % only read MRM-data
            spectraData.noChromatograms = i-2;
            break;
        end

        spectraData.precIon(i-1) = xmlData.mzML.run.chromatogramList.chromatogram(i).precursor.isolationWindow.cvParam.valueAttribute;
        spectraData.prodIon(i-1) = xmlData.mzML.run.chromatogramList.chromatogram(i).product.isolationWindow.cvParam.valueAttribute;
        spectraData.CE(i-1) = xmlData.mzML.run.chromatogramList.chromatogram(i).precursor.activation.cvParam(2).valueAttribute;
    
        tmpTimeArr = xmlData.mzML.run.chromatogramList.chromatogram(i).binaryDataArrayList.binaryDataArray(1).binary;
        tmpIntArr = xmlData.mzML.run.chromatogramList.chromatogram(i).binaryDataArrayList.binaryDataArray(2).binary;
    
        %convert binary base64 data to arrays
        tmpT = typecast(matlab.net.base64decode(tmpTimeArr), 'double');
        tmpI = typecast(matlab.net.base64decode(tmpIntArr), 'single');
        spectraData.estRt(i-1) = median(tmpT);
        spectraData.baseLine(i-1) = min(tmpI);
   
        % scale data to 1200 data points per min
        deltaT = tmpT(end)-tmpT(1);
        noDataPoints = floor(1200 * deltaT) + 1;
        tmpIntArr = interp1(tmpT, tmpI, linspace(tmpT(1), tmpT(end), noDataPoints), "linear");
        tmpTimeArr = linspace(tmpT(1), tmpT(end), noDataPoints);    
        
        deltaT = tmpTimeArr(2)-tmpTimeArr(1);
        startPos = floor(tmpTimeArr(1)/deltaT) + 1;
        endPos = startPos + noDataPoints - 1;
    
        spectraData.timeArr(i-1, startPos:endPos) = tmpTimeArr;
        spectraData.intArr(i-1, startPos:endPos) = tmpIntArr - spectraData.baseLine(i-1);  % remove baseline 
        
        if endPos > spectraData.TICsize
            spectraData.TICsize = endPos;
        end

        spectraData.TIC(startPos:endPos) = spectraData.TIC(startPos:endPos) + transpose(spectraData.intArr(i-1, startPos:endPos));        
    end

end