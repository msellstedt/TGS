classdef MSMSdata
    properties
        name char = ""
        collisionEnergy int32 = 0
        theoreticalMass double = 0
        foundParentMass1 double = 0
        foundParentMass2 double = 0
        retentionTimePeak1 double = 0               % retention time of MSMS-data in dataset 1
        retentionTimePeak2 double = 0               % retention time of MSMS-data in dataset 2
        fragmentMassPeak1(12, 1) double = 0         % mass of the twelve most intense fragment peaks in dataset 1 at the selected retention time
        fragmentIntensityPeak1(12, 1) double = 0    % intensities of the corresponing fragment peaks in dataset 1
        fragmentMassPeak2(12, 1) double = 0
        fragmentIntensityPeak2(12, 1) double = 0
     
    end
end