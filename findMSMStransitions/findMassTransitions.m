% ______________________________________________________________________
%
% MANUALLY CURATE THE RETENTION TIME DATA BEFORE PROCEEDING WITH THE MSMS-ANALYSIS !
% ______________________________________________________________________

load MSvariables.mat

combineMSMSData      % Combines and saves MSMS-data into one MSMS-spectra per product and dataset.   

MSMSAnalysis;        % Plots MSMS-data and outputs four MSMS-transitions per compound if MSMS-data is present
