
clearvars
global massTolerance;
global noCollisionEnergies;
global productName;

massTolerance = 0.003;            % mass tolerance (delta m/z) between recorded and found mass when searching for MS-data
noCollisionEnergies = 3;          % number of collision energies collected in the data
massAdduct = 1.00728;             % search for [M+H]+ ions (data in positive mode ionization)
mzMLpath = "mzML\";               % path to mzML data files. 

% The mzML files should contain data from a reaction mixture of one builidning block
% combined with all building blocks from the complementary reaction set. 
% The naming of the files should be as exemplified by A1Z.mzML - it containes
% the reaction mixture when builing block A1 has been reacted with all
% building blocks Z.

[productMass, productName, bb1, bb2, bbMasses] = generateProductMasses("Buildingblock masses.xlsx", massAdduct);    % generates names and M+massAdduct masses for all combination of builidingblocks defined in the given file

triazoleMSMSProd1(size(bb1, 1), size(bb2, 1)) = triazoleMSMS; % declare MSMSdata variables for products from dataset1 (each buildingblock1 combined with all buildingblocks2)
triazoleMSMSProd2(size(bb1, 1), size(bb2, 1)) = triazoleMSMS; % declare MSMSdata variables for product from dataset2 (each buildingblock2 combined with all buildingblocks1)

triazoleMSProd1(size(bb1, 1), size(bb2, 1)) = triazoleMS; % declare MSdata variables for products from dataset1 (each buildingblock1 combined with all builigblocks2)
triazoleMSProd2(size(bb1, 1), size(bb2, 1)) = triazoleMS; % declare MSdata variables for products from dataset2 (each buildingblock2 combined with all builigblocks1)

bbLett = cell2mat(bb2(1));

for m = 1:size(bb1, 1)
    filename = strcat(mzMLpath, bb1(m), bbLett(1), ".mzml");
    
    disp(strcat("Reading data ", filename));

    [triazoleMSMSProd1(m, 1:size(bb2, 1)), triazoleMSProd1(m, 1:size(bb2, 1))] = extractMSMSDataMZML(filename, m, 2, productName, productMass);    
end

bbLett = cell2mat(bb1(1));
for n = 1:size(bb2, 1)
    filename = strcat(mzMLpath, bb2(n), bbLett(1), ".mzml");

    disp(strcat("Reading data ", filename));
    
    [triazoleMSMSProd2(1:size(bb1, 1), n), triazoleMSProd2(1:size(bb1, 1), n)] = extractMSMSDataMZML(filename, n, 1, productName, productMass);    
end

RTAnalysis;             % plot and save data for verification of retentions times for all products

save MSvariables.mat;
