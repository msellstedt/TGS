function [masses, names, bb1, bb2, bbMasses] = generateProductMasses(f, massAdduct)
 
    M = readtable(f);
    masses1 = table2array(M(1:size(M,1),"ExactMass1"));
    masses2 = table2array(M(1:size(M,1),"ExactMass2"));
    names1 = table2array(M(1:size(M,1),"Buildingblock1"));
    names2 = table2array(M(1:size(M,1),"Buildingblock2"));
    
    size1 = size(M,1);
    size2 = size(M,1);

    for i = 1:size(M,1)
        if isnan(masses1(i))
            size1 = i-1;
            break;
        end
    end

    for i = 1:size(M,1)
        if isnan(masses2(i))
            size2 = i-1;
            break;
        end
    end
    
    masses(1:size1, 1:size2) = 0;
    names(1:size1, 1:size2) = "";
    bbMasses = zeros(size1+size2, 1);
    bbMasses(1:size1) = masses1(1:size1) + massAdduct;
    bbMasses(size1 + 1:size1 + size2) = masses2(1:size2) + massAdduct;
    
    for i = 1:size1
        for j = 1:size2
            masses(i, j) = masses1(i) + masses2(j) + massAdduct; 
            names(i, j) = strcat(names1(i), names2(j));
        end
    end

    bb1 = names1(1:size1);
    bb2 = names2(1:size2);

end