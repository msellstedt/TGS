function [pName, pMass, pos] = searchedMass(m, mainBB, dim, productName, productMass)
   
    global massTolerance;

    pName = "";
    pMass = 0;
    pos = 0;

    if dim == 1
        for i=1:size(productMass, 1)
            if m + massTolerance > productMass(i, mainBB) & m - massTolerance < productMass(i, mainBB)
                pName = productName(i, mainBB);
                pMass = productMass(i, mainBB);
                pos = i;
                return;
            end
        end
    end

    if dim == 2
        for j=1:size(productMass, 2)
            if m + massTolerance > productMass(mainBB, j) & m - massTolerance < productMass(mainBB, j)
                pName = productName(mainBB, j);
                pMass = productMass(mainBB, j);
                pos = j;
                return;
            end
        end
    end
end
