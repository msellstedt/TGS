function delay = align(arr1, arr2, maxLag)

    minDiff = sum((arr1 - arr2).^2, 1);
    delay = 0;

    for i = 1:maxLag
        tmpArr = circshift(arr1, i);
        tmpDiff = sum((tmpArr - arr2).^2, 1);
      
        if tmpDiff < minDiff
            minDiff = tmpDiff;
            delay = -i;
        end
     
        tmpArr = circshift(arr1, -i);
        tmpDiff = sum((tmpArr - arr2).^2, 1);     
      
        if tmpDiff < minDiff
            minDiff = tmpDiff;
            delay = i;
        end
    
    end
end