function [peaks, intensities] = addFragment(peaks1, intensities1, peaks2, intensities2)
    
    tmpIntensities = zeros(size(peaks1, 1) + size(peaks2, 1), 1);
    
    intensities = cat(1, intensities1, intensities2);
    peaks = cat(1, peaks1, peaks2);
    peaks = round(peaks, 1);                        % round masses to 1 decimal

    [tmpPeaks, ~, itmp] = unique(peaks, 'stable');  % combine peaks with same mass (with 1 decimal) from different recordings
    
    for a = 1:size(tmpPeaks, 1)
        tmpIntensities(a) = sum(intensities(itmp == a));
    end
    
    %sort peaks by intensity and keep those with highest intesity
    [tintensities, tpeaks] = sort(tmpIntensities, 'descend');
    
    peaks = zeros(size(peaks1, 1), 1);
    peaks(1:min(size(tmpPeaks, 1), size(peaks1, 1))) = tmpPeaks(tpeaks(1:min(size(tmpPeaks, 1), size(peaks1, 1))));
    
    intensities = tintensities(1:size(intensities1,1));
    
end