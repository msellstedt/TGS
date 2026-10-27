function [name, RT] = readRT(filename)
    fileID = fopen(filename, 'r');
    if fileID == -1
        msgbox(strcat("No such file: ", filename));
        return;
    end
    
    data = textscan(fileID, "%s%f%f", 'Delimiter', '\t');
    RT(:,1) = double(data{2});
    RT(1:size(data{3}, 1), 2) = double(data{3});
    RT(isnan(RT))=0;
    name = string(data{1});
    fclose(fileID);
end