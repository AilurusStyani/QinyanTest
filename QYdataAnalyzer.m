function QYdataAnalyzer(path)
if nargin < 1
    path = fullfile(pwd,'data');
end

files = dir(fullfile(path,'*.txt'));

for filei = 1:length(files)
    fileSaveName = strrep(files(filei).name, '.txt', '.mat');
    fileH = fopen(fullfile(files(filei).folder, files(filei).name),'r');
    validResult = {};
    pointNum = 0;
    fixData = {};
    while feof(fileH)~=-1
        lineData = fgetl(fileH);
        if lineData == -1
            break
        end
        if contains(lineData,'get validate values')
            resultNum = str2num(lineData(strfind(lineData, 'values')+6:end));
            validResult = cat(1,validResult,resultNum);
        elseif contains(lineData,'mark: point in')
            pointNum = pointNum+1;
            pointLoc = str2num(lineData(strfind(lineData,'point in')+8:end));
            fixData{pointNum,2} = pointLoc;
        end
        if pointNum>0
            fixData{pointNum,1} = cat(1,fixData{pointNum,1},{lineData});
        end
    end
    fclose(fileH);
    save(['extracted_' fileSaveName], 'fixData','validResult');
end
            