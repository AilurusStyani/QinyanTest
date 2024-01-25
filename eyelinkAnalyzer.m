clear all; clc; close all

edf2asc = false; % convert edf file to asc file?
asc2mat = true; % convert asc file to mat file?

%% convert .edf file to .asc file
if edf2asc
    edf2ascConvert(fullfile(pwd,'data'),pwd);
end

%% convert .asc file to .mat file
if asc2mat
    asc2matConvert(fullfile(pwd,'data'),true);
end

%% remove eyelink converted data from mat file list
file = dir(fullfile(pwd,'data','*.mat'));
eyeDataFile = dir(fullfile(pwd,'data','Converted_*.mat'));
delIndex = [];
for k = 1:length(file)
    for j = 1:length(eyeDataFile)
        if isequal(eyeDataFile(j).name,file(k).name)
            delIndex = cat(1,delIndex,k);
        end
    end
end
file(delIndex) = [];

%% start analysis
CR = cell(length(file),3);
color = {[0 0 0] [1 0 0] [0 1 0] [0 0 1]};
for filei = 1:length(file)
    data = load(fullfile(file(filei).folder, file(filei).name));
    eyeData = load(fullfile(file(filei).folder, ['Converted_' file(filei).name]));
    
    % shuffle for trial validation
    data.choice = [data.choice,reshape(1:data.info.trialNum,[],1)];
    behaviourDataValid = data.choice(~isnan(data.choice(:,1)),2);
    validIndex = intersect(behaviourDataValid,eyeData.validIndex);
    
    data.RT = [data.RT(validIndex),validIndex];
    data.choice = [data.choice(validIndex),validIndex];
    data.info.conditionList = [data.info.conditionList(validIndex),validIndex];
    
    %     if ishandle(2*filei-1); close (2*filei-1); end;  figure(2*filei-1);   set(gcf,'color','white');
    conditionIndex = unique(data.info.conditionList(:,1));
    
    conIndex = cell(1,length(conditionIndex));
    chIndex = cell(1,length(conditionIndex));
    Ch = cell(1,length(conditionIndex));
    correctRate = nan(1,length(conditionIndex));
    
    for j = 1:numel(conditionIndex)
        conIndex{j} = data.info.conditionList(:,1) == conditionIndex(j);
        Ch{j} = data.choice(conIndex{j},1);
        if ismember(conditionIndex(j),[1,2])
            correctRate(j) = sum(Ch{j}==1)/length(Ch{j});
        elseif ismember(conditionIndex(j),[3,4])
            correctRate(j) = sum(Ch{j}==2)/length(Ch{j});
        end
    end
    
    % behavior result for correct rate and repetition for each conditions
    CR{filei,1} = correctRate;
    CR{filei,2} = cellfun(@length,Ch);
    CR{filei,3} = file(filei).name;
    
    if ishandle(2*filei); close (2*filei); end;  figure(2*filei);   set(gcf,'color','white');hold on;
    
    %% analysis eye data
    for coni = 1:5
        if ~ismember(coni,[2 5])
            con = data.info.conditionList(:,1) == coni;
        elseif coni == 2
            con = data.info.conditionList(:,1) == 2;
            con = and(con,data.choice(:,1) == 1);
        elseif coni == 5
            con = data.info.conditionList(:,1) == 2;
            con = and(con,data.choice(:,1) == 2);
        end
        St = eyeData.stimulusStart; % where to start analysis
        St = St(ismember(St(:,2),validIndex),:);
        
        End = eyeData.choiceMade; % where to end analysis
%         End = eyeData.stimulusFin; % where to end analysis
        End = End(ismember(End(:,2),validIndex),:);
        
        if ~isequal(St(:,2),End(:,2))
            error('Please check, the trials for start points and end points are not the same.');
        end
        
        St = St(con,:);
        End = End(con,:);
        eyeDataIndex = nan(size(St,1), max(End(:,1)-St(:,1))+1);
        
        for i = 1:size(St,1)
            normSt = find(eyeData.eyeData(:,1) >= (St(i,1)-100),1,'first');
            normEd = find(eyeData.eyeData(:,1) >= St(i,1),1,'first');
            normalize = mean(eyeData.eyeData(normSt:normEd,4));
            
            stIndex = find(eyeData.eyeData(:,1) >= St(i,1),1,'first');
            enIndex = find(eyeData.eyeData(:,1) >= End(i,1),1,'first');
            eyeDataIndex(i, 1:End(i,1)-St(i,1)+1) = reshape(eyeData.eyeData(stIndex:enIndex,4)/normalize,1,[]);
        end
        
        % remove invalid data
        while sum(isnan(eyeDataIndex(:,end)),'all')
            eyeDataIndex(:,end) = [];
        end
        
        % calculate for mean and standard error
        errorBar = nanstd(eyeDataIndex,0,1)./sqrt(nansum(eyeDataIndex,1)); % standard error
        eyeMean = nanmean(eyeDataIndex,1);
        
        % draw plot
        if coni ~=5
            shadedErrorBar([],eyeMean,errorBar,'lineprops',{'color',color{coni}},'transparent',1);
        else
            shadedErrorBar([],eyeMean,errorBar,'lineprops',{'-*r'},'transparent',1);
        end
        
        % plot flash duration
        if ismember(coni,[1 2 5])
            plot([0 0],[min(eyeMean,[],'all') max(eyeMean,[],'all')],'-.r');
            plot([20 20],[min(eyeMean,[],'all') max(eyeMean,[],'all')],'--r');
        else
            plot([0 0],[min(eyeMean,[],'all') max(eyeMean,[],'all')],'-.r');
            plot([20 20],[min(eyeMean,[],'all') max(eyeMean,[],'all')],'--r');
            plot([110 110],[min(eyeMean,[],'all') max(eyeMean,[],'all')],'-.r');
            plot([130 130],[min(eyeMean,[],'all') max(eyeMean,[],'all')],'--r');
        end
        
        % plot beep duration
        if ismember(coni,[1,3])
            plot([0 0],[min(eyeMean,[],'all') max(eyeMean,[],'all')],'-.k');
            plot([20 20],[min(eyeMean,[],'all') max(eyeMean,[],'all')],'--k');
        else
            plot([0 0],[min(eyeMean,[],'all') max(eyeMean,[],'all')],'-.k');
            plot([20 20],[min(eyeMean,[],'all') max(eyeMean,[],'all')],'--k');
            plot([110 110],[min(eyeMean,[],'all') max(eyeMean,[],'all')],'-.k');
            plot([130 130],[min(eyeMean,[],'all') max(eyeMean,[],'all')],'--k');
        end
    end
end

%% functions below
function edf2ascConvert(dataPath,exePath,overWrite)
% This script can convert EYELINK DATA FILE(.EDF) to ASCII(.asc) files.
% The processing might need a period of time. If you have a significant
% number of files need to processing, please well arrange your time
% before running this script.
%
% overWrite: how to deal the existed asc files? true - overwrite; false - skip.
% default as true.
%
% Write by BYC June,2019
% Modified  by BYC Aug,2021

if nargin==2
    overWrite = true;
elseif nargin >3
    error('Invalid input number');
end

if ~exist('dataPath','var')
    error('There is no file in this path or there is a wrong path.');
end

if ~exist('exePath','var')
    error('There is no file in this path or there is a wrong path.');
elseif contains(exePath,'edf2asc.exe')
    exeFullPath = [exePath ' -ntime_check'];
else
    exeFullPath = [fullfile(exePath,'edf2asc.exe') ' -ntime_check'];
end
curPath = pwd;
cd(dataPath);
dataFile= dir(fullfile(dataPath, '*.EDF'));

for i = 1:length(dataFile)
    originFilePath = fullfile(dataPath,dataFile(i).name);
    
    if contains(originFilePath,'.EDF')
        ascName = strrep(dataFile(i).name,'.EDF','.asc');
    elseif contains(originFilePath,'.edf')
        ascName = strrep(dataFile(i).name,'.edf','.asc');
    else
        error([dataFile(i).name ' is not a EDF file']);
    end
    
    if exist(fullfile(dataPath,ascName),'file')
        if overWrite
            delete(fullfile(dataPath,ascName));
        else
            continue
        end
    end
    
    cmd = [exeFullPath 32 originFilePath];
    [~,log] = system(cmd);
    changeLine = strfind(log,char(13));
    disp(log(changeLine(end-1):end));
end
cd(curPath);
end

function asc2matConvert(dataPath,overWrite)
% This script can convert ASCII formed EYELINK DATA FILE (.asc) to matlab data file (.mat) files.
% This script need to be modified for different EYELINK MESSAGE MARKERS
% The processing might need a period of time. If you have a significant
% number of files need to processing, please well arrange your time
% before running this script.
%
% Write by BYC Aug,2021
if nargin==1
    overWrite = true;
elseif nargin >2
    error('Invalid input number');
end

curDir = pwd;
cd(dataPath);
ascFile = dir('*.asc');
for i = 1:length(ascFile)
    ascNamei = ascFile(i).name;
    matNamei = ['Converted_' strrep(ascNamei,'.asc','.mat')];
    
    if exist(matNamei,'file')
        if overWrite
            delete(matNamei);
        else
            continue
        end
    end
    
    fid = fopen(ascNamei);
    fseek(fid,0,'eof');
    numline = ftell(fid);
    fclose(fid);
    rawData = importdata(ascNamei,' ',numline);
    
    pIStr = rawData(contains(rawData,'point in','IgnoreCase',true));
    fixationPoints = nan(length(pIStr),3);
    for j = 1:length(pIStr)
        fixationPoints(j,:) = cell2mat(cellfun(@str2num,regexp(pIStr{j},'\d*\.?\d*','match'),'UniformOutput',false));
    end
%     delIndex = diff(fixationPoints(:,2))==0;
%     fixationPoints(delIndex,:) = [];
    
%     pPStr = rawData(contains(rawData,'pursuit','IgnoreCase',true));
%     pursuitPoints = nan(length(pPStr),2);
%     for j = 1:length(pPStr)
%         pursuitPoints(j,:) = cell2mat(cellfun(@str2num,regexp(pPStr{j},'\d*\.?\d*','match'),'UniformOutput',false));
%     end
%     delIndex = diff(pursuitPoints(:,2))==0;
%     pursuitPoints(delIndex,:) = [];
%     
%     sfStr = rawData(contains(rawData,'stimulus finish','IgnoreCase',true));
%     stimulusFin = nan(length(sfStr),2);
%     for j = 1:length(sfStr)
%         stimulusFin(j,:) = cell2mat(cellfun(@str2num,regexp(sfStr{j},'\d*\.?\d*','match'),'UniformOutput',false));
%     end
%     delIndex = diff(stimulusFin(:,2))==0;
%     stimulusFin(delIndex,:) = [];
%     
%     scStr = rawData(contains(rawData,'start choice','IgnoreCase',true));
%     startChoice = nan(length(scStr),2);
%     for j = 1:length(scStr)
%         startChoice(j,:) = cell2mat(cellfun(@str2num,regexp(scStr{j},'\d*\.?\d*','match'),'UniformOutput',false));
%     end
%     delIndex = diff(startChoice(:,2))==0;
%     startChoice(delIndex,:) = [];
%     
%     cmStr = rawData(contains(rawData,'choice made','IgnoreCase',true));
%     choiceMade = nan(length(cmStr),2);
%     for j = 1:length(cmStr)
%         choiceMade(j,:) = cell2mat(cellfun(@str2num,regexp(cmStr{j},'\d*\.?\d*','match'),'UniformOutput',false));
%     end
%     delIndex = diff(choiceMade(:,2))==0;
%     choiceMade(delIndex,:) = [];
    
%     % find valid trials
%     Ind1 = intersect(fixationPoints(:,2),pursuitPoints(:,2));
%     Ind2 = intersect(Ind1,stimulusFin(:,2));
%     Ind3 = intersect(Ind2,startChoice(:,2));
%     validIndex = intersect(Ind3,choiceMade(:,2));

%     % extract valid trials
%     fixationPoints = fixationPoints(ismember(fixationPoints(:,2),validIndex),:);
%     pursuitPoints = pursuitPoints(ismember(pursuitPoints(:,2),validIndex),:);
%     stimulusFin = stimulusFin(ismember(stimulusFin(:,2),validIndex),:);
%     startChoice = startChoice(ismember(startChoice(:,2),validIndex),:);
%     choiceMade = choiceMade(ismember(choiceMade(:,2),validIndex),:);
    
    eyeDataIndex = contains(rawData,'...');
    eyeData = rawData(eyeDataIndex);
    for linei = 1:length(eyeData)
        delIndex = strfind(eyeData{linei},'...');
        eyeData{linei} = eyeData{linei}(1:delIndex-1);
    end
    eyeData = strrep(eyeData,char(9),'');
    eyeData = cellfun(@str2num,eyeData,'UniformOutput',false);
    
    delIndex = cell2mat(cellfun(@isempty,eyeData,'UniformOutput',false));
    eyeData(delIndex) = [];
    
    % add nan value to blink and invalid data
    temp = cell2mat(eyeData(1:10));
    dt = mode(diff(temp(:,1)));
    addList = find(cellfun(@isempty,eyeData));
    for addi = reshape(addList,1,[])
        t = cell2mat(eyeData(addi-1));
        eyeData{addi} = [t(1)+dt,nan,nan,nan];
    end
    eyeData = cell2mat(eyeData);
    
    save(matNamei,'fixationPoints','eyeData');
end
cd(curDir)
end