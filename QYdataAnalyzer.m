function QYdataAnalyzer(path)
if nargin < 1
    path = fullfile(pwd,'data');
end

% QYdataExtract(path);

files = dir(fullfile(path,'extracted_*.mat'));
frameD = cell(length(files),1);

%% set config for subplots
set(0,'defaultfigurecolor','w');
set(figure(1),'pos',[20 20 1849 900],'Name','End-to-End sample delay');clf;
subplot1 = tight_subplot(4,length(files),[0.1 0.035],[0.1 0.1]);

for filei = 1:length(files)
    load(fullfile(files(filei).folder, files(filei).name));
    param = load(fullfile(files(filei).folder, strrep(files(filei).name,'extracted_','')));
    
    frameD{filei,1} = cell(size(fixData,1),1);
    frameD{filei,2} = cell(size(fixData,1),1);
    for i = 1:size(fixData,1)
        t=cell2mat(fixData{i,1}(:,1));
        eyeData = cell2mat(fixData{i,1}(:,2));
        frameD{filei,1}{i,1} = round(diff(t)/datenum(milliseconds(1)),2);
        a = [1 3 5];b=[2 4 6];
        frameD{filei,2}{i,1} = eyeData(:,1:6);
        frameD{filei,2}{i,1}(:,a) = eyeData(:,a)*param.SCREEN.widthPix/2+param.SCREEN.widthPix/2;
        frameD{filei,2}{i,1}(:,b) = param.SCREEN.heightPix-(-eyeData(:,b)*param.SCREEN.heightPix/2+param.SCREEN.heightPix/2);
        frameD{filei,2}{i,1}(:,7:8) = repmat(([0, param.SCREEN.heightPix] - fixData{i,2}).*[-1,1],size(eyeData,1),1);
    end
    
    axes(subplot1(filei));
    hold on
    sampleDelay = cell2mat(frameD{filei,1}(:,1));
    plot(sampleDelay,'k','LineWidth',2);
    
    lengthNum = cellfun(@length,frameD{filei,1});
    lengthPos = cumsum(lengthNum);
    
    for i = 1:length(lengthPos)
        plot([lengthPos(i),lengthPos(i)], [min(sampleDelay), max(sampleDelay)], '--k');
    end
    
    textStr = ['mean ¡À sd = ' num2str(mean(sampleDelay)) ' ¡À ' num2str(std(sampleDelay))];
    text(20,33,textStr);
%     title(['End-to-End sample delay for QY-I in sample ' num2str(filei) ' (ms)']);
    xticklabels('auto');
    yticklabels('auto');
    ylim([0 35]);    
    
    axes(subplot1(filei+length(files)));
    hold on
    eyePos = cell2mat(frameD{filei,2});
    eyePos(sum(eyePos<0,2)>0,:) = nan; % eliminate blinks as nan
    plot(eyePos(:,1),eyePos(:,2),'k');
    plot(eyePos(:,3),eyePos(:,4),'r');
    plot(eyePos(:,5),eyePos(:,6),'g');
    xticklabels('auto');
    yticklabels('auto');
%     ylim([-1.2 1.2]);

    axes(subplot1(filei));    
    sampleDelay(sum(eyePos<0,2)>0,:) = nan; % eliminate blinks as nan
    plot(sampleDelay,'y','LineWidth',1);
    
    axes(subplot1(filei+length(files)*2));
    hold on
    meanEye = atand(sqrt(((eyePos(:,1)-eyePos(:,7))/param.SCREEN.widthPix*param.SCREEN.widthCM).^2 + ...
        ((eyePos(:,2)-eyePos(:,8))/param.SCREEN.heightPix*param.SCREEN.heightCM).^2) / param.SCREEN.distance);
    lEye = atand(sqrt(((eyePos(:,3)-eyePos(:,7))/param.SCREEN.widthPix*param.SCREEN.widthCM).^2 + ...
        ((eyePos(:,4)-eyePos(:,8))/param.SCREEN.heightPix*param.SCREEN.heightCM).^2) / param.SCREEN.distance);
    rEye = atand(sqrt(((eyePos(:,5)-eyePos(:,7))/param.SCREEN.widthPix*param.SCREEN.widthCM).^2 + ...
        ((eyePos(:,6)-eyePos(:,8))/param.SCREEN.heightPix*param.SCREEN.heightCM).^2) / param.SCREEN.distance);
    
    for i = 1:length(lengthPos)
        plot([lengthPos(i),lengthPos(i)], [min(sampleDelay), max(sampleDelay)], '--k');
    end
    
    plot(meanEye,'k');
    plot(lEye,'r');
    plot(rEye,'g');
    xticklabels('auto');
    yticklabels('auto');
    textStr = {'Validation result: ', num2str(validResult{end,end})};
    text(20,23,textStr,'FontSize',7);
    ylim([0,20])
    
    axes(subplot1(filei+length(files)*3));
    hold on
    
    for i = 1:length(lengthPos)
        plot([lengthPos(i),lengthPos(i)], [min(sampleDelay), max(sampleDelay)], '--k');
    end
    
    plot(meanEye,'k');
    plot(lEye,'r');
    plot(rEye,'g');
    xticklabels('auto');
    yticklabels('auto');
    textStr = {['Median for two-eye average: ' num2str(nanmedian(meanEye))], ['Median for left eye: ' num2str(nanmedian(lEye))], ['Median for right eye: ' num2str(nanmedian(rEye))]};
    text(23,2.5,textStr,'FontSize',9);
    ylim([0 2]);
end


end

function [ha, pos] = tight_subplot(Nh, Nw, gap, marg_h, marg_w)
% tight_subplot creates "subplot" axes with adjustable gaps and margins
%
% [ha, pos] = tight_subplot(Nh, Nw, gap, marg_h, marg_w)
%
%   in:  Nh      number of axes in hight (vertical direction)
%        Nw      number of axes in width (horizontaldirection)
%        gap     gaps between the axes in normalized units (0...1)
%                   or [gap_h gap_w] for different gaps in height and width 
%        marg_h  margins in height in normalized units (0...1)
%                   or [lower upper] for different lower and upper margins 
%        marg_w  margins in width in normalized units (0...1)
%                   or [left right] for different left and right margins 
%
%  out:  ha     array of handles of the axes objects
%                   starting from upper left corner, going row-wise as in
%                   subplot
%        pos    positions of the axes objects
%
%  Example: ha = tight_subplot(3,2,[.01 .03],[.1 .01],[.01 .01])
%           for ii = 1:6; axes(ha(ii)); plot(randn(10,ii)); end
%           set(ha(1:4),'XTickLabel',''); set(ha,'YTickLabel','')

% Pekka Kumpulainen 21.5.2012   @tut.fi
% Tampere University of Technology / Automation Science and Engineering


if nargin<3; gap = .02; end
if nargin<4 || isempty(marg_h); marg_h = .05; end
if nargin<5; marg_w = .05; end

if numel(gap)==1
    gap = [gap gap];
end
if numel(marg_w)==1
    marg_w = [marg_w marg_w];
end
if numel(marg_h)==1
    marg_h = [marg_h marg_h];
end

axh = (1-sum(marg_h)-(Nh-1)*gap(1))/Nh; 
axw = (1-sum(marg_w)-(Nw-1)*gap(2))/Nw;

py = 1-marg_h(2)-axh; 

% ha = zeros(Nh*Nw,1);
ii = 0;
for ih = 1:Nh
    px = marg_w(1);
    
    for ix = 1:Nw
        ii = ii+1;
        ha(ii) = axes('Units','normalized', ...
            'Position',[px py axw axh], ...
            'XTickLabel','', ...
            'YTickLabel','');
        px = px+axw+gap(2);
    end
    py = py-axh-gap(1);
end
if nargout > 1
    pos = get(ha,'Position');
end
ha = ha(:);
end

function QYdataExtract(path)
if nargin < 1
    path = fullfile(pwd,'data');
end

files = dir(fullfile(path,'*.txt'));

for filei = 1:length(files)
    fileSaveName = strrep(files(filei).name, '.txt', '.mat');
    fileH = fopen(fullfile(files(filei).folder, files(filei).name),'r');
    dateStemp = strfind(files(filei).name,'_');
    dateNum = fileSaveName(dateStemp(2)+1:dateStemp(2)+6);
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
        elseif pointNum>0
            dotLoc = strfind(lineData,',');
            dotLoc = dotLoc(1);
            timeStemp = datenum(datetime([dateNum,'-', lineData(1:dotLoc-1)],'InputFormat','yyMMdd-HH:mm:ss.SSS'));
            eyeData = str2num(lineData(dotLoc+1:end));
            fixData{pointNum,1} = cat(1,fixData{pointNum,1},{timeStemp, eyeData});
        end
    end
    fclose(fileH);
    save(fullfile(files(filei).folder, ['extracted_' fileSaveName]), 'fixData','validResult');
end
end