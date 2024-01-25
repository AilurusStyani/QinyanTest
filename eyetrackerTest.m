function eyetrackerTest(eyetracker,task,distance)
% essential environment: psychtoolbox-3
delete TRIALINFO SCREEN
global SCREEN
global TRIALINFO

if nargin <1
    TRIALINFO.eyetracker = 0; % 0 for no TRIALINFO.eyetracker mode; 1 for eyelink; 2 for QY-1
else
    TRIALINFO.eyetracker = eyetracker;
end

if nargin <2
    task = 'fixation'; % fixation: one point fixation task without distractor;
                                  % pursuit: linear pursuit from left to right and top to down
                                  % scan: points from the mid cycle to the border
end

if nargin<3
    SCREEN.distance = 60; % subject to screen distance, unit cm
else
    SCREEN.distance = distance; % cm
end

% fixaiton
TRIALINFO.fixationDuration = 10; % second
TRIALINFO.fixationSizeD = 1.5; % degree
TRIALINFO.fixationCenterSizeD = 0.3; % pixel

% pursuit
TRIALINFO.LRDegree = 20; % degree
TRIALINFO.TDDegree = 15; % degree
TRIALINFO.pursuitDuration = 2; % second

% scan
TRIALINFO.scanXPointNum = 5;
TRIALINFO.scanYPointNum = 3;
TRIALINFO.scanFixDuration = 1.5; % second

saveDir = fullfile(pwd,'data');
mkdir(saveDir);
curdir = pwd;

subjectName = inputdlg({'Please input participant''s initials.'},'Subject Name',1,{''},'on');
if isempty(subjectName)
    return
end
subjectName = cell2mat(subjectName);

if TRIALINFO.eyetracker == 0
    resultFileName = ['DemoTestResult_' subjectName '_' datestr(now,'yymmddHHMM')];
elseif TRIALINFO.eyetracker == 1
    resultFileName = ['EyelinkTestResult_' subjectName '_' datestr(now,'yymmddHHMM')];
elseif TRIALINFO.eyetracker == 2
    resultFileName = ['QYITestResult_' subjectName '_' datestr(now,'yymmddHHMM')];
    QY.testtimes = 500;
    QY.dataFrequency = 200;
end

% set keyboard
KbName('UnifyKeyNames');
escape    = KbName('ESCAPE');
enter     = KbName('Return');
spaceKey = KbName('space');

upArror   = KbName('UpArrow');
leftKey   = KbName('LeftArrow');
rightKey  = KbName('RightArrow');

cKey = KbName('c');
vKey = KbName('v');
dKey = KbName('d');

escFlag = false;

%% initial QY-I
if TRIALINFO.eyetracker == 2
    % udp
    u=udp('127.0.0.1',9999);
    fopen(u);
    fwrite(u,['savepath:' fullfile(saveDir,[resultFileName '.txt']) ]); %set path
    pause(0.5)
    fwrite(u,'init'); % initial eye tracker
    pause(1);
    caseId = [];
    readData = false;
    escFlag = false;
    % calibration and validation
    while true
        [keyDown, ~, keyCode]=KbCheck;
        if keyCode(escape)
            escFlag = true;
            break
        elseif keyCode(cKey) % c to start calibration
            caseId = 'cal';
        elseif keyCode(vKey) % v to start validation
            caseId = 'val';
            readData = true;
        elseif keyCode(dKey) % d to start drift correction
            caseId = 'singlecal';
        elseif keyCode(spaceKey) % skip calibration and validation
            break
        end
        
        if keyDown
            fwrite(u,caseId);
            caseId = [];
        end
        
        if readData
            returnStr = fread(u,10);
            if ~isempty(returnStr)
                strTurn = char(returnStr');
                disp(strTurn);
                if contains(strTurn,'get validate values')
                    break
                end
            end
        end

        while keyDown % check unless key release
            [keyDown, ~, ~]=KbCheck;
            pause(0.1);
        end
    end
    
    if escFlag
        fwrite(u,'stop');
        pause(0.5);
        fwrite(u,'close');
        fclose(u);
        return
    end
    
    pause(1);
    
    disp('Is everything going smoothly? Or press ESC to terminate.')
    tic
    while toc<3
        [~, ~, keyCode]=KbCheck;
        if keyCode(escape)
            return
        end
        pause(0.3)
    end
    
    fwrite(u,'start');
    pause(0.5);
    fwrite(u,['getdata:' num2str(QY.dataFrequency)]);
    pause(0.2)
    fwrite(u,['mark:test start ' task]);
end


%% initial PTB
Screen('Preference', 'SkipSyncTests', 2);
AssertOpenGL;
oldVisualDebugLevel = Screen('Preference','VisualDebugLevel',3);
oldSupressAllWarnings = Screen('Preference','SuppressAllWarnings',1);

if max(Screen('Screens')) > 1
    SCREEN.screenId = max(Screen('Screens'))-1;
else
    SCREEN.screenId = max(Screen('Screens'));
end
PsychImaging('PrepareConfiguration');

% Define background color:
whiteBackground = WhiteIndex(SCREEN.screenId);
blackBackground = BlackIndex(SCREEN.screenId);

% Open a double-buffered full-screen window on the main displays screen.
[win , winRect] = PsychImaging('OpenWindow', SCREEN.screenId, blackBackground);
SCREEN.widthPix = winRect(3);
SCREEN.heightPix = winRect(4);
SCREEN.center = [SCREEN.widthPix/2, SCREEN.heightPix/2];

[width, height] = Screen('DisplaySize', SCREEN.screenId);
SCREEN.widthCM = width/10; % mm to cm
SCREEN.heightCM = height/10; % mm to cm

TRIALINFO.fixationSizeP = degree2pix(TRIALINFO.fixationSizeD/2);
TRIALINFO.fixationCenterSizeP = degree2pix(TRIALINFO.fixationCenterSizeD/2);

TRIALINFO.fixationPosition = [SCREEN.widthPix/2,SCREEN.heightPix/2];

SCREEN.refreshRate = Screen('NominalFrameRate', SCREEN.screenId);

Screen('FillRect', win ,blackBackground,[0 0 SCREEN.widthPix SCREEN.heightPix]);

%% initial eyelink
if TRIALINFO.eyetracker == 1
    tempName = 'TEMP1'; % need temp name because Eyelink only know hows to save names with 8 chars or less. Will change name using matlab's moveFile later.
    dummymode=0;
    
    el=EyelinkInitDefaults(win);
    %     el.backgroundcolour = BlackIndex(el.window);
    %     el.foregroundcolour = GrayIndex(el.window);
    %     el.msgfontcolour    = WhiteIndex(el.window);
    %     el.imgtitlecolour   = WhiteIndex(el.window);
    el.calibrationtargetsize=1;  % size of calibration target as percentage of screen
    el.calibrationtargetwidth=0.5; % width of calibration target's border as percentage of screen
    
    if ~EyelinkInit(dummymode)
        fprintf('Eyelink Init aborted.\n');
        cleanup;  % cleanup function
        Eyelink('ShutDown');
        Screen('CloseAll');
        return
    end
    
    testi = Eyelink('Openfile', tempName);
    if testi~=0
        fprintf('Cannot create EDF file.\n');
        cleanup;
        Eyelink('ShutDown');
        Screen('CloseAll');
        return
    end
    
    %   SET UP TRACKER CONFIGURATION
    Eyelink('command', 'calibration_type = HV9');
    %	set parser (conservative saccade thresholds)
    Eyelink('command', 'saccade_velocity_threshold = 35');
    Eyelink('command', 'saccade_acceleration_threshold = 9500');
    Eyelink('command', 'link_event_filter = LEFT,RIGHT,FIXATION,SACCADE,BLINK,MESSAGE,BUTTON,FIXUPDATE,INPUT');
    Eyelink('command', 'link_sample_data  = LEFT,RIGHT,GAZE,HREF,GAZERES,AREA,STATUS,INPUT,HTARGET');
    Eyelink('command', 'online_dcorr_refposn = %1d, %1d', SCREEN.center(1), SCREEN.center(2));
    Eyelink('command', 'online_dcorr_maxangle = %1d', 30.0);
    
    % you must call this function to apply the changes from above
    EyelinkUpdateDefaults(el);
    
    % Calibrate the eye tracker
    EyelinkDoTrackerSetup(el);
    
    % do a final check of calibration using driftcorrection
%     EyelinkDoDriftCorrection(el);
    
    Eyelink('StartRecording');
    
    Eyelink('message', 'SYNCTIME');	 	 % zero-plot time for EDFVIEW
    eye_used = Eyelink('EyeAvailable'); % get eye that's tracked
    if eye_used ~= el.BINOCULAR % if both eyes are tracked
        if eye_used == el.LEFTEYE % use left eye
            fprintf('Left eye is in recording.\n');
        elseif eye_used == el.RIGHTEYE % use left eye
            fprintf('Right eye is in recording.\n');
        end
    else
        fprintf('Both eyes are in recording.\n');
    end
    
    errorCheck=Eyelink('checkrecording'); 		% Check recording status */
    if(errorCheck~=0)
        fprintf('Eyelink checked wrong status.\n');
        cleanup;  % cleanup function
        Eyelink('ShutDown');
        Screen('CloseAll');
    end
    
    [keyDown, ~, ~]=KbCheck;
    while keyDown % check unless key release
        [keyDown, ~, ~]=KbCheck;
        pause(0.1);
    end
%     pause(1); % wait a little bit, in case the key press during calibration influence the following keyboard check
end

%% tasks
% fixation
if contains(task,'fixation')
    Screen('Flip', win);
    drawPoint(TRIALINFO.fixationPosition,TRIALINFO.fixationSizeP, TRIALINFO.fixationCenterSizeP,win);
    Screen('Flip', win);
    if TRIALINFO.eyetracker == 1
        Eyelink('message',['point in ' num2str(TRIALINFO.fixationPosition)]);
    elseif TRIALINFO.eyetracker == 2
        fwrite(u,['mark: point in ' num2str(TRIALINFO.fixationPosition)]);
    end
    pause(TRIALINFO.fixationDuration);
    Screen('Flip', win);
end

% pursuit
if contains(task,'pursuit')
    Screen('Flip', win);

    lrSTP = TRIALINFO.fixationPosition - [degree2pix(TRIALINFO.LRDegree)/2,0];
    lrEDP = TRIALINFO.fixationPosition + [degree2pix(TRIALINFO.LRDegree)/2,0];
    tdSTP = TRIALINFO.fixationPosition + [0,degree2pix(TRIALINFO.TDDegree)/2];
    tdEDP = TRIALINFO.fixationPosition - [0,degree2pix(TRIALINFO.TDDegree)/2];
    
    frameNum = round(TRIALINFO.pursuitDuration * SCREEN.refreshRate);
    lr = linspace(SCREEN.widthPix/2-degree2pix(TRIALINFO.LRDegree)/2, SCREEN.widthPix/2+degree2pix(TRIALINFO.LRDegree)/2, frameNum);
    td = linspace(SCREEN.heightPix/2-degree2pix(TRIALINFO.TDDegree)/2, SCREEN.heightPix/2+degree2pix(TRIALINFO.TDDegree)/2, frameNum);
    
    lrPos = [lr', ones(frameNum,1)*SCREEN.heightPix/2];
    tdPos = [ones(frameNum,1)*SCREEN.widthPix/2, td'];
    
    % left-right pursuit
    for i = 1:frameNum
        position = lrPos(i,:);
        drawPoint( position,TRIALINFO.fixationSizeP, TRIALINFO.fixationCenterSizeP,win);
        Screen('Flip', win);
        if TRIALINFO.eyetracker == 1
            Eyelink('message',['point in ' num2str(position)]);
        elseif TRIALINFO.eyetracker == 2
            fwrite(u,['mark: point in ' num2str(position)]);
        end
        if i == 1
            pause(2);
        end
    end
    pause(2);
    
    % top-down pursuit
    for i = 1:frameNum
        position = tdPos(i,:);
        drawPoint( position,TRIALINFO.fixationSizeP, TRIALINFO.fixationCenterSizeP,win);
        Screen('Flip', win);
        if TRIALINFO.eyetracker == 1
            Eyelink('message',['point in ' num2str(position)]);
        elseif TRIALINFO.eyetracker == 2
            fwrite(u,['mark: point in ' num2str(position)]);
        end
        if i == 1
            pause(2);
        end
    end
    pause(2);
end

% scan
if contains(task,'scan')
    if TRIALINFO.eyetracker == 2
        ScanFunction(win,u);
    else
        ScanFunction(win);
    end
end

%% terminate
if TRIALINFO.eyetracker == 1 % terminate eyelink
    Eyelink('StopRecording');
    Eyelink('CloseFile');
    try
        fprintf('Receiving data file ''%s''\n',resultFileName);
        status=Eyelink('ReceiveFile',tempName ,saveDir,1);
        if status > 0
            fprintf('ReceiveFile status %d\n ', status);
        end
        if exist(resultFileName, 'file')==2
            fprintf('Data file ''%s'' can be found in '' %s\n',resultFileName, saveDir);
        end
    catch
        fprintf('Problem receiving data file ''%s''\n',resultFileName);
    end
    movefile(fullfile(saveDir,[tempName '.edf']),fullfile(saveDir,[resultFileName '.edf']));
    pause(2);
    Eyelink('ShutDown');
elseif TRIALINFO.eyetracker == 2 % terminate QY-I
    fwrite(u,'mark:test finished');
    fwrite(u,'stopgetdata');
    pause(0.5);
    fwrite(u,'stop');
    pause(0.5);
    fwrite(u,'close');
    fclose(u);
end
Screen('CloseAll');
save(fullfile(saveDir,resultFileName));
end
%% functions
function pixel = degree2pix(degree,dir)
% this function convert degree value to pixel value
% it is better been used to calculte the pixel length from the central point
% On X axis: dir = 1; On Y axis: dir = 2
global SCREEN

if nargin == 1
    a=SCREEN.widthPix / SCREEN.widthCM;
    b=SCREEN.heightPix / SCREEN.heightCM;
    if abs(a-b)/min(a,b) < 0.05
        length = tand(degree) * SCREEN.distance;
        pixel = length / SCREEN.widthCM * SCREEN.widthPix;
    else
        error('Error in screen parameter or screen config, or you should define it is horizontal(1) / vertical(2).')
    end
elseif nargin == 2
    length = tand(degree) * SCREEN.distance;
    if dir == 1
        pixel = length / SCREEN.widthCM * SCREEN.widthPix;
    elseif dir == 2
        pixel = length / SCREEN.heightCM * SCREEN.heightPix;
    else
        error('Invalid value for dir.')
    end
end
end

function drawPoint(position,sizeP,centerSizeP,win)
color1 = [255 255 255];
color2 = [255 0 0];
pointMetrix1 = [position(1)-sizeP/2, position(2)-sizeP/2, position(1)+sizeP/2, position(2)+sizeP/2];
pointMetrix2 = [position(1)-centerSizeP/2, position(2)-centerSizeP/2, position(1)+centerSizeP/2, position(2)+centerSizeP/2];
Screen('FillOval', win, color1, pointMetrix1);
Screen('FillOval', win, color2, pointMetrix2);
end

function ScanFunction(win,u)
global TRIALINFO SCREEN
if nargin
end

Screen('Flip', win);

% calculate for points
xMetrix = linspace(0,SCREEN.widthPix,TRIALINFO.scanXPointNum+2);
yMetrix = linspace(0,SCREEN.heightPix,TRIALINFO.scanYPointNum+2);
if or(isempty(xMetrix),isempty(yMetrix))
    warning('Invalid number for scanXPointNum or scanYPointNum, they must be integer numbers.')
    return
end

[X,Y] = meshgrid(xMetrix(2:end-1),yMetrix(2:end-1));
X = reshape(X, 1, numel(X));
Y = reshape(Y, 1, numel(Y));

for i = 1:length(X)
    position = [X(i), Y(i)];
    drawPoint( position,TRIALINFO.fixationSizeP, TRIALINFO.fixationCenterSizeP,win);
    Screen('Flip', win);
    if TRIALINFO.eyetracker == 1
        Eyelink('message',['point in ' num2str(position)]);
    elseif TRIALINFO.eyetracker == 2
        fwrite(u,['mark: point in ' num2str(position)]);
    end
    pause(TRIALINFO.scanFixDuration)
end
end