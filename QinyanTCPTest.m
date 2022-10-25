% writeline() and readline() were published in R2020b
function QinyanTCPTest()
path = pwd;
testtimes = 500;
dataFrequency = 200;
resultFileName = ['TCPTestResult_' datestr(now,'yymmddHHMM') ];

%% tcp
t = tcpip('127.0.0.1',9998);
t.OutputBuffersize=100000;
fopen(t);
fwrite(t,double(['savepath:' fullfile(path,resultFileName)])); %set path
pause(0.5)
fwrite(t,'init'); % initial eye tracker
pause(1);

fwrite(t,'cal'); % calibration
for i = 1:5
    A = fread(t)
    if ~isempty(A)
        break
    elseif i == 5
        fwrite(t,'stopcal');
    end
end

fwrite(t,'val'); % validation
for i = 1:5
    B = fread(t)
    if ~isempty(B)
        break
    elseif i ==5
        fwrite(t,'stopval');
    end
end

fwrite(t,'singlecal'); % drift correction
for i = 1:3
    C = fread(t)
    if ~isempty(C)
        break
    elseif i ==3
        fwrite(t,'stopsinglecal');
    end
end
pause(0.5);
fwrite(t,'start');
pause(0.5);
fwrite(t,['getdata:' num2str(dataFrequency)]);
text = cell(testtimes,1);
time = nan(testtimes,1);
T = tic;
fwrite(t,'mark:test start');
for i = 1:testtimes
    str = fread(t);
%     str = fread(t,10); % set size as 10
%     char(str')
    text{i} = char(str');
    time(i) = toc(T);
end
fwrite(t,'mark:test finish');
fwrite(t,'stopgetdata');
pause(0.5);
fwrite(t,'stop');
pause(0.5);
fwrite(t,'close');
fclose(t);
save([resultFileName '.mat'],'text','time','testtimes','dataFrequency');
disp('TCP test finished')