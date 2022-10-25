function QinyanUDPTest()
path = pwd;
testtimes = 500;
dataFrequency = 200;
resultFileName = ['UDPTestResult_' datestr(now,'yymmddHHMM')];

%% udp
u=udp('127.0.0.1',9999);
fopen(u);
fwrite(u,['savepath:' fullfile(path,resultFileName) ]); %set path
pause(0.5)
fwrite(u,'init'); % initial eye tracker
pause(1);

fwrite(u,'cal'); % calibration
for i = 1:5
    A = fread(u,10)
    if ~isempty(A)
        break
    elseif i == 5
        fwrite(u,'stopcal');
    end
end

fwrite(u,'val'); % validation
for i = 1:5
    B = fread(u,10)
    if ~isempty(B)
        break
    elseif i ==5
        fwrite(u,'stopval');
    end
end

fwrite(u,'singlecal'); % drift correction
for i = 1:3
    C = fread(u,10)
    if ~isempty(C)
        break
    elseif i ==3
        fwrite(u,'stopsinglecal');
    end
end
pause(0.5);
fwrite(u,'start');
pause(0.5);
fwrite(u,['getdata:' num2str(dataFrequency)]);
text = cell(testtimes,1);
time = nan(testtimes,1);
T = tic;
fwrite(u,'mark:test start');
for i = 1:testtimes
    str = fread(u,10);
%     char(str')
    text{i} = char(str');
    time(i) = toc(T);
end
fwrite(u,'mark:test finish');
fwrite(u,'stopgetdata');
pause(0.5);
fwrite(u,'stop');
pause(0.5);
fwrite(u,'close');
fclose(u);
save([resultFileName '.mat'],'text','time','testtimes','dataFrequency');
disp('UDP test finished')