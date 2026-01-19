%% Active torque analysis from maximal voluntary contractions (MVCs)
% Run getPasTorque before this script!
% Initialise workspace
close all; clearvars -except files

%% Required functions
% filtTAE vline

%% Inputs
% Note: 0 = no; 1 = yes
reanalyse = 0;
pasFitFile = 'TQangfit.mat';
angScale = 0;   % change scaling of angle data
emgCh = 1;      % EMG channel to use if there are multiple
fps = 2000;     % samples per second
zz = 1;         % muscle compartment to analyse: 1=superficial; 2=deep
check = 0;      % check individual trials before continuing
zci = @(v) find(v(:).*circshift(v(:), [-1 0]) <= 0); % returns zero-crossing indices
winDur = 500;   % 500 ms overlapping window for EMG RMS amplitude calculation

%% Select files to analyse
if reanalyse == 0 && isfile('TQmvc.mat')
    load('TQmvc.mat','files')
end

if ~exist('files','var') || ischar(files) || ~iscell(files) && files == 0
    clear files
    [files(:,1), pname] = uigetfile('*.mat','MultiSelect','On');
end

%% Catch if one file is selected
if ~iscell(files)
    error('Error: Select files to analyse.')
end

if ~exist('mvc','var')
    %% Load necessary files to calculate active torque and active force
    % Load passive torque-angle fit:
    % P refers to steady state torque;
    % P1 to torque during passive lengthening;
    % P2 to torque during passive shortening
    load(pasFitFile,'P'); %TQangFit
    % Load TA moment arm-angle fit derived from Maganaris et al. 1999 data
    %load('C:\Users\brent\Documents\MATLAB\GitHub\DFG\maganaris1999_new.mat','MArest');
    %load('C:\Users\brent\Documents\MATLAB\GitHub\DFG\maganaris1999_new.mat','MAmvc');
    %MA = MAmvc; clear MAmvc

    %% Loop through files
    for jj = 1:length(files)
        if isfile(['T',files{jj}])
            [tq,ang,emg,us,mot,time,fname] = filtTAE(fps,['T',files{jj}]); % load combined fascicle data
            load(['T',files{jj}],'Fdat','Time') %timeUS
        else
            [tq,ang,emg,us,mot,time,fname] = filtTAE(fps,files{jj}); % load Spike2 data only
        end

        % Scaling of angle data was incorrect
        if angScale
            ang = ang/2;
        end

        %% Filter the rectified EMG signal with a dual-pass 2nd-order filter
        [emgRMS,~] = rmsDC(emg(:,emgCh), winDur, 1, fps);
        %emgAbsFilt = Wfilt(abs(emg(:,emgCh)),10,'low',fps);

        %% Calculate active torque
        tqPassive = polyval(P,ang);
        tqActive = tq-tqPassive;

        %% Calculate values of interest
        [mvc(jj,1),tqMaxX] = max(tqActive); % maximum torque value
        [mvc(jj,2),emgMaxX] = max(emgRMS); % maximum  RMS amplitude
        mvc(jj,3) = ang(tqMaxX); % angle at maximum torque

        if exist('Fdat','var')

            % Load tracked fascicle data
            FLx(:,1) = Fdat.Region(zz).FL;
            FAx(:,1) = rad2deg(Fdat.Region(zz).PEN);

            if ~exist('Time','var') && exist('timeUS','var')
                Time = timeUS;
            end

            % Upsample fascicle data rate to 2000 Hz
            if length(FLx) > length(Time)
                FLx(end) = [];
                FAx(end) = [];
            elseif length(FLx) < length(Time)
                Time(end) = [];
            end

            FL = interp1(Time,FLx,time,'linear');
            FA = interp1(Time,FAx,time,'linear');

            %% Calculate active fascicle force & length
            % Determine literature-based moment arms based on crank-arm
            % angle (0 deg = footplate perpendicular to tibia)
            %MAang = polyval(MA,ang)/100; % cm to m

            % Determine active tendon force
            %Ftendon = tqActive./MAang*0.5;

            % Determine active fascicle force
            %Ffas = Ftendon./cosd(FA);
            %[mvc(jj,4),fMaxX] = max(Ffas); % maximum fasicle force

            % Determine active fascicle length
            %mvc(jj,5) = FL(fMaxX);
            %mvc(jj,6) = min(FLx);
        end

        %% Plot MVC trials
        if check
            h = figure(jj); set(h,'Name',files{jj});
            subplot(311); hold on; plot(tq,'-'); plot(tqPassive,':'); plot(tqActive,'.');
            vline(tqMaxX,'k'); vline(emgMaxX,'r');
            subplot(312); hold on; plot(emg);
            vline(tqMaxX,'k'); vline(emgMaxX,'r');
            subplot(313); plot(ang);
            vline(tqMaxX,'k'); vline(emgMaxX,'r');

            waitfor(h)
        end

        clearvars -except reanalyse pasFitFile angScale emgCh fps zz check zci files P* MA mvc*
    end

    %% Select maximum values
    figure;
    hPts = plot(mvc(:,3),mvc(:,1),'o');
    brush on
    title('Select maximum values by highlighting them and clicking confirm')
    hPos = get(gcf);
    h = uicontrol('String', 'Confirm', 'Position', [hPos.OuterPosition(1,3)-120 10 100 50], ...
        'Callback', 'uiresume(gcbf)');
    uiwait(gcf);
    ptsSelected = logical(hPts.BrushData.');
    MVC = [hPts.XData(ptsSelected).' ...
        hPts.YData(ptsSelected).'];
    names = {'use'};
    assignin('base',names{1},MVC)
    close(gcf);

end

%% Plots of interest
figure; subplot(121); plot(mvc(:,3),mvc(:,2),'o');
xlabel('Angle (deg)'); ylabel('Maximum EMG RMS amplitude (V)')
subplot(122); plot(mvc(:,2),mvc(:,1),'o');
xlabel('Maximum EMG RMS amplitude (V)'); ylabel('Maximum active torque (Nm)')

% if size(mvc,2) > 3
% 
%     figure; hold on;
%     plot(mvc(:,5),mvc(:,4),'o'); plot(mvc(:,6),mvc(:,4),'rx');
%     xlabel('Active fascicle length (mm)'); ylabel('Maximum active fascicle force (N)')
% 
% end

%% Save analysed data
if isfile('TQmvc.mat') && reanalyse == 1
    ques = input('Overwrite previously saved MVC data y/n: ','s');
    if strcmpi(ques,'n')
        disp('TQmvc not overwritten.')
    elseif isempty(ques) || strcmpi(ques,'y')
        disp('TQmvc overwritten.')
        save('TQmvc.mat','files','mvc','MVC')
    end
else
    disp('TQmvc saved in current path.')
    save('TQmvc.mat','files','mvc','MVC')
end