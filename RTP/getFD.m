%% Active torque analysis from dynamic submaximal voluntary contractions
% Run getPasTorque then getMVC before this script!
% Initialise workspace
close all; clearvars -except files

%% Required functions / files
% filtTAE vline Wvel vline

%% Inputs
pasFitFile = 'TQfitSS.mat';
mvcFile = 'TQmax.mat';
emgCh = 1;
fps = 2000;     % samples per second
ssStart = 11;    % start of steady state in sec
dur = 13;       % contraction duration - 15 sec is typical
ht = 181;       % participant height
thrPL = 10;     % PL threshold
zz = 1;         % muscle compartment to analyse: 1=superficial; 2=deep
check = 1;      % check individual trials before continuing
workCheck = 0;  % check work calculation
zci = @(v) find(v(:).*circshift(v(:), [-1 0]) <= 0); % returns zero-crossing indices

ssStart = ssStart*fps+1;
dur = dur*fps;
cnt = [0 0]; % initialise counters

%% Select files
if ~exist('files','var') || ischar(files) || ~iscell(files) && files == 0
    clear files
    [files(:,1), pname] = uigetfile('*.mat','MultiSelect','On');
end

%% Catch if one file is selected
if ~iscell(files)
    error('Error: Select files to analyse.')
end

%% Remove files that aren't of interest
for kk = 1:length(files)
    % Check if the filename contains the specified string
    if contains(files{kk}, 'TQ') || contains(files{kk}, 'mvc') ||...
            contains(files{kk}, 'passive') || contains(files{kk}, 'X') 
            %contains(files{kk}, '--5')
        keep(kk,1) = 0;
    else
        keep(kk,1) = 1;
    end
end
% Remove files with specified string from above
files(keep==0) = [];

%% %% Load necessary files to calculate active torque and active force
% Passive torque fit
load(pasFitFile,'P'); %TQangFit
% TA moment arm data
load('C:\Users\Brent\Desktop\Research\Masters\Nogol\maganaris1999_new.mat','MArest');
MA = MArest; clear MArest
% MVC data
load(mvcFile,'mvc'); %TQmvc
% Find max MVC at 35 deg
mvcFiles = find(mvc(:,3)>=32.5);
[mvcMax,mvcMaxIdx] = max(mvc(mvcFiles,1));
%mvcAngle = mvc(mvcMaxIdx,3);

for jj = 1:length(files)
    if isfile(['T',files{jj}])
        [tq,ang,emg,us,mot,time,fname] = filtTAE(fps,['T',files{jj}]); % load combined fascicle data
        load(['T',files{jj}],'Fdat','Time')
        if ~exist('Time','var')
            load(['T',files{jj}],'timeUS')
            Time = timeUS;
        end
    else
        [tq,ang,emg,us,mot,time,fname] = filtTAE(fps,files{jj}); % load Spike2 data only
    end

    %% Filter the rectified EMG signal with a dual-pass 2nd-order filter
    emg = emg(:,emgCh);
    if istable(emg)
        emg = table2array(emg);
    end
    %emgAbsFilt = Wfilt(abs(emg),10,'low',fps);
    emgAbsFilt = rmsDC(emg,250,1,fps);

    %% Calculate active torque
    tqPassive = polyval(P,ang);
    tqActive = tq-tqPassive;

    %% Determine the type of trial
    angVel = Wvel(ang,fps);
    idx = zci(angVel);
    [~,maxDiff] = max(diff(idx));
    angStart = idx(maxDiff);
    angEnd = idx(maxDiff+1);
    angVelMedian = median(angVel(angStart:angEnd));
    angDis = ang(angEnd)-ang(angStart);

    % If there was an angular rotation, determine the preload used
    if abs(angDis) < 10
        PL = NaN;
    else
        % Rotation started at index 8001 (time: 4 s)
        PL = tqActive(angStart);
        %PLang = mean(ang(1:angStart));
        %[~,MVCang] = min(abs(mvc(:,3)-PLang));
        PL = PL/mvcMax*100;
    end

    %% Store desired data
    % Filename
    file{1,jj} = files{jj};

    % Condition
    cond(1,jj) = round(angDis,0);
    cond(2,jj) = round(angVelMedian,0);
    cond(3,jj) = round(ang(angStart),0);
    cond(4,jj) = round(ang(angEnd),0);
    cond(5,jj) = round(PL,0);

    % Active torque, angle, EMG
    TQ(:,jj) = tqActive(1:dur);
    ANG(:,jj) = ang(1:dur);
    EMG(:,jj) = emgAbsFilt(1:dur);

    % Calculate mean torque over steady state
    SS_TQ(1,jj) = mean(tqActive(ssStart:ssStart+fps-1));
    SS_TQ(2,jj) = mean(tqActive(ssStart+fps:ssStart+fps*2-1));
    SS_EMG(1,jj) = mean(emgAbsFilt(ssStart:ssStart+fps-1));
    SS_EMG(2,jj) = mean(emgAbsFilt(ssStart+fps:ssStart+fps*2-1));
    SS_ANG(1,jj) = mean(ang(ssStart:ssStart+fps-1));
    SS_ANG(2,jj) = mean(ang(ssStart+fps:ssStart+fps*2-1));

    if exist('Fdat','var')

        % Load tracked fascicle data
        FLx(:,1) = Fdat.Region(zz).FL;
        FAx(:,1) = rad2deg(Fdat.Region(zz).PEN);
        FVx(:,1) = Wvel(FLx(1:end-1),Time);

        % Upsample fascicle data rate to fps rate
        if length(FLx) > length(Time)
            FLx(end) = [];
            FAx(end) = [];
        elseif length(FLx) < length(Time)
            Time(end) = [];
        end

        FL = interp1(Time,FLx,time,'linear','extrap');
        FA = interp1(Time,FAx,time,'linear','extrap');

        if length(FVx) < length(Time)
            FV = interp1(Time(1:end-1),FVx,time,'linear','extrap');
        else
            FV = interp1(Time,FVx,time,'linear','extrap');
        end

        %% Calculate active fascicle force & length
        % Determine literature-based moment arms based on crank-arm
        % angle (0 deg = footplate perpendicular to tibia)
        MAang = polyval(MA,ang)/100; % cm to m

        % Determine active tendon force
        Ftendon = tqActive./MAang*0.5;

        % Determine active fascicle force
        Ffas = Ftendon./cosd(FA);

        %% Calculate active fascicle work
        % Determine fascicle displacement relative to a common start time
        tqStart = 1;
        if abs(FL(tqStart)-FL(tqStart+fps-1)) > 1
            error(['Person contracted at the start of trial ' files{jj} '.'])
        end

        %% Calculate active fascicle work relative to a common end time
        tqEnd = ssStart; % end of steady state
        if angEnd > ssStart && angDis > 10
            error(['Steady state definition started at ' num2str(round(ssStart/fps,1)) 's but rotation ended at ' num2str(round(angEnd/fps,1)) 's for trial ' files{jj} '.'])
        end

        FLdis = (FL-FL(tqStart))/1000;

        [WFAS(1,jj),WFAS(2,jj),WFAS(3,jj)] = work_calc(-FLdis(tqStart:tqEnd),Ffas(tqStart:tqEnd));

        %% Calculate TA MTU work from common start to end time
        % Determine MTU length
        MTU = 0.715+(-0.00130*gnegate(ang(tqStart:tqEnd))); % TA equation

        % Scale MTU length based on participant height
        MTU = MTU*((ht*0.25)/100);

        % Determine MTU velocity
        MTUV = Wvel(MTU,time);

        % Calculate MTU displacement
        MTUdis = MTU-MTU(1,1);

        % Use trapzoidal integration to calculate net MTU work
        [WMTU(1,jj),WMTU(2,jj),WMTU(3,jj)] = work_calc(-MTUdis(tqStart:tqEnd),Ftendon(tqStart:tqEnd));

        %% Calculate net ankle joint work from common start to end time
        % Calculate angular displacement
        angRad = deg2rad(ang);
        angDisRad = angRad-angRad(tqStart);

        % Use trapzoidal integration to calculate net joint work
        [WJNT(1,jj),WJNT(2,jj),WJNT(3,jj)] = work_calc(-angDisRad(tqStart:tqEnd),tqActive(tqStart:tqEnd));

        %% Check work calculation by plotting
        if workCheck

            h = figure(1);
            set(h,'Name',files{jj})

            % Joint work
            subplot(211); plot(tqActive); hold on; plot(tqStart,tqActive(tqStart),'ro'); plot(tqEnd+tqStart-1,tqActive(tqEnd+tqStart-1),'ro');

            % Angular rotation
            subplot(212); plot(ang); hold on; plot(angStart,ang(angStart),'bo'); plot(angEnd,ang(angEnd),'bo');

            waitfor(h);

        end

        %% Calculate fascicle stretch during rotation
        if abs(angDis) > 10
            % Use the actual fascicle length trace to calculate stretch
            [~,timeUSstart] = min(abs(Time-time(1)));
            FLx = FLx(timeUSstart:end);

            % Downsample crank-arm angle data to US frame rate
            angUS = interp1(time,ang,Time(timeUSstart:end));

            % Determine angular velocity
            angVelUS = Wvel(angUS,Time);

            % Find zero-crossing points
            idx = zci(angVelUS);

            % Determine maximum time difference in zero-crossing points
            [~,maxDiff] = max(diff(idx));

            % Determine start and end of rotation
            angStartUS = idx(maxDiff);
            angEndUS = idx(maxDiff+1);

            % Calculate local and net fascicle stretch
            fasStretch = fas_stretch(FLx(angStartUS:angEndUS));

            % Store fascicle stretch data
            FAS_STR(1,jj) = fasStretch(1,1); % local fascicle stretch
            FAS_STR(2,jj) = fasStretch(1,2); % cumulative fascicle stretch
            FAS_STR(3,jj) = fasStretch(1,3); % net fascicle stretch

            if FAS_STR(3,jj) >= 1
                error(['>1 mm of fascicle stretch during rotation so remove trial ' files{jj} '.'])
            end
        end

        %% Calculate fascicle shortening
        % Find the maximum shortening velocity
        FVy = min(FVx);
        FVX = find(FVx<=FVy*0.1,1,'first');

        % Determine the zero-crossing points until max. sho. vel.
        idx = zci(FVx(1:FVX));

        % The last zero-crossing point before max. sho. vel. should
        % roughly correspond to the instant of fascicle shortening
        FLstartShoX = idx(end-1);
        [~,FLstartSho] = max(FLx(FLstartShoX:FVX));
        FLstartSho = FLstartSho+FLstartShoX-1;

        % Store fascicle length at instant of fascicle shortening
        FLstart = FLx(FLstartSho);

        % Check fascicle shortening instant
        if abs(mean(FLx(1:17))-FLstart) > 0.5
            disp(['Check instant of fascicle shortening for trial ' files{jj} '.']);
        end

        % Store fascicle shortening data
        FAS_SHO(1,jj) = min(FLx)-FLstart; % max. fascicle shortening
        FAS_SHO(2,jj) = min(FLx); % min. fascicle length
        FAS_SHO(3,jj) = FLstart; % fascicle length before contraction

        %% Store fascicle data
        % Active fascicle force, length, and velocity
        FAS_FOR(:,jj) = Ffas(1:dur);
        FAS_LEN(:,jj) = FL(1:dur);
        FAS_VEL(:,jj) = FV(1:dur);

        %% Store data related to rotation
        if abs(angDis) > 10

            %ROTFASFOR(1,jj) = mean(Ffas(angStart:angEnd)); % mean fascicle force during the rotation
            ROTMTUFOR(1,jj) = mean(Ftendon(angStart:angEnd)); % mean MTU force during rotation
            %ROTFASVEL(1,jj) = mean(angVelUS(angStartUS:angEndUS)); % mean fasicle velocity during the rotation
            ROTMTUVEL(1,jj) = mean(MTUV(angStartUS:angEndUS)); % mean MTU velocity during the rotation
            ROTVEL(1,jj) = angVelMedian;
            ROTDIS(1,jj) = angDis;
            ROTMTUDIS(1,jj) = (max(MTU(1:ssStart))-min(MTU(1:ssStart)))*1000;

            %figure(2); hold on; plot(FL)

        else

            ROTVEL(1,jj) = angVelMedian;
            ROTDIS(1,jj) = angDis;

            %figure(2); hold on; plot(FL,':')

        end

    end

    %% Calculate TA MTU work from common start to end time
            % Determine MTU length
            % MTU = 0.715+(-0.00130*gnegate(ang)); % TA equation
            % 
            % % Scale MTU length based on participant height
            % MTU = MTU*((ht*0.25)/100);
            % 
            % % Determine MTU velocity
            % MTUV = Wvel(MTU,time);
            % 
            % % Calculate MTU displacement
            % MTUdis = MTU-MTU(1,1);
            % 
            % % Use trapzoidal integration to calculate net MTU work
            % [WMTU(1,jj),WMTU(2,jj),WMTU(3,jj)] = work_calc(-MTUdis(tqStart:tqEnd),Ftendon(tqStart:tqEnd));

            %% Calculate net ankle joint work from common start to end time
            % Calculate angular displacement
            angRad = deg2rad(ang);
            %angDisRad = angRad-angRad(tqStart);
            angDisRad = angRad;

            % Use trapzoidal integration to calculate net joint work
            %[WJNT(1,jj),WJNT(2,jj),WJNT(3,jj)] = work_calc(-angDisRad(tqStart:tqEnd),tqActive(tqStart:tqEnd));
            %[WJNT(1,jj),WJNT(2,jj),WJNT(3,jj)] = work_calc(-angDisRad,tqActive);
            %ROTMTUFOR(1,jj) = mean(Ftendon(angStart:angEnd)); % mean MTU force during rotation
            %ROTMTUVEL(1,jj) = mean(MTUV(angStartUS:angEndUS)); % mean MTU velocity during the rotation
            %ROTMTUDIS(1,jj) = (max(MTU(1:ssStart))-min(MTU(1:ssStart)))*1000;

            %% Check work calculation by plotting
            if workCheck && ~isnan(cond(5,jj))

                h = figure(1);
                set(h,'Name',files{jj})

                % Joint work
                subplot(211); plot(tqActive); 

                % Angular rotation
                subplot(212); plot(ang); hold on; plot(angStart,ang(angStart),'bo'); plot(angEnd,ang(angEnd),'bo');

                waitfor(h);

            end

    ROTEMG(1,jj) = mean(emgAbsFilt(angStart:angEnd)); % mean EMG during the rotation

    %% Plot reference trials
    if check && abs(angDis) < 10

        cnt(1) = cnt(1)+1;
        figure(3); subplot(311); hold on;
        sgtitle('Reference trials')
        plot(tqActive);
        %vline(ssStart,'k:'); vline(ssStart+fps,'k:'); vline(ssStart+fps*2-1,'k:')
        subplot(312); hold on;
        plot(emgAbsFilt);
        %vline(ssStart,'k:'); vline(ssStart+fps,'k:'); vline(ssStart+fps*2-1,'k:')
        subplot(313); hold on;
        plot(ang);
        %vline(ssStart,'k:'); vline(ssStart+fps,'k:'); vline(ssStart+fps*2-1,'k:')

    elseif check
        %% Plot dynamic trials

        cnt(2) = cnt(2)+1;
        figure(1); subplot(311); hold on;
        sgtitle('Dynamic trials')
        plot(tqActive);
        %vline(ssStart,'w:'); vline(ssStart+fps,'w:'); vline(ssStart+fps*2-1,'w:')
        subplot(312); hold on;
        plot(emgAbsFilt);
        %vline(ssStart,'w:'); vline(ssStart+fps,'w:'); vline(ssStart+fps*2-1,'w:')
        subplot(313); hold on;
        plot(ang);
        %vline(ssStart,'w:'); vline(ssStart+fps,'w:'); vline(ssStart+fps*2-1,'w:')
        figure(2); hold on; plot(angVel);
        sgtitle('Angular velocities')

    end

    clearvars -except thrPL pasFitFile mvcFile mvcMax emgCh fps ssStart dur ht zz check workCheck zci cnt files P MA mvc file cond TQ ANG EMG SS* W* FAS* ROT*
end

%% Calculate force depression relative to preload of 0% (rot_0)
% Condition (cond) variable stores the following:
% 1) Angular displacement in degrees
% 2) Mean angular velocity in degrees
% 3) Starting angle
% 4) Ending angle
% 5) Preload

valR = min(cond(3,:)); % maximum angle among all trials
colR = find(cond(3,:)<=valR+3); % allow for 3 deg of dynamometer variation

%% Remove trials if there steady-state torque varied by >= 1 Nm
% Loop through REF trials if >1
if length(colR) == 1
    % Remove REF trials if there was no steady state (torque within 10 Nm)
    if round(abs(SS_TQ(1,colR) -SS_TQ(2,colR)),1) >= 1
        error(['Remove reference trial ' file{colR} '.']);
    end
elseif length(colR) > 1 % See above for description
    for aa = 1:length(colR)
        if round(abs(SS_TQ(1,colR(aa)) -SS_TQ(2,colR(aa))),1) >= 1
            error(['Remove reference trial ' file{colR(aa)} '.']);
        end
    end
else
    error('No reference trial found.')
end

cnt(1,1) = sum(isnan(cond(5,:)));
cnt(1,2) = sum(~isnan(cond(5,:)));

% Find different conditions: low(L), mid(M), high(H)
PLmax = max(cond(5,:),[],'omitnan');
colH = find(cond(5,:)>=PLmax-thrPL);
PLmin = min(cond(5,:),[],'omitnan');
colL = find(cond(5,:)<=PLmin+thrPL);
dynFound = sort([colL colH]);
dyn = ~isnan(cond(5,:));
dyn = find(dyn==1);
dynLost = ismember(dyn,dynFound);
dynLost = find(dynLost==0);
dynLost = dyn(dynLost);
if length(dynLost)>1 && max(abs(diff(cond(5,dynLost)))) < 5
    colM = dynLost;
elseif length(dynLost) == 1
    colM = dynLost;
else
    error(['Remaining trials for preloads between ' num2str(PLmin+5) '-' num2str(PLmax-5) '%  are not matched with 10%.']);
end

% Average REF and DYN trials (column 1 = REF)
if exist('FAS_FOR','var')

    data = {TQ EMG ANG FAS_FOR FAS_LEN FAS_VEL}; % columns 1 to 6

else

    data = {TQ EMG ANG}; % columns 1 to 3

end

for aa = 1:size(data,2)
    data{2,aa}(:,1) = mean(data{1,aa}(:,colR),2); 
    data{2,aa}(:,2) = mean(data{1,aa}(:,colL),2);
    data{2,aa}(:,3) = mean(data{1,aa}(:,colM),2);
    data{2,aa}(:,4) = mean(data{1,aa}(:,colH),2);
end

%% Calculate rFD
for aa = 1:size(data,2)

    for bb = 1:4
        FD{aa}(1,bb) = mean(data{2,aa}(ssStart:ssStart+fps*2-1,bb),1);
    end

end

%% Calculate dFD
%P = plateau; A = ascending limb of force-length relation
% Perform time-matched comparisons
%refP = find(cond(3,:)>=13&cond(3,:)<=17);
%refA = find(cond(3,:)>=3&cond(3,:)<=7);

for aa = 1:size(data,2)
    %data{3,aa}(:,1) = mean(data{1,aa}(:,refP),2);
    data{3,aa}(:,1) = mean(data{1,aa}(:,colL),2);
    data{3,aa}(:,2) = mean(data{1,aa}(:,colM),2);
    data{3,aa}(:,3) = mean(data{1,aa}(:,colH),2);
    %data{4,aa}(:,1) = mean(data{1,aa}(:,refA),2);
    data{4,aa}(:,1) = mean(data{1,aa}(:,colL),2);
    data{4,aa}(:,2) = mean(data{1,aa}(:,colM),2);
    data{4,aa}(:,3) = mean(data{1,aa}(:,colH),2);
end

%% Determine time taken to reach REF angle (i.e., 17.5-12.5 or 7.5-2.5 deg)
t1 = round(median([find(data{2,3}(:,2)<=15,1,'first')...
    find(data{2,3}(:,3)<=15,1,'first')...
    find(data{2,3}(:,4)<=15,1,'first')]),-1);
t2 = round(median([find(data{2,3}(:,2)<=5,1,'first')...
    find(data{2,3}(:,3)<=5,1,'first')...
    find(data{2,3}(:,4)<=5,1,'first')]),-1);

%12500
%12250-12750 is 250 ms (2.5 deg)
%12400-12600 is 200 ms (+-1 deg)
%14650
if round(t1,-2) < t1
    t1 = round(t1,-2);
else
    t1 = round(t1,-2)-50;
end
if round(t2,-2) < t2
    t2 = round(t2,-2);
else
    t2 = round(t2,-2)-50;
end
s1 = t1-250+1;
s2 = t1+250;
s3 = t2-250+1;
s4 = t2+250;

for aa = 1:size(data,2)

    for bb = 1:3

        % 12250 is < 14400 & 15 deg is reached before 5 deg
        FD{3,aa}(1,bb) = mean(data{3,aa}(s1:s2,bb),1);
        FD{2,aa}(1,bb) = mean(data{4,aa}(s3:s4,bb),1);

    end
end

% FD in percent
for aa = 1:3

    for bb = 1:3

        for cc = 4:6

            if aa < 3 && cc-3 == 1 % Torque

                FD{cc,aa}(1,bb) = round((FD{cc-3,aa}(1,bb+1)-FD{cc-3,aa}(1,1))/FD{cc-3,aa}(1,1)*100,2);

            elseif aa < 3 && cc-3 > 1

                FD{cc,aa}(1,bb) = round((FD{cc-3,aa}(1,bb)-FD{cc-3,aa}(1,1))/FD{cc-3,aa}(1,1)*100,2);

            elseif aa == 4 && size(FD,2) > 3 % EMG

                FD{cc,aa}(1,bb) = round((FD{cc-3,aa}(1,bb)-FD{cc-3,aa}(1,1))/FD{cc-3,aa}(1,1)*100,2);

            elseif aa == 3 % Angle

                FD{cc,aa}(1,bb) = round(FD{cc-3,aa}(1,bb)-FD{cc-3,aa}(1,1),1);

            end

        end
    end
end

%% Check angle matching
%thrANG = 6;
thrEMG = 20.5;

% if max(abs(FD{4,3})) > thrANG
%     error(['>' num2str(thrANG) ' deg difference in angles: ' num2str(FD{1,3})])
% elseif max(abs(FD{5,3})) > thrANG
%     error(['>' num2str(thrANG) ' deg difference in angles: ' num2str(FD{2,3})])
% elseif max(abs(FD{6,3})) > thrANG
%     error(['>' num2str(thrANG) ' deg difference in angles: ' num2str(FD{3,3})])
% end

dynOrdered = [colL colM colH];

%% Check EMG matching
for cc = 1:length(dynOrdered)

    tqPreload(cc,1) = cond(5,dynOrdered(cc));
    emgMatch(cc,1) = round((mean(data{1,2}(ssStart:ssStart+fps*2-1,dynOrdered(cc)),1)-FD{1,2}(1,1))/FD{1,2}(1,1)*100,2);
    emgMatch(cc,3) = round((mean(data{1,2}(s1:s2,dynOrdered(cc)),1)-FD{3,2}(1,1))/FD{3,2}(1,1)*100,2);
    emgMatch(cc,2) = round((mean(data{1,2}(s3:s4,dynOrdered(cc)),1)-FD{2,2}(1,1))/FD{2,2}(1,1)*100,2);
    if max(abs(emgMatch(cc,:))) > thrEMG
        emgRemove(cc,1) = 1;
    else
        emgRemove(cc,1) = 0;
    end
    FDcheck(cc,1) = round((mean(data{1,1}(ssStart:ssStart+fps*2-1,dynOrdered(cc)),1)-FD{1,1}(1,1))/FD{1,1}(1,1)*100,2);
    FDcheck(cc,3) = round((mean(data{1,1}(s1:s2,dynOrdered(cc)),1)-FD{3,1}(1,1))/FD{3,1}(1,1)*100,2);
    FDcheck(cc,2) = round((mean(data{1,1}(s3:s4,dynOrdered(cc)),1)-FD{2,1}(1,1))/FD{2,1}(1,1)*100,2);
    emgFile{1,cc} = file{1,dynOrdered(cc)};

end

emgCheck = table(emgMatch(:,1),emgMatch(:,2),emgMatch(:,3),'VariableNames',{'-5','5','15'});
FDcheck = table(FDcheck(:,1),FDcheck(:,2),FDcheck(:,3),'VariableNames',{'-5','5','15'});

tEMG = table(tqPreload,emgRemove,emgCheck,FDcheck);
tEMG.Properties.RowNames = emgFile;
tEMG

%% Create FD summary table
cc = 0;
for aa = 4:6
    cc = cc+1;
    for bb = 1:3
        TQpcentFD(cc,bb) = FD{aa,1}(1,bb);
        EMGpcentFD(cc,bb) = FD{aa,2}(1,bb);
    end
end

TQpcentFD = table(TQpcentFD(:,1),TQpcentFD(:,2),TQpcentFD(:,3),'VariableNames',{'PL0','PL25','PL50'});
EMGpcentFD = table(EMGpcentFD(:,1),EMGpcentFD(:,2),EMGpcentFD(:,3),'VariableNames',{'PL0','PL25','PL50'});
tFD = table(TQpcentFD,EMGpcentFD);
tFDRows = {'-5 deg','5 deg','15 deg'};
tFD.Properties.RowNames = tFDRows;
tFD

%% Remove trials that were not EMG matched within 20%
if sum(emgRemove) > 0
    ques = input(['Rerun code after removing trials that were not EMG matched within ' num2str(thrEMG) '%? y/n: '],'s');
    if strcmpi(ques,'n')
        disp(['Trials that were not EMG matched within ' num2str(thrEMG) '% are included in the TQanalysed output.'])
        %save('TQanalysed.mat')
    elseif isempty(ques) || strcmpi(ques,'y')
        files(dynOrdered(emgRemove==1)) = [];
        getFD
    end
else
    %save('TQanalysed.mat')
    %error('Something went wrong with saving the file.')
end

%% Plot time traces for rFD-related trials
figure(4); hold on
sgtitle('rFD comparisons')

for cc = 1:3

    subplot(3,1,cc); hold on; plot(data{2,cc});
    vline(ssStart,'w:'); vline(ssStart+fps*2-1,'w:')

end

%% Plot force-velocity points for dFD-related trials
figure(5);
sgtitle('dFD comparisons')
tcl = tiledlayout(1,3);
conds = {'-5 deg';'5 deg';'15 deg'};

for cc = 1:3

    nexttile(tcl);

    for dd = 1:3

        if size(FD,2) > 3

            plot(FD{cc,6}(dd)*-1,FD{cc+3,1}(dd),'o'); hold on

        else

            plot(0,FD{cc+3,1}(dd),'o'); hold on

        end

    end

    axis square
    title(conds{cc})

end

hL = legend({'PL0','PL25','PL50'},'location','eastoutside');
legend boxoff

% Variables you should check
%mean(emgCheck,2)
%cond(5,:)

% Export relevant data to excel spreadsheet
isNested = varfun(@istable, tFD, 'OutputFormat','uniform');   % true for table variables
nestedNames = tFD.Properties.VariableNames(isNested);
t2 = splitvars(tFD);
%writetable(t2,'FDresults.xlsx')

% Update vertical lines in figures
figure(1);
subplot(311); vline(ssStart,'w:'); vline(ssStart+fps,'w:'); vline(ssStart+fps*2-1,'w:')
subplot(312); vline(ssStart,'w:'); vline(ssStart+fps,'w:'); vline(ssStart+fps*2-1,'w:')
subplot(313); vline(ssStart,'w:'); vline(ssStart+fps,'w:'); vline(ssStart+fps*2-1,'w:')

figure(2);

figure(3);
subplot(311); vline(ssStart,'w:'); vline(ssStart+fps,'w:'); vline(ssStart+fps*2-1,'w:')
subplot(312); vline(ssStart,'w:'); vline(ssStart+fps,'w:'); vline(ssStart+fps*2-1,'w:')
subplot(313); vline(ssStart,'w:'); vline(ssStart+fps,'w:'); vline(ssStart+fps*2-1,'w:')

figure(4);

figure(5);