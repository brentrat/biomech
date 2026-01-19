function [tq,ang,emg,us,mot,time,fname] = filtTAE(fs,fname)
%*Filter signals with a dual-pass 2nd-order Butterworth filter*
%filtTAE(sampling_frequency,file_name)
%*Requires Wfilt & Wcu functions.

%% Inputs
tqCU = 20;
angCU = 6;

%% Load spike data exported as .mat
if nargin < 2
    [fname,pname] = uigetfile('*.mat','Choose the Spike file to analyse.');
    cd(pname);
end
data = load(fname);

%% Get variable names
vars = fieldnames(data);
fns = {'Angle';'Torque'};
match = contains(vars,fns);
vars(match==1) = [];

%% Check fps
if sum(match) == 2
    if 1/fs ~= data.Torque.interval && 1/fs ~= data.Angle.interval
        error('Wrong fps input.')
    end
else
    error('Error processing Spike2 file. Angle or torque variable does not exist.')
end

%% Sync times
time(:,1) = 0:1/fs:1/fs*(length(data.Torque.values)-1);
if isfield(data,'Keyboard') && ~isempty(data.Keyboard.times)
    [~,start] = min(abs(time-data.Keyboard.times(1)));
else
    start = 1;
end

%% Filter, interpolate, and crop data
% Torque
tq(:,1) = Wfilt(data.Torque.values,tqCU,'low',fs);
tq = interp1(data.Torque.times,tq,time);
tq = tq(start:end);

% Angle
ang(:,1) = Wfilt(data.Angle.values,angCU,'low',fs);
ang = interp1(data.Angle.times,ang,time);
ang = ang(start:end);

% Other waveforms e.g., EMG
matchEMG = contains(vars,'EMG');
varsEMG = vars;
varsEMG(matchEMG==0) = [];

if ~isempty(varsEMG)
    [b,a] = butter(4,[20 400]/(fs/2));
    for jj = 1:length(varsEMG)
        emg(:,jj) = filtfilt(b,a,data.(varsEMG{jj}).values-mean(data.(varsEMG{jj}).values)); % remove DC offset
        emgT(:,jj) = data.(varsEMG{jj}).times;
        emg(:,jj) = interp1(emgT(:,jj),emg(:,jj),time);
    end

    emg = emg(start:end,:);
    %     emgAbs = abs(emg);
    %     emgAbsFilt = Wfilt(emgAbs,10,'low',fps);

else
    matchEMG = contains(vars,'TA');
    varsEMG = vars;
    varsEMG(matchEMG==0) = [];

    if ~isempty(varsEMG)
        [b,a] = butter(4,[20 400]/(fs/2));
        for jj = 1:length(varsEMG)
            emg(:,jj) = filtfilt(b,a,data.(varsEMG{jj}).values-mean(data.(varsEMG{jj}).values)); %remove DC offset
            emgT(:,jj) = data.(varsEMG{jj}).times;
            emg(:,jj) = interp1(emgT(:,jj),emg(:,jj),time);
        end

        emg = emg(start:end,:);
    else
        emg = [];
    end
end

emg = array2table(emg);
emg.Properties.VariableNames = varsEMG;

% Events
if isfield(data,'US') && isfield(data,'Keyboard')
    us(:,1) = data.US.times-data.Keyboard.times(1);
elseif isfield(data,'US')
    us(:,1) = data.US.times;
else
    us = [];
end

if isfield(data,'Motive') && isfield(data,'Keyboard')
    mot(:,1) = data.Motive.times-data.Keyboard.times(1);
elseif isfield(data,'Motive')
    mot(:,1) = data.Motive.times;
else
    mot = [];
end

% Crop time vector
if isfield(data,'Keyboard') && ~isempty(data.Keyboard.times)
    time = time(start:end)-data.Keyboard.times(1);
end
