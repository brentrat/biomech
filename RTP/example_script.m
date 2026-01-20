clear; close all

% Analyse dynamometry data from Tim's project
% Load the dynamometry data
pname = 'C:\Users\Brent\Desktop\Research\Tim\dynamometry\p1\MAT data';
cd(pname)
folders = dir;

for ff = 1:length(folders)
    if contains(folders(ff).name(1),'p') ||  contains(folders(ff).name(1),'1')
        keep(ff,1) = 1;
    else
        keep(ff,1) = 0;
    end
end

folders(keep==0) = [];

for ff = 1:length(folders)

    cd(pname)
    cd(folders(ff).name)
    cd("MAT data")

    files = dir;

    for trial = 1:length(files)
        if size(files(trial).name,2) >= 4
            if contains(files(trial).name(1:4),'sine')
                keepFile(trial,1) = 1;
            else
                keepFile(trial,1) = 0;
            end
        end
    end

    files(keepFile==0) = [];

    for trial = 1:length(files)

        % Load data
        load(files(trial).name)
        
        % Get sampling frequency + plot power spectrum
        fs = 1/Torque.interval;
        pSpec(Torque.values,fs)

        % Filter torque and plot the raw and filtered signals
        tqRaw = Torque.values;
        tq = Wfilt(tqRaw,10,'low',fs);
        figure(1); plot(Torque.times,tqRaw);
        hold on; plot(Torque.times,tq)

        % Sync data
        time(:,1) = Torque.times;
        first = find(time>=US.times(1),1,'first');
        last = find(time<=US.times(end),1,'last');
        time = time(first:last);
        time = time-US.times(1);
        tq = tq(first:last);

        % Linear interpolation
        % x = current time points
        % v = current amplitudes
        % xq = desired time points
        ang = interp1(Angle.times,Angle.values,Torque.times,'linear','extrap');

        % I could also just use filtTAE
        [tq,ang,emg,us,mot,time,fname] = filtTAE(fs,files(trial).name);

        % RMS amplitude for raw EMG channels
        for emgCh = 1:size(emg,2)
            emgRMS(:,emgCh) = rmsDC(EMG(:,emgCh),250,1,fs);
        end

    end

end