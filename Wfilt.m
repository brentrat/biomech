function yF = Wfilt(y,fc,type,fs)
%*Filter signal with a 2nd-order Butterworth filter*
%Wfilt(signal, cut-off_frequency_in_Hz, filter_type, sampling_frequency)
%*Requires Wcu function.

% Author:
% BJ Raiteri, April 2021, if you find errors pls email brent.raiteri@rub.de
% tested in R2019b

% Change signal to column vector
if isrow(y)
    y = y';
end

% Remove NaNs from signal
if sum(isnan(y)) > 0 
    tf = isnan(y);
    idx = 1:numel(y);
    y(tf) = interp1(idx(~tf),y(~tf),idx(tf));
    if ~isempty(isnan(y))
        cut = find(isnan(y),1,'first');
        y = y(1:cut-1);
    end      
end

% Change filter type based on type
if strcmp(type,'low')
    [cu,order] = Wcu(4, fc, 'low', fs);   % 1st input = 4 for 2nd-order filter
    [B,A] = butter(order,cu/(fs/2),'low');
elseif strcmp(type,'high')
    [cu,order] = Wcu(4, fc, 'high', fs);   % 1st input = 4 for 2nd-order filter
    [B,A] = butter(order,cu/(fs/2),'high');
elseif strcmp(type,'bandpass')
    if size(fc,2) ~= 2
        error('Provide two cut-off frequencies for a bandpass filter.')
    end
    cu = fc;
    order = 2;  
    [B,A] = butter(order,cu/(fs/2),'bandpass');
else
    error('Enter a valid filter type [low high bandpass stop].')
end

% Filter signal 
yF = filtfilt(B,A,y);