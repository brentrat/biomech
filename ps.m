function ps(signal,fps)
%*Compute power spectrum in dB over signal frequency*
%ps(signal,fps)
%Signal demo if no inputs are provided.

%   Code updated by Brent J. Raiteri on 21.10.2026 

if nargin == 0
    signal = zeros(1000,1);
    signal([1 51:100 151:200 251:300 351:400 451:500 551:600 651:700 751:800 851:900 951:1000]) = 1;
    signal = signal-mean(signal);
    fps = length(signal);
    figure; plot(signal)
    % Because this is a harsh square-like wave, the power spectrum graph 
    % will show a  dominant spike at 10 Hz (the fundamental frequency), 
    % followed by smaller spikes at its odd harmonic frequencies 
    % (30 Hz, 50 Hz, 70 Hz, etc.).
end
%signal = signal(1:fps);

% Transpose to column vector if it is a row vector
if isrow(signal)
    signal = signal';
end

% Remove DC offset so 0 Hz doesn't dominate the peak search
% That is, fPeak (see below) can be at 0 Hz due to a DC offset, which is
% a non-zero average, representing constant unmoving energy in the signal.
signal = signal - mean(signal); 

% Define other aspects of the signal
N = length(signal);
dt = 1/fps;
time = (dt:dt:dt*N)';

%% Convert to timetable for modern MATLAB pspectrum usage
% see pspectrum documentation
xTable = timetable(seconds(time),signal);

%% Open figure and plot the spectrum
%[pxx,f] = pspectrum(xTable);
%f = f/1000;
%plot(f,pow2db(pxx))
%xlabel('Frequency (Hz)')
%ylabel('Power Spectrum (dB)')
figure;
pspectrum(xTable,'FrequencyResolution', 1); % Forced 1 Hz resolution (1-s signal duration required!)
hold on;

% Style axes labels
fig = gca;
tit = fig.Title.String;
%title(tit,'fontweight','normal');
xTit=get(get(fig,'xlabel'),'string');
yTit=get(get(fig,'ylabel'),'string');
xlabel(xTit,'fontweight','bold');
ylabel(yTit,'fontweight','bold');

% Extract frequency (XData) and power (YData) from the plotted line
h = findobj(gca,'Type','line');
%h(1).LineStyle = ':'; % uncomment to change from solid line to dots
x = h(1).XData; % frequency vector
y = h(1).YData; % power vector (pspectrum plots in dB by default)

% Convert dB back to linear power for cumulative energy math
% Cumulative sums (cumsum) must be computed on linear power (\(W\) or 
% \(V^{2}\)), not logarithmic power (dB). 
% Summing dB values yields mathematically invalid results. 
yLinear = db2pow(y);
yCum = cumsum(yLinear);
yCumNorm = yCum / yCum(end);

% Find index where 99% of total cumulative power is reached
% idx99 is mapped  to the frequency axis instead of using the raw amplitude
idx99 = find(yCumNorm >= 0.99, 1, 'first');
f99 = x(idx99); % actual frequency value

% Find index of the peak dominant frequency (highest power)
[~, idxPeak] = max(yLinear);
fPeak = x(idxPeak);

% --- Visual markers --- %
% Get unit scaling factor from x-axis label
xAxisLabel = get(get(fig, 'xlabel'), 'string');
if contains(xAxisLabel, '(kHz)', 'IgnoreCase', true)
    unitScale = 1000; % Convert kHz to Hz
elseif contains(xAxisLabel, '(MHz)', 'IgnoreCase', true)
    unitScale = 1000000; % Convert MHz to Hz
else
    unitScale = 1; % Already in Hz
end

% Convert calculated frequencies strictly to Hz for the text labels
f99_Hz = f99 * unitScale;
fPeak_Hz = fPeak * unitScale;

% Vertical dashed line for 95% power bandwidth
xline(f99, 'b--', {['99% Power Cutoff: ', num2str(round(f99_Hz, 1)), ' Hz']}, ...
    'LabelVerticalAlignment', 'bottom', 'LabelHorizontalAlignment', 'right');

% Vertical dotted line for peak frequency
xline(fPeak, 'k:', {['Peak Freq: ', num2str(round(fPeak_Hz, 1)), ' Hz']}, ...
    'LabelVerticalAlignment', 'top', 'LabelHorizontalAlignment', 'left');

hold off;
end
