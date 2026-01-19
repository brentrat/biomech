function pSpec(signal,fps)
%*Calculate single-sided power spectrum of signal*
%pSpec(signal,fps)

%   Code updated by Brent J. Raiteri on 19.01.2026 
%   from Electromyogram analysis by William Rose
%   KAAP686 Mathematics and Signal Processing for Biomechanics

% Create column vector of frequencies up to Nyquist frequency
N = length(signal);
freqs(:,1) =0:fps/N:(fps/2);

% Compute fft
xfft = fft(signal-mean(signal));

% Plot fft up to Nyquist frequency
%plot(freqs,abs(xfft(1:N/2+1)));

% Compute and plot the power spectrum up to Nyquist frequency
%Pxx = xfft.*conj(xfft);

% Ensure peaks represent average power contributed by each frequency
% and that power does not increase for longer signals 
Pxx = (abs(xfft)/N).^2;
plot(freqs,abs(Pxx(1:N/2+1)));
xlabel('Frequency (Hz)','fontweight','bold')
ylabel('Power (Magnitude^2)','fontweight','bold')

end