function pk = fft_plot(signal,fps)
%*Calculate single-sided amplitude spectrum of signal*
%fft_plot(signal,fps)

%https://au.mathworks.com/help/matlab/ref/fft.html

y = fft(signal);
n = length(signal);          % number of samples
P2 = abs(y/n);          % compute two-sided spectrum P2
P1 = P2(1:n/2+1);       % compute single-sided spectrum
P1(2:end-1) = 2*P1(2:end-1);
f = fps*(0:(n/2))/n;
%fig = figure(500);
plot(f,P1)
[~,pkX] = max(P1);
pk = f(pkX);
xlabel('Frequency (Hz)','fontweight','bold')
ylabel('Amplitude','fontweight','bold')
%ax = findall(fig, 'type', 'axes');
%ax = ax(strcmp('', get(ax, 'Tag')));
%set(ax, 'box', 'off','FontName' ,'Arial','TickDir','out','TickLength',[.0075 .0075])
%set(ax,'linewidth',1)