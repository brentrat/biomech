function xRect = rectify(x,fsamp)
%*Rectify and low-pass EMG signal at 10 Hz
%rectify(x, fsamp)

% Author:
% BJ Raiteri, 08/2024, if you find errors pls email brent.raiteri@rub.de
% tested in R2022a
xAbs = abs(x-mean(x));
xRect = Wfilt(xAbs,10,'low',fsamp);

end