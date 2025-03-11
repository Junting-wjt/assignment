%-------------------------------------------------
% PBMMI_Assignment01 - KarplusStrong
% 
% Coding the basic Karplus-Strong algorithm
% 
% Junting Wang 30/01/24
%-------------------------------------------------


% Clear the command window, workspace and close all plots -----------------
clc;                                                 % clear the command window
clear;                                               % clear workspace
close all;                                           % close all plots


% Set the governing parameters for script ---------------------------------
Fs = 44.1e3;                                         % the sampling rate in Hz
dur = 2;                                             % duration of simulation in seconds
f0 = 880;                                            % the (desired) fundamental frequency of the string in Hz
rho = 0.998;                                         % the loss parameter, ρ
R = 0.95;                                            % the dynamics parameter


% Calculate the derived parameters ----------------------------------------
M = round(dur * Fs);                                 % duration of simulation in samples
N = floor(Fs/f0 - 1/2);                              % the (truncated) delay line length


% Initialise 2 vectors ----------------------------------------------------
%rng(0)
v = 2*rand(1,N+1) - 1;                               % a vector of white noise 
y = zeros(1,M);                                      % the ouput of the KS algorithm


% Implement the dynamics filter -------------------------------------------
x1 = 0;                                              % initialise a (scalar) state variable                                                  
for n = 0:N
    x0 = (1-R)*v(n+1) + R*x1;                        % read from white noise vector
    y(n+1) = x0;                                     % write x0 into the output vector
    x1 = x0;                                         % update the state variable by copying x0 into x1
end

 
% Main Karplus-Strong algorithm -------------------------------------------
for n = N+1:M-1
    y(n+1) = (y(n - N + 1) + y(n - N))/2*rho;
end


% Listen to the results ---------------------------------------------------
soundsc(y, Fs);


% Plot the output waveform and the spectrum -------------------------------
t = (0:M-1)/Fs;                                       % time vector
f = (0:round(M/2))*Fs/M;                              % frequency vector
Y = fft(y);

subplot(2,1,1);
plot(t, y);
title('Output Waveform');
xlabel('Time (s)');
ylabel('Amplitude');

subplot(2,1,2);
plot(f, abs(Y(1:round(M/2)+1)));
title('Spectrum');
xlabel('Frequency (Hz)');
ylabel('Magnitude');
xlim([0 Fs/2]);
vertical_line = xline(f0, '--');
legend(vertical_line, 'f0');

