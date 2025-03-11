%-------------------------------------------------
% PBMMI_Assignment01 - KarplusStrong
% 
% Glissandi and rho related to T60
% 
% Junting Wang 01/02/24
%-------------------------------------------------


% Clear the command window, workspace and close all plots -----------------
clc;                                                 % clear the command window
clear;                                               % clear workspace
close all;                                           % close all plots


% Set the governing parameters for script ---------------------------------
Fs = 44.1e3;                                         % the sampling rate in Hz
T60 = 2;                                             % duration of simulation in seconds
f0 = 880;                                            % the (desired) fundamental frequency of the string in Hz
R = 0.95;                                            % the dynamics parameter


% Calculate the derived parameters ----------------------------------------
M = round(T60 * Fs);                                 % duration of simulation in samples
Nexact = Fs/f0 - 0.5;                                % the ideal number of samples in the delay line
N= floor(Nexact);                                    % the integer part of the delay


% Set and calculate parameters for glissandi ------------------------------
f1 = 440;                                            % the target frequency
Ntarget = Fs/f1 - 0.5;                               % the ideal number of samples in target delayline

% Parameters for controlling glissandi speed
glissandi_dur = 0.5;                                 % duration of glissandi
glissandi_sample = glissandi_dur * Fs;               % samples of glissandi
% delay line changed for each sample
delayline_change = (Ntarget - Nexact) / glissandi_sample; 


% Initialise 2 vectors ----------------------------------------------------
v = 2*rand(1,N+1) - 1;
y = zeros(1,M);                                      % the output of the KS algorithm
yp0 = zeros(1,M);                                    % the output of the allpass filter


% Implement the dynamics filter -------------------------------------------
x1 = 0;                                              % initialise a (scalar) state variable                                                  
for n = 0:N
    x0 = (1-R)*v(n+1) + R*x1;                        % read from white noise vector
    y(n+1) = x0;                                     % write x0 into the output vector
    x1 = x0;                                         % update the state variable by copying x0 into x1
end


% Main Karplus-Strong algorithm -------------------------------------------
yp1 = 0;
for n = N + 1:M - 1
    if n <= Nexact + glissandi_sample
        % change delay line
        Nexact = Nexact + delayline_change;
        N = floor(Nexact);
        % the fractional delay of each delay line
        P = Nexact - N;
        C = (1-P) / (1+P); 
        % fundamental frequency of each delay line
        f = Fs/N;
        % calculate rho with f and user-selected T60
        rho = exp(-6.91/(T60*f))/cos(pi*f/Fs);
    end
    yp0 = C*y(n-N+1) + y(n-N) - C*yp1;
    y(n+1) = (yp0 + yp1)/2*rho;
    yp1 = yp0;
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
vertical_line_f0 = xline(f0, '--');
vertical_line_f1 = xline(f1,'--','Color','r');
legend([vertical_line_f0, vertical_line_f1], {'f0', 'f1'});

