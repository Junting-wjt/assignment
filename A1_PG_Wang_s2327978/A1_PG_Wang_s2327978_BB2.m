%-------------------------------------------------
% PBMMI_Assignment01 - KarplusStrong
% 
% Create chords
% 
% Junting Wang 02/02/24
%-------------------------------------------------


% Clear the command window, workspace and close all plots -----------------
clc;                                                 % clear the command window
clear;                                               % clear workspace
close all;                                           % close all plots


% Set the governing parameters for script ---------------------------------
Fs = 44.1e3;                                         % the sampling rate in Hz
T60 = 2;                                             % duration of each note in seconds
overlap_dur = 1.92;                                  % duration of overlap between each note in seconds
f0 = 880;                                            % the (desired) fundamental frequency of the string in Hz
f1 = f0*5/4;
f2 = f0*3/2;
f0_chord = [f0, f1, f2];                             % fundamental frequencies for a major chord
rho = exp(-6.91/(T60*f0))/cos(pi*f0/Fs);             % the loss parameter, ρ
R = 0.95;                                            % the dynamics parameter

% Calculate the chords duration with overlaps
chords_dur = T60 + (length(f0_chord) - 1) * (T60 - overlap_dur);
y_chords = zeros(1, round(chords_dur*Fs));           % the output of the chords


% Implement algorithm for each note ---------------------------------------
start_index = 1;                                     % start for the first note
for f0 = f0_chord
   % Calculate the derived parameters
   M = round(T60 * Fs);                                 % duration of simulation in samples
   Nexact = Fs/f0 - 0.5;                                % the ideal number of samples in the delay line
   N = floor(Nexact);                                   % the integer part of the delay
   P = Nexact - N;                                      % the fractional delay
   C = (1-P)/(1+P);                                     % allpass filter coefficient

   % Initialise 2 vectors
   v = 2*rand(1,N+1) - 1;
   y = zeros(1,M);                                      % the output of the KS algorithm
   yp0 = zeros(1,M);                                    % the output of the allpass filter

   % Implement the dynamics filter 
   x1 = 0;                                              % initialise a (scalar) state variable                                                  
   for n = 0:N
       x0 = (1-R)*v(n+1) + R*x1;                        % read from white noise vector
       y(n+1) = x0;                                     % write x0 into the output vector
       x1 = x0;                                         % update the state variable by copying x0 into x1
   end
 
   % Main Karplus-Strong algorithm 
   yp1 = 0;
   for n = N+1:M-1
       yp0 = C*y(n-N+1) + y(n-N) - C*yp1;
       y(n+1) = (yp0 + yp1)/2*rho;
       yp1 = yp0;
   end

   % Normalize the output vector
   y = y / max(abs(y));

   % Calculate the end index 
   end_index = start_index + M - 1;
   end_index = min(end_index, length(y_chords));

   % Mix each currente note
   y_chords(start_index:end_index) = y_chords(start_index:end_index) + y(1:length(start_index:end_index));

   % Update start index for the next note
   start_index = start_index + round((T60 - overlap_dur) * Fs);
end


% Listen to the results ---------------------------------------------------
soundsc(y_chords, Fs);

