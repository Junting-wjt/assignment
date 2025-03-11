%++++++++++++++++++++++++++++++++++++++++
% Moog VCF BB2
%
% Wang s2327978
% 24 February 2024
%++++++++++++++++++++++++++++++++++++++++

clc
clear
close all

%++++++++++++++++++++++++++++++++++++++++
% input parameters
SR = 44100 ;               % sample rate [Hz]
Tf = 20;                   % total simulation time [s]
f0 = 120 ;                 % resonant filter frequency [Hz]
r = 0.7 ;                  % feedback coeff [choose 0 \leq r \leq 1]

% derived parameters
om0 = 2*pi*f0;             % angular resonant frequency
Nf = floor(Tf*SR);         % total number of samples for the simulation
k = 1/SR;                  % time step
a = 2^0.5*r^0.25;          % derived parameter to check stability condition 

%++++++++++++++++++++++++++++++++++++++++
% check if the stability condition for Forward Euler is satisfied
if k >= (2^0.5*a+2) / (om0*(a^2+2^0.5*a+1))...
   || k >= (2-2^0.5*a) / (om0*(a^2-2^0.5*a+1))
    error('Stability condition is violated');
end

% initialise the matrices and vectors 
I = eye(4);                % 4 by 4 identity matrix
xf = zeros(4, 1);          % initialize 4 by 1 vector for state
yf = zeros(Nf, 1);         % initialise Nf by 1 vector to hold output y
[u, SR] = audioread("PourIrina.wav");
tvec = (0:Nf-1)'*k;        % time vector for plots
fvec = (0:Nf-1)'*SR/Nf;    % frequency vector for plots

%++++++++++++++++++++++++++++++++++++++++
% main loop
for n = 1 : Nf
    % Update state xf with the nonlinear state space function
    x_state = nonlinear_state_space(xf, u(n), om0, r);
    xf = xf + k*x_state;   
    % Write sample to the output vector
    yf(n) = xf(4);
end
%++++++++++++++++++++++++++++++++++++++++

% shift the value of the obtained yf vector to the right as a whole
yf = circshift(yf,1);

% discrete transfer functions
Hf = fft(yf);

% Plot the output signal
subplot(2,1,1)
plot(tvec, yf);
title('Output Signal');
xlabel('Time (s)');
ylabel('Amplitude');

% plot the magnitudes of the transfer functions 
subplot(2,1,2)
loglog(fvec, abs(Hf))
title('Transfer Function Magnitudes');
xlabel('Frequency (Hz)');
ylabel('Magnitude');

% play the output signal
soundsc(yf,SR);

% Nonlinear state space function
function x_state = nonlinear_state_space(x, u, om0, r)
    x_state = zeros(4,1);
    x_state(1) = om0*(-tanh(x(1)) - tanh(4*r*x(4) + u));
    x_state(2) = om0*(-tanh(x(2)) + tanh(x(1)));
    x_state(3) = om0*(-tanh(x(3)) + tanh(x(2)));
    x_state(4) = om0*(-tanh(x(4)) + tanh(x(3)));
end
