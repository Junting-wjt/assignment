%-------------------------------------------
% Simple harmonic oscillator (basic FD scheme)
% Junting Wang 06/02/2024
%-------------------------------------------

clc
clear
close all

% input parameters

SR = 44100;     % sample rate (Hz)
f0 = 1e3;       % frequency (Hz)
Tf = 1;         % duration (s)
u0 = 1;         % initial displacement
v0 = 0;         % initial velocity
alpha = 0.8;    % scheme free parameter

% derived quantities

k = 1/SR;               % time step (s)
w0 = 2*pi*f0;           % angular frequency (rad./s)
Nf = floor(Tf*SR);      % total number of time steps
b = (2-alpha*w0^2*k^2)/(1+((1-alpha)*w0^2*k^2)/2);  % scheme parameter

% initialize

u2 = u0;                % set initial displacement u2 : u ^{n-2}
u1 = u0+k*v0;           % set second displacement  u1  : u^{n-1} 

out = zeros(Nf,1);      % output vector

% Verify stability condition

if alpha >= 1/2 && k >= 2/(w0*(sqrt(2*alpha - 1)))
    error('stability condition is violated');
end

% main loop

tic
for n=1:Nf
    u = b*u1-u2;        % scheme update u^{n} u^{n} = b * u^{n-1} - u^{n-2}
    out(n) = u2;        % write output
    
    u2 = u1;            % shift state
    u1 = u;
end
toc

% plot

tax = [0:Nf-1]'*k;
plot(tax,out);
xlabel('t');
ylabel('u');
title('Simple Harmonic Oscillator');

% play sound

soundsc(out,SR)
