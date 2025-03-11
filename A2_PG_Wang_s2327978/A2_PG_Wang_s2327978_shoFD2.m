%-------------------------------------------
% Simple harmonic oscillator (exact FD scheme)
% Junting Wang 06/02/2024
%-------------------------------------------

clc
clear 
close all

% input parameters

SR = 44100;     % sample rate (Hz)
f0 = 1e3;       % frequency (Hz)
Tf = 1;         % duration (s)
u0 = 0.3;       % initial displacement
v0 = 0;         % initial velocity

% Check that stability condition is satisfied

if(SR<=pi*f0)
    error ('Stability condition violated')
end

% derived quantities

k = 1/SR;               % time step (s)
w0 = 2*pi*f0;           % angular frequency (rad./s)
Nf = floor(Tf*SR);      % total number of time steps
b = 2*cos(w0*k);        % scheme parameter

% initialize

u2 = u0;                % set initial displacement u2 : u ^{n-2}
% initial conditions such that an exact output, set second displacement  u1  : u^{n-1} 
u1 = u0 * cos(w0 * k) + (v0 / w0) * sin(w0 * k); 

out = zeros(Nf,1);      % output vector

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

subplot(2,1,1);
tax = [0:Nf-1]'*k;

% calculate analytic solutions at time tax
out_an = u0 * cos(w0 * tax) + (v0 / w0) * sin(w0 * tax);

% calculate error
error = out - out_an;

plot(tax,out);
xlabel('t');
ylabel('u');
title('Simple Harmonic Oscillator');

subplot(2,1,2);
plot(tax,error);
xlabel('Time (s)');
ylabel('Error');
title('Error between FDTD and Analytical Solutions');

% play sound

soundsc(out,SR)
