%-------------------------------------------
% Nonlinear cubic oscillator 
% Junting Wang 10/02/2024
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

% Check that stability condition is satisfied

if(SR<=pi*f0)
    error ('Stability condition violated')
end

% derived quantities

k = 1/SR;               % time step (s)
w0 = 2*pi*f0;           % angular frequency (rad./s)
Nf = floor(Tf*SR);      % total number of time steps

% initialize

u2 = u0;                % set initial displacement u2 : u ^{n-2}
u1 = u0+k*v0;           % set second displacement  u1  : u^{n-1} 

out = zeros(Nf,1);      % output vector

% main loop

tic
for n=1:Nf
    b = 2/(1+(w0^4*k^2*(u1^2))/2);
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
title('Nonlinear cubic oscillator');

% play sound

soundsc(out,SR)
