%-------------------------------------------
% Energy for simple harmonic oscillator 
% Junting Wang 08/02/2024
%-------------------------------------------

clc
clear
close all

% input parameters

SR = 44100;         % sample rate (Hz)
f0 = 1e3;           % frequency (Hz)
Tf = 1;             % duration (s)
u0 = 1;             % initial displacement
v0 = 0;             % initial velocity
alpha = 0.8;        % scheme free parameter

% derived quantities

k = 1/SR;                % time step (s)
w0 = 2*pi*f0;            % angular frequency (rad./s)
Nf = floor(Tf*SR);       % total number of time steps
b = (2-alpha*w0^2*k^2)/(1+((1-alpha)*w0^2*k^2)/2);  % scheme parameter

% initialize

u2 = u0;                 % set initial displacement u2 : u ^{n-2}
u1 = u0+k*v0;            % set second displacement  u1  : u^{n-1} 
out = zeros(Nf,1);       % output vector
energy = zeros(Nf,1);    % energy vector
energy_u0 = 1/2*v0^2+1/2*w0^2*u0^2;  % initial energy

% Verify stability condition

if alpha >= 1/2 && k >= 2/(w0*(sqrt(2*alpha - 1)))
    error('stability condition is violated');
end

% main loop

for n=1:Nf
     u = b*u1-u2;   % scheme update u^{n} u^{n} = b * u^{n-1} - u^{n-2}
     out(n) = u2;   % write output

     if n==1
         energy(n) =energy_u0;
     else 
          if n~=Nf
          energy(n+1) = (1/(2*k^2)) * (u1^2 + (w0^2 * k^2 - 2)*u2*u1 + u2^2);
          end
     end

     u2 = u1;       % shift state
     u1 = u;
end

energy(2)= (1/(2*k^2)) * (out(2)^2 + (w0^2 * k^2 - 2)*out(2)*u0 + u0^2);

% plot

subplot(3,1,1);
tax = [0:Nf-1]'*k;
plot(tax,out);
xlabel('t');
ylabel('u');
title('Simple Harmonic Oscillator');

subplot(3,1,2);
plot(tax, energy);
xlabel('t');
ylabel('h');
title('Numerical Energy of the SHO');

% Calculate energy error

energy_error = (energy - energy_u0)/energy_u0;
subplot(3,1,3); 
plot(tax, energy_error);
xlabel('t');
ylabel('Energy Error');
title('Error Energy with Time Step');

% play sound

soundsc(out,SR)