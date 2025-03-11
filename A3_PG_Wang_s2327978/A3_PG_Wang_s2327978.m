%++++++++++++++++++++++++++++++++++++++++
% Moog VCF
%
% Wang s2327978
% 20 February 2024
%++++++++++++++++++++++++++++++++++++++++

clc
clear
close all

%++++++++++++++++++++++++++++++++++++++++
% input parameters
SR = 44100 ;               % sample rate [Hz]
Tf = 0.2;                  % total simulation time [s]
f0 = 120 ;                 % resonant filter frequency [Hz]
r = 0.7 ;                  % feedback coeff [choose 0 \leq r \leq 1]

% derived parameters
om0 = 2*pi*f0;             % angular resonant frequency
Nf = floor(Tf*SR);         % total number of samples for the simulation
k = 1/SR;                  % time step
a = 2^0.5*r^0.25;          % derived parameter to check stability condition 

%++++++++++++++++++++++++++++++++++++++++
%check if the stability condition for Forward Euler is satisfied
if k >= (2^0.5*a+2) / (om0*(a^2+2^0.5*a+1))...
   || k >= (2-2^0.5*a) / (om0*(a^2-2^0.5*a+1))
    error('Stability condition is violated');
end

%++++++++++++++++++++++++++++++++++++++++
% initialise the matrices and vectors 
I = eye(4);                 % 4 by 4 identity matrix
A = om0*[-1 0 0 -4*r;1 -1 0 0;0 1 -1 0;0 0 1 -1]; % 4 by 4 matrix A
b = om0*[1;0;0;0];          % 4 by 1 vector b
c = [0 0 0 1];              % 1 by 4 vector c
Bf = I + k*A;               % 4 by 4 matrix used for matrix multiplication in FE
Bb = I - k*A;               % 4 by 4 matrix used for matrix inversion or linear system solution in BE
xf = zeros(4, 1);           % initialize 4 by 1 vector for FE state
xb = zeros(4, 1);           % initialize 4 by 1 vector for BE state
yf = zeros(Nf, 1);          % initialise Nf by 1 vector to hold output y from FE
yb = zeros(Nf, 1);          % initialise Nf by 1 vector to hold output y from BE
u = [1; zeros(Nf-1, 1)];    % initialise Nf by 1 vector to hold output input sequence
tvec = (0:Nf-1)'*k;         % time vector for plots
fvec = (0:Nf-1)'*SR/Nf;     % frequency vector for plots
 
%++++++++++++++++++++++++++++++++++++++++
% main loop
tic;
for n = 1 : Nf
    % update state xf from n to n+1, with current sample of u as input (FE)
    xf = Bf*xf + k*b*u(n);
    % update state xb from n to n+1, with current sample of u as input (BE)
    xb = Bb\(xb + k*b*u(n));
% write sample to output vector yf (FE)
yf(n) = c*xf;
% write sample to output vector yb (BE)
yb(n) = c*xb;
end
simTime = toc ;
%++++++++++++++++++++++++++++++++++++++++

% shift the value of the obtained yf vector to the right as a whole
yf = circshift(yf,1);
% reset the initial value for the yf vector
yf(1)=yb(1);

% discrete transfer functions
Hf = fft(yf);
Hb = fft(yb);

% compute the exact transfer function of the continuous time system
Hc = zeros(Nf,1) ;
for n=1:Nf
    Hc(n) = c*((1i*2*pi*fvec(n)*I-A)\b);
end

% plot the magnitudes of the transfer functions 
loglog(fvec, abs(Hf), 'r', fvec, abs(Hb), 'b', fvec, abs(Hc), 'g');
title('Transfer Function Magnitudes');
xlabel('Frequency (Hz)');
ylabel('Magnitude');
legend('FE', 'BE', 'Exact');
grid on;
