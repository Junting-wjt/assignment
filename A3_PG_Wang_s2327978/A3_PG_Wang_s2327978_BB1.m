%++++++++++++++++++++++++++++++++++++++++
% Moog VCF BB1
%
% Wang s2327978
% 22 February 2024
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

%++++++++++++++++++++++++++++++++++++++++
% initialise the matrices and vectors 
I = eye(4);                 % 4 by 4 identity matrix
A = om0*[-1 0 0 -4*r;1 -1 0 0;0 1 -1 0;0 0 1 -1]; % 4 by 4 matrix A
b = om0*[1;0;0;0];          % 4 by 1 vector b
c = [0 0 0 1];              % 1 by 4 vector c
Lt = (I - k/2*A);           % matrix used in Trapezoidal rule for left-hand side
Rt = (I + k/2*A);           % matrix used in Trapezoidal rule for right-hand side
xt = zeros(4, 1);           % initialize 4 by 1 vector for Trapezoidal state
yt = zeros(Nf, 1);          % initialise Nf by 1 vector to hold output y from Trapezoidal
u = [1; zeros(Nf-1, 1)];    % initialise Nf by 1 vector to hold output input sequence
tvec = (0:Nf-1)'*k;         % time vector for plots
fvec = (0:Nf-1)'*SR/Nf;     % frequency vector for plots

% set the initial value for yt vector
xt = Lt \ (Rt*xt + k/2*b*(u(1) + u(2)));
yt(1) = c*xt;

%++++++++++++++++++++++++++++++++++++++++
% main loop
tic;
for n = 2:Nf-1
    % update state xt from n to n+1, with current and next sample of u as input (Trapezoidal)
    xt = Lt \ (Rt*xt + k/2*b*(u(n-1) + u(n)));
    % write sample to output vector yt (Trapezoidal)
    yt(n) = c*xt;
end
simTime = toc;
%++++++++++++++++++++++++++++++++++++++++

% discrete transfer functions
Ht = fft(yt);

% compute the exact transfer function of the continuous time system
Hc = zeros(Nf,1) ;
for n=1:Nf
    Hc(n) = c*((1i*2*pi*fvec(n)*I-A)\b);
end

% plot the magnitudes of the transfer functions 
loglog(fvec, abs(Ht), 'r', fvec, abs(Hc), 'g');
title('Transfer Function Magnitudes');
xlabel('Frequency (Hz)');
ylabel('Magnitude');
legend('Trapezoidal', 'Exact');
grid on;
