%-----------------------------------------------------
% PBMMI Modal Plate and Plate Reverberation
% Junting Wang 03/22/2024
%-----------------------------------------------------

clc
clear all
close all

% sample rate
SR = 44.1e3;


%% physical parameters

Lx = 2;                   % dimensions [m] 
Ly = 1;                   % "
H = 5e-4;                 % thickness [m]
T = 700;                  % tension [N/m]
rho = 7.87e3;             % density [kg/m^3]
E = 200e9;                % youngs modulus [N/m^2]
v = 0.29;                 % poisons ratio
T60_low = 4;              % T60 for first mode
T60_high = 1;             % T60 for highest mode


%%  derived paramaters

k = 1/SR;                          % sample period
K = sqrt(E*H^2/(12*rho*(1-v^2)));  % kappa
c = sqrt(T/(rho*H));               % wavespeed
Nf = floor(T60_low * SR);          % number of samples



%% stability condition

% the maximum modal frequency
omega_max = 2/k;
% maximum modal wavenumber
beta_max = sqrt((sqrt(c^4+4*K^2*omega_max^2)-c^2)/(2*K^2)); 

Mx = floor(sqrt(beta_max^2 - (pi/Ly)^2) * Lx/pi);
My = floor(sqrt(beta_max^2 - (pi/Lx)^2) * Ly/pi);


%% create two column vectors

M = Mx * My;
[x, y] = meshgrid(1:Mx, 1:My); 
mx = reshape(x, [M, 1]);
my = reshape(y, [M, 1]);


%% calculate the modal wavenumbers

beta = sqrt((mx*pi / Lx).^2 + (my*pi / Ly).^2);

% remove problematic wavenumbers
stable = beta >= 0 & beta < beta_max;
beta = beta(stable);


%% calculate the modal frequencies

omega = sqrt(c^2*beta.^2 + K^2*beta.^4);


%% calculate the modal damping coefficients

% damping coefficients of lowest and highest mode
sigma_max = 6*log(10) / T60_low;
sigma_min = 6*log(10) / T60_high;

beta_min = sqrt((pi / Lx)^2 + (pi / Ly)^2);

% loss constants
sigma1 = (sigma_max - sigma_min) / ((beta_max)^2 - (beta_min)^2);
sigma0 = sigma_min - sigma1*beta_min^2;

% damping coefficients for each mode
sigma = sigma0 + sigma1*beta.^2;


%% positions for input and outputs

input = [0.8, 0.5];
output = [0.7, 0.6];


% calculate modal shape functions
phi_in = 2*sin(mx*pi*input(1) / Lx).*sin(my*pi*input(2) / Ly) / sqrt(Lx*Ly);
phi_out = 2*sin(mx*pi*output(1) / Lx).*sin(my*pi*output(2) / Ly) / sqrt(Lx*Ly); 
%Phi = 2 * sin(mx * pi * input(1) / Lx) .* sin(my * pi * input(2) / Ly) / sqrt(Lx * Ly);
%Phi_out = 2 * sin(mx * pi * output(1) / Lx) .* sin(my * pi * output(2) / Ly) / sqrt(Lx * Ly);

% stability
phi_in = phi_in(stable);
phi_out = phi_out(stable);


%% calculate coefficients and initialise input and output

B = (2 - omega.^2*k^2) ./ (1 + sigma*k);
C = (k*sigma - 1) ./ (1 + k*sigma);
D = k^2 ./ (1 + sigma*k);

% initialise
mode_number = length(phi_in);
p1 = zeros(mode_number, 1);
p2 = zeros(mode_number, 1);

% input
u = zeros(Nf, 1);
u(1) = 1;
num_samples = length(u);

% output
y = zeros(num_samples, 1);


%% main loop

for n = 1:num_samples
    % compute current state
    p0 = B.*p1 + C.*p2 + D.*phi_in*u(n);
    % shift states
    p2 = p1;
    p1 = p0;
    % weight sum to get audio output
    y(n) = phi_out'*p0;
end


%% normalise and differentiate

y_max = max(abs(y));
y = y / y_max;

y = diff(y);

soundsc(y, SR);

