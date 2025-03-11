%-----------------------------------------------------
% PBMMI Modal Plate and Plate Reverberation - BB1
%
% -apply reverb to custom file
% -allow stereo input and output
% -allow users to pick materials
%
% Junting Wang 03/30/2024
%-----------------------------------------------------

clc
clear all
close all

% sample rate
SR = 44.1e3;
% pick materials
choose_materials = 1;
% choose input
choose_input = 2;


%% physical parameters

switch choose_materials

% steel
    case 1
Lx = 2;                   % dimensions [m] 
Ly = 1;                   % "
H = 5e-4;                 % thickness [m]
T = 700;                  % tension [N/m]
rho = 7.87e3;             % density [kg/m^3]
E = 200e9;                % youngs modulus [N/m^2]
v = 0.29;                 % poisons ratio
T60_low = 4;              % T60 for first mode
T60_high = 1;             % T60 for highest mode

% glass
    case 2
Lx = 2;                   % dimensions [m] 
Ly = 1;                   % "
H = 5e-4;                 % thickness [m]
T = 700;                  % tension [N/m]
rho = 2.5e3;              % density [kg/m^3]
E = 70e9;                 % youngs modulus [N/m^2]
v = 0.25;                 % poisons ratio
T60_low = 4;              % T60 for first mode
T60_high = 1;             % T60 for highest mode
end


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

sigma_max = 6*log(10) / T60_low;
sigma_min = 6*log(10) / T60_high;
beta_min = sqrt((pi / Lx)^2 + (pi / Ly)^2);
sigma1 = (sigma_max - sigma_min) / ((beta_max)^2 - (beta_min)^2);
sigma0 = sigma_min - sigma1*beta_min^2;
sigma = sigma0 + sigma1*beta.^2;


%% positions for input and outputs

input = [0.8, 0.5];
output_left = [0.7, 0.6];
output_right = [0.5, 0.2]; 

% calculate modal shape functions
phi_in = 2*sin(mx*pi*input(1) / Lx).*sin(my*pi*input(2) / Ly) / sqrt(Lx*Ly);
phi_left = 2*sin(mx*pi*output_left(1) / Lx).*sin(my*pi*output_left(2) / Ly) / sqrt(Lx*Ly); 
Phi_right = 2 * sin(mx * pi * output_right(1) / Lx) .* sin(my * pi * output_right(2) / Ly) / sqrt(Lx * Ly);

% stability
phi_in = phi_in(stable);
phi_left = phi_left(stable);
Phi_right = Phi_right(stable);


%% calculate coefficients and initialise input and output

B = (2 - omega.^2*k^2) ./ (1 + sigma*k);
C = (k*sigma - 1) ./ (1 + k*sigma);
D = k^2 ./ (1 + sigma*k);

mode_number = length(phi_in);
p1 = zeros(mode_number, 1);
p2 = zeros(mode_number, 1);

% input
switch choose_input

    % impulse
    case 1
        u = zeros(Nf, 1);
        u(1) = 1;

    % audio file
    case 2
        [u, SR] = audioread("Godin4_44.wav");
        % average left and right channels to mono          
        [r_number, c_number] = size(x);
        if c_number == 2
            x = (x(:,1) + x(:,2)) / 2;                   
        end
end

num_samples = length(u);

% output
y = zeros(num_samples, 2);


%% main loop

for n = 1:num_samples
    % compute current state
    p0 = B.*p1 + C.*p2 + D.*phi_in*u(n);
    % shift states
    p2 = p1;
    p1 = p0;
    % weight sum to get audio output
    y(n,1) = phi_left'*p0;
    y(n,2) = Phi_right'*p0;
end


%% normalise and differentiate

y_max = max(max(abs(y)));
y = y / y_max;
y = diff(y);
soundsc(y, SR);

