%-----------------------------------------------------
% PBMMI stiff string assignment template - BB1
% find parameters for acoustic guitar
% Junting Wang 03/15/2024
%-----------------------------------------------------

clc
clear all
close all

%%%%% flags

plot_on = 0;                % in-loop plotting on (1) or off (0)
itype = 1;                  % type of input: 1: pluck, 2: strike
string_num = 6;             % choose different notes

%%%%% parameters

% physical string parameters

switch string_num
    % E
    case 1
        T = 90;                     % tension (N)
        r = 0.00083;                % string radius (m)
        E = 2e11;                   % Young's modulus (Pa) for steel
        rho = 7850;                 % density (kg/m^3)
        T60 = 3;                    % T60 (s)
        L = 0.65;                   % length (m)
    % A
    case 2 
        T = 87;                     % tension (N)
        r = 0.00060;                % string radius (m)
        E = 2e11;                   % Young's modulus (Pa) for steel
        rho = 7850;                 % density (kg/m^3)
        T60 = 3;                    % T60 (s)
        L = 0.65;                   % length (m)
    % D
    case 3 
        T = 80;                     % tension (N)
        r = 0.00043;                % string radius (m)
        E = 2e11;                   % Young's modulus (Pa) for steel
        rho = 7850;                 % density (kg/m^3)
        T60 = 3;                    % T60 (s)
        L = 0.65;                   % length (m)
    % G
    case 4
        T = 74;                     % tension (N)
        r = 0.00031;                % string radius (m)
        E = 2e11;                   % Young's modulus (Pa) for steel
        rho = 7850;                 % density (kg/m^3)
        T60 = 3;                    % T60 (s)
        L = 0.65;                   % length (m)
    % B
    case 5
        T = 71;                     % tension (N)
        r = 0.00024;                % string radius (m)
        E = 2e11;                   % Young's modulus (Pa) for steel
        rho = 7850;                 % density (kg/m^3)
        T60 = 3;                    % T60 (s)
        L = 0.65;                   % length (m)
    % E
    case 6
        T = 70;                     % tension (N)
        r = 0.00018;                % string radius (m)
        E = 2e11;                   % Young's modulus (Pa) for steel
        rho = 7850;                 % density (kg/m^3)
        T60 = 3;                    % T60 (s)
        L = 0.65;                   % length (m)
end

% I/O

SR = 44100;                 % sample rate (Hz)
Tf = 5;                     % duration of simulation (s)

xi = 0.9;                   % coordinate of excitation (normalised, 0-1)
famp = 5;                   % peak amplitude of excitation (N)
dur = 0.001;                % duration of excitation (s)
exc_st = 0.001;             % start time of excitation (s)

xo = 0.1;                   % coordinate of output (normalised, 0-1)

%% checking on all these parameters

if T<0 || r<0 || E<0 || rho<0 || T60<0 || L<0 || SR<0 || Tf<0 || xi<0 ...
    || xi>1 || famp<0 || dur<0 || exc_st<0 || xo<0 || xo>1
    error ('all the parameters need to be non-negativity,' + ...
        'xi and xo should between 0 and 1');
end

if exc_st + dur > Tf
    error('excitation should between duration of simulation')
end

%% derived parameters

A = pi*r^2;                 % string cross-sectional area
I = 0.25*pi*r^4;            % string moment of inertia

c = sqrt(T/(rho*A));        % wave speed
K = sqrt(E*I/(rho*A));      % stiffness constant 
sig = 6*log(10)/T60;        % frequency independent loss parameter

k = 1/SR;                   % time step
Nf = floor(SR*Tf);          % number of time steps

%%

% determine a minimal grid spacing hmin
hmin = sqrt((c^2*k^2 + (sqrt((c^4*k^4) + (16*K^2*k^2))))/2);

% determine the number of grid spacings N
N = floor(L/hmin);

% reset h
h = L/N;


%% 

% ensure the number of segments, N, is less than 10000.
if N > 10000
    error('N should less than 10000');
end

if xi*L < h || (1-xi)*L < h || xo*L < h || (1-xo)*L < h
    error('xo or xi should be at least h metres away from either endpoint of the string');
end

%% representing the input force in Newtons. 

t_input = (0:Nf-1) / SR;         % time vector of force
f = zeros(1, Nf);                % initialise force vector

for n = 1:Nf
% check the input falls within the time span
   if t_input(n) >= exc_st && t_input(n) <= exc_st + dur 
       f(n) = 1/2 * famp * (1 - cos(itype*pi*(t_input(n) - exc_st) / dur));
   end
end

%% create matrix

I_mtr = speye(N-1);                % identity matrix
i = ones(N-1,1);

% create Dxx
Dxx_v = [i -2*i i];                % diagonal value
Dxx_p = [-1 0 1];                  % diagonal position
Dxx = (1/h^2)*spdiags(Dxx_v, Dxx_p, N-1, N-1);

% create Dxxxx
Dxxxx_v = [i -4*i 6*i -4*i i];     % diagonal value
Dxxxx_p = [-2 -1 0 1 2];           % diagonal value
Dxxxx = spdiags(Dxxxx_v, Dxxxx_p, N-1, N-1);
Dxxxx(1, 1) = 5;
Dxxxx(N - 1, N - 1) = 5;
Dxxxx = (1/h^4) * Dxxxx;

% approximated the Dirac delta using an indicator function J
li = floor(L*xi/h);
lo = floor(L*xo/h);
J = zeros(N-1,1);
J(li) = 1;

B = 1/(1+sig*k) * (2*I_mtr + c^2*k^2*Dxx - K^2*k^2*Dxxxx);
C =(-1+sig*k) / (1+sig*k) * I_mtr;
J = k^2/(rho*A*h) * J;

c_sel = zeros(N-1,1);              % initialize vector for selecting
c_sel(lo) = 1;                     % vector selecting output location

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% initialise scheme variables

u2 = zeros(N-1,1);              % state
u1 = u2;                        % state
u = u2;                         % state

y = zeros(Nf,1);                % output
xax = (1:N-1)'*h;               % x-axis for plotting

%% main loop

tic
for n=1:Nf
    
    % update state, and insert current value of f.
    
    u = B*u1 + C*u2 + J*f(n);
    
    % read output
    y(n) = c_sel'*u;
    
    % plot
    
    if(plot_on==1)
        % draw current state
        figure(1)
        plot(xax, u, 'k');
        xlabel('x (m)')
        ylabel('u (m)')
        axis([0 L -0.005 0.005])
        drawnow
    end
    
    % shift state
    
    u2 = u1;
    u1 = u;
    
end
toc

%% play sound

soundsc(y,SR);

%% plot spectrum

figure(2)

yfft = 10*log10(abs(fft(y)));
plot([0:Nf-1]'*(SR/Nf), yfft, 'k')
xlim([0 SR/2])
xlabel('Frequency (Hz)')
ylabel('Magnitude')
title('Transform of output')