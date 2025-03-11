%-----------------------------------------------------
% PBMMI State-space VA Modelling
% Junting Wang 04/06/2024
%-----------------------------------------------------

clc
clear all
close all

% Simulation parameters
SR = 88.2e3;               % sample rate [Hz]
dur = 2;                   % duration of sim [s]
Nf = round(SR*dur);        % number of samples in sim
k = 1 / SR;                % time step

% Settings
plotting = true;           % plotting on/off
audio = true;              % play output sounds on/off

% Input - create sine wave input
f0 = 220;                  % frequency [Hz]
amp = 5;                   % amplitude of input
tvec = dur*(0:Nf-1)'/Nf;   % time vector [s]
u = amp*sin(2*pi*f0*tvec); % input vector

if audio
    % read in the audio file
    [input, SR] = audioread('A6_PG_Wang_s2327978_Input.wav');

    % average left and right channels to mono
    [r_number, c_number] = size(input); 
    if c_number == 2
       u = (input(:,1) + input(:,2)) / 2; 

    else
       u = input;
    end

    u = amp/max(abs(u))*u;
    Nf = length(u);
    dur = round(Nf/SR);
    tvec = dur*(0:Nf-1)'/Nf;
end

% Physical parameters
r = 1e3;                   % resistance of R [Ohms]
c = 33e-9;                 % capacitance of C [F]
Is = 2.52e-9;              % diode saturation current (A)
Vt = 25.3e-3;              % diode thermal voltage (V)
Ni = 1.752;                % diode ideality factor

% Newton-Raphson parameters
tol = 1e-9;                % convergence tolerance
max_iter = 50;             % maximum allowed iterations per time-step

% Simulate the diode clipper system
I = eye(1);
A = -1/(r*c);
B = 1/(r*c);
C = -1/c;

D = 1;
E = 0;
F = 0;

L = 1;
M = 0;
N = 0;

% Derived parameter
H_minus = 2/k * I - A;
H_plus = 2/k * I + A;
K = D * H_minus \ C + F;

% Initialize variables
y = zeros(Nf, 1);
x = 0;       % x(n+1)

x1 = 0;      % x(n)
f1 = 0;      % f(n)

f = zeros(Nf, 1);
v = zeros(Nf, 1);
u1 = 0;

tic
% main loop
for n = 1 : Nf
    
    vn = v(n);
    v1 = D*x1 + E*u1 + F*f1;
    f1 = 2 * Is * sinh(v1/(Ni*Vt));
  
    % calculate p(n)
    p = D*H_minus\(H_plus*x1+B*u1+C*f1) + (D*H_minus\B+E)*u(n);


    % find v(n+1)
    step = 5;  
    iters = 0;
   
    while(abs(step) > tol)&&(iters < max_iter)

        g = p + K *  2*Is*sinh(vn/(Ni*Vt)) - vn;
        j = K*2*Is*cosh(vn/(Ni*Vt))/(Ni*Vt)-1;

        step = g/j;
    
        vn = vn-step;
        iters = iters+1;
    end


    % compute state update
    f =  2*Is*sinh(vn/(Ni*Vt));    % f(n+1)
    x = H_minus \ (H_plus*x1 + B*(u1+u(n)) + C*(f1+f)); % x(n+1)

    % compute output
    y(n) = L*x1 +M*u1 + N*f;  

    % update state
    x1 = x;
    u1 = u(n);
end
toc


if audio
    soundsc(y,SR);     
    audiowrite("A6_PG_Wang_s2327978_Output.wav", y, SR);
end

if plotting
    % plot the input and output in time domain
    subplot(2,1,1);
    plot(tvec, u, 'b', tvec, y, 'r');
    title('Input and Output Signals');
    xlabel('Time [s]');
    ylabel('Voltage [V]');
    legend('Input', 'Output');

    % plot the frequency domain
    subplot(2,1,2);
    Y = fft(y)/Nf;
    frequency = SR*(0:(Nf/2))/Nf;
    plot(frequency, 20*log10(2*abs(Y(1:Nf/2 +1))));
    title('Frequency Domain Plot');
    xlabel('Frequency(Hz)');
    ylabel('dB');
end

