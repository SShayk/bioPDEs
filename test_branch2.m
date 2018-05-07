dt = .001;          % timestep [ms]
dx = .001;          % distance between nodes [cm]
L = .01;            % length of branch [cm]
T = .05;            % duration [ms]

N = 1:2;            % junction indices
R = 0:1;            % junction rootward neighbors
FC = 1:-1:0;        % junction first children
NS = [0 0];         % junction next siblings

% Hodgkin-Huxley model parameters from Griffith and Peskin, "Electrophysiology",
% Communications on Pure and Applied Mathematics (2013)
Cm      = 1;        % capacitance per unit area [uF/cm^2]
g_Na    = 120;      % conductance (Na+ channel) [(uA/mV)/cm^2]
g_K     = 36;       % conductance (K+ channel) [(uA/mV)/cm^2]
g_L     = 0.3;      % conductance (leak) [(uA/mV)/cm^2]
g_tot = g_Na + g_K + g_L;

E_Na    = 45;       % rest potential (Na+ channel) [mV]
E_K     = -82;      % rest potential (K+ channel) [mV]
E_L     = -59;      % rest potential (leak channel) [mV]
rad     = .0238;    % axon radius [cm]
rho     = .0354;    % electrical resistivity [(mV/uA)cm]

E = -70;
Vr = -70;

n = 1:(L/dx+1);     % node indices
r = [0 n(1:end-1)]; % node rootward neighbors
fc = [n(2:end) 0];  % node first children
ns = zeros(size(n));% node next siblings

t = 0:dt:T;                     % time vector
v = zeros(length(n),length(t)); % voltages at each point in space/time

%% 

%begin at resting voltage
v(:,1) = Vr; 

% clamp ends
v(1,:) = Vr;
v(end,:) = Vr;



I_stim = 0;

for tp = 2:length(t)


end