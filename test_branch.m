dt = .001;          % timestep [ms]
dx = .0001*10^-4;   % distance between nodes [cm]
L  = 10*10^-4;      % length of  entire branch [cm]
T  = .05;           % duration [ms]

N = 1:2;    % junction indices (for one branch, there are 2 nodes - one at the start and one at the end)
R = 0:1;    % rootward neighbor of each junction
NS = [0,0]; % next sibling of each junction (none because one branch)
FC = 1:-1:0;% first child of each junction

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

% v_star  = E_L;      % voltage [mV] (solution for U0(0))
% tau_n = 2;          % n-gate time constant [ms]
% tau_h = 5;          % h-gate time constant [ms]


sigma = 1;

H1 = .8; N1 = .3;

E = (g_Na*H1*E_Na+g_K*(N1^4)*E_K + g_L*E_L)/(g_Na*H1+g_K*(N1^4)+g_L);
% E = 70;

% HH paper

                    % extracellular conductivity
Vr = -65;           % rest voltage
%% 
n = 1:(L/dx+1);         % indices of nodes along branch
r =  [0 n(1:end-1)];    % rootward neighbor of each node
fc = [n(2:end) 0];      % first child of each node
ns = zeros(size(n));    % next sibling of each node 

t = 0:dt:T;            % timepoints [ms]

% voltages at each node at each timepoint
v = zeros(length(n), length(t));
v(:,1) = ones(1,length(n))*Vr;      % starting voltage

% first node is passive soma
% other nodes are part of passive dendrite


I_stim = 0;         % applied current (uA)

stim_coeff = zeros(length(n), length(t));
% stim_coeff = I_stim/(4*pi*



for tp = 2:length(t)
a = ones(size(n))*(rad/(2*rho))*(dt/(2*dx^2)); a(1) = 0;
b = ones(size(n))*(Cm + g_tot*dt/2+2*(rad/(2*rho))*(dt/(2*dx^2)));
c = -ones(size(n))*((rad/(2*rho))*(dt/(2*dx^2)));
% solving tridiagonal system
for i = max(n):(-1):1
% b_p(i) = b(i);
% a_p(i) = a(i);
    w = 0;
%     if fc(i), w = w +v(fc(i),tp-1)*(rad/(2*rho))*(dt/(2*dx^2)); end
    if r(i),  w = w + v(r(i),tp-1)*(rad/(2*rho))*(dt/(2*dx^2)); end
    w = w + v(i,tp-1)*(Cm - 2*(rad/(2*rho))*(dt/(2*dx^2))-g_tot*dt/2);
%     if fc(i), w = v(fc(i),tp-1)*a(i); end
%     if r(i),  w = w + v(r(i),tp-1)*c(i); end
%     w = w + v(i,tp-1)*b(i);
    w = w + g_tot*E*dt/2 - I_stim/(4*pi*sigma);
    
    j = fc(i);
    while j 
        b(i) = b(i) - c(j)*a(j);
        v(i,tp) = v(i,tp) + c(j)*v(j,tp);
        w = w + v(j,tp-1)*(rad/(2*rho))*(dt/(2*dx^2));
        j = ns(j);
    end
    a(i) = a(i)/b(i);
    v(i,tp) = (v(i,tp)+w)/b(i);
    
end
for i = 2:max(n)
    v(i,tp) = v(i,tp) + a(i)*v(r(i),tp);
end
end

figure, hold on
plot(t,v(1,:))
plot(t,v(2,:))
plot(t,v(100,:),'--')
