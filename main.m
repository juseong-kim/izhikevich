%% 1.A Regular Spiking (RS) Neuron
close all
clear
% Length
figure
t = 0:0.1:200;
[v0, u0] = deal(-65, -65 * 0.2);
[a, b, c, d] = deal(0.02, 0.2, -65, 8);
I = [zeros(1,50), 10 * ones(1,size(t,2)-50)];
[v, u] = izhikevich(t,v0,u0,I,a,b,c,d);
plot(t,v,'DisplayName','v')
hold on
plot(t,u,'DisplayName','u')
hold on
plot(t,I,'DisplayName','I')
legend()
title('RS Neuron - Izhikevich Model',"a = "+a+", b = "+b+", c = "+c+", d = "+d)
xlabel('Time (ms)')
ylabel('Membrane Potential (mV)')
hold off
saveas(gcf,'Figures/1A_RS.png')

%% 1.B Low-threshold (LTS) Neuron
% Length
figure
t = 0:0.1:200;
[v0, u0] = deal(-65, -65*0.2);
[a, b, c, d] = deal(0.02, 0.25, -65, 2);
I = [zeros(1,50), 10 * ones(1,size(t,2)-50)];
[v, u] = izhikevich(t,v0,u0,I,a,b,c,d);
plot(t,v,'DisplayName','v')
hold on
plot(t,u,'DisplayName','u')
hold on
plot(t,I,'DisplayName','I')
legend()
title('LTS Neuron - Izhikevich Model',"a = "+a+", b = "+b+", c = "+c+", d = "+d)
xlabel('Time (ms)')
ylabel('Membrane Potential (mV)')
hold off
saveas(gcf,'Figures/1B_LTS.png')

%% 1.C TC Neuron in Both Firing Regimes
figure
t = 0:0.1:100;
[a, b, c, d] = deal(0.02, 0.25, -65, 0.05);

% Depolarization
[v0_d, u0_d] = deal(-63, -63*b);
I_d = [zeros(1,400), 10 * ones(1,size(t,2)-400)];
[v_d, u_d] = izhikevich(t,v0_d,u0_d,I_d,a,b,c,d);

% Hyperpolarization
[v0_h, u0_h] = deal(-87, -87*b);
I_h = [-10 * ones(1,400), zeros(1,size(t,2)-400)];
[v_h, u_h] = izhikevich(t,v0_h,u0_h,I_h,a,b,c,d);

% Plot
plot(t,v_d,'DisplayName','v_{depolarization}')
hold on
plot(t,v_h,'DisplayName','v_{hyperpolarization}')
hold on
plot(t,I_d,'DisplayName','I_{depolarization}')
hold on
plot(t,I_h,'DisplayName','I_{hyperpolarization}')
legend()
title('TC Neuron - Izhikevich Model',"a = "+a+", b = "+b+", c = "+c+", d = "+d)
xlabel('Time (ms)')
ylabel('Membrane Potential (mV)')
hold off
saveas(gcf,'Figures/1C_TC.png')


%% 2.A Two RS Neurons Connected - Excitatory / Inhibitory
close all
figure
t = 0:0.5:500;
[a, b, c, d] = deal(0.02, 0.2, -65, 8);
[v0, u0, i0, g0] = deal(-65, -65 * 0.2, 0, 0);

ei = 0; % 1 - excitatory, 0 = inhibitory
if ei
    [tau,gbar,vr] = deal(2,0.1,0); % Excitatory
else
    [tau,gbar,vr] = deal(7,0.4,-80); % Inhibitory
end

% ddelta = [zeros(1,149) 1 zeros(1,size(t,2)-150)];
% [v1, v2, I1, I2, u1, u2, g1, g2] = izhikevich_two(t,v0,u0,i0,a,b,c,d,tau,gbar,ddelta,vr);
[v1, v2, I1, I2, u1, u2, g1, g2] = izhi_two(t,v0,u0,i0,g0,a,b,c,d,tau,gbar,vr,ei);

plot(t,v1,'DisplayName','v1')
hold on
plot(t,v2,'DisplayName','v2')
hold on
plot(t,I1,'DisplayName','I1')
hold on
plot(t,I2,'DisplayName','I2')
hold on
plot(t,g1,'DisplayName','g1')
hold on
plot(t,g2,'DisplayName','g2')
legend()
if ei; ei_label = "Excitatory"; else; ei_label = "Inhibitory"; end
title('Connected RS Neurons - '+ei_label,"a = "+a+", b = "+b+", c = "+c+", d = "+d+", V0 = "+v0+", gbar = "+gbar)
xlabel('Time (ms)')
ylabel('Membrane Potential (mV)')
hold off
saveas(gcf,'Figures/2A_'+ei_label+'.png')


%% 3.A Large Network
close all
% Network Parameters
N = 1000;
Ne = 0.8 * N;       Ni = 0.2 * N;
re = rand(Ne, 1);   ri = rand(Ni, 1);
net_params.N = N;
net_params.Izhi_params.a = [0.02*ones(Ne,1);    0.02+0.08*ri];
net_params.Izhi_params.b = [0.2*ones(Ne,1);     0.25-0.05*ri];
net_params.Izhi_params.c = [-65+15*re.^2;       -65*ones(Ni,1)];
net_params.Izhi_params.d = [8-6*re.^2;          2*ones(Ni,1)];
net_params.v0 = -65*ones(N,1);
net_params.I = [5*randn(Ne,1); 2*randn(Ni,1)];
net_params.cell_types = [ ones(Ne,1); zeros(Ni,1) ]; % 1 - excitatory, 0 - inhibitory
net_params.post_syn = mat2cell(randi(N,N),N,ones(N,1));
net_params.weights = mat2cell(rand(N,N),N,ones(N,1)); % num2cell([0.5*rand(N,Ne) -rand(N,Ni)]);

% Simulation Parameters
sim_params.dt = 0.1; % (ms) simulation time step
sim_params.tstop = 1000; % (ms) simulation end time
sim_params.record_dt = sim_params.dt; % (ms) simulation recording time step

% Synaptic Parameters
syn_params.E_exc = 0;
syn_params.E_inh = -80;
syn_params.gmax_exc = 0.1;
syn_params.gmax_inh = 0.4; % 0.4
syn_params.tau_exc = 2; % 2
syn_params.tau_inh = 7;

% Simulate
Data = simNetwork(net_params, sim_params, syn_params);
plotData(Data, net_params, 1, 1, 1, 1);
saveas(figure(1),strcat('Figures/3a/vm',int2str(N),'.png'));
saveas(figure(2),strcat('Figures/3a/raster',int2str(N),'.png'));
saveas(figure(3),strcat('Figures/3a/meanvm',int2str(N),'.png'));
saveas(figure(4),strcat('Figures/3a/firingrate',int2str(N),'.png'));

%% 3.B Large Network - Synaptic Weights E > I
close all
dweight = 1;
% Network Parameters
N = 100;
Ne = 0.8 * N;       Ni = 0.2 * N;
re = rand(Ne, 1);   ri = rand(Ni, 1);
net_params.N = N;
net_params.Izhi_params.a = [0.02*ones(Ne,1);    0.02+0.08*ri];
net_params.Izhi_params.b = [0.2*ones(Ne,1);     0.25-0.05*ri];
net_params.Izhi_params.c = [-65+15*re.^2;       -65*ones(Ni,1)];
net_params.Izhi_params.d = [8-6*re.^2;          2*ones(Ni,1)];
net_params.v0 = -65*ones(N,1);
net_params.I = [5*randn(Ne,1); 2*randn(Ni,1)];
net_params.cell_types = [ ones(Ne,1); zeros(Ni,1) ]; % 1 - excitatory, 0 - inhibitory
net_params.post_syn = mat2cell(randi(N,N),N,ones(N,1));
net_params.weights = mat2cell([dweight+0.1*rand(N,Ne) 0.1*rand(N,Ni)],N,ones(N,1));

% Simulation Parameters
sim_params.dt = 0.1; % (ms) simulation time step
sim_params.tstop = 1000; % (ms) simulation end time
sim_params.record_dt = sim_params.dt; % (ms) simulation recording time step

% Synaptic Parameters
syn_params.E_exc = 0;
syn_params.E_inh = -80;
syn_params.gmax_exc = 0.1;
syn_params.gmax_inh = 0.4; % 0.4
syn_params.tau_exc = 2; % 2
syn_params.tau_inh = 7;

% Simulate
Data = simNetwork(net_params, sim_params, syn_params);
plotData(Data, net_params, 1, 1, 1, 1);
saveas(figure(1),strcat('Figures/3b/vm',int2str(N),"_dw",num2str(dweight),'.png'));
saveas(figure(2),strcat('Figures/3b/raster',int2str(N),"_dw",num2str(dweight),'.png'));
saveas(figure(3),strcat('Figures/3b/meanvm',int2str(N),"_dw",num2str(dweight),'.png'));
saveas(figure(4),strcat('Figures/3b/firingrate',int2str(N),"_dw",num2str(dweight),'.png'));

%% 3.C Large Network - Increased Network Size
close all
% Network Parameters
N = 10000;
Ne = 0.8 * N;       Ni = 0.2 * N;
re = rand(Ne, 1);   ri = rand(Ni, 1);
net_params.N = N;
net_params.Izhi_params.a = [0.02*ones(Ne,1);    0.02+0.08*ri];
net_params.Izhi_params.b = [0.2*ones(Ne,1);     0.25-0.05*ri];
net_params.Izhi_params.c = [-65+15*re.^2;       -65*ones(Ni,1)];
net_params.Izhi_params.d = [8-6*re.^2;          2*ones(Ni,1)];
net_params.v0 = -65*ones(N,1);
net_params.I = [5*randn(Ne,1); 2*randn(Ni,1)];
net_params.cell_types = [ ones(Ne,1); zeros(Ni,1) ]; % 1 - excitatory, 0 - inhibitory
net_params.post_syn = mat2cell(randi(N,N),N,ones(N,1));
net_params.weights = mat2cell(rand(N,N),N,ones(N,1));

% Simulation Parameters
sim_params.dt = 0.1; % (ms) simulation time step
sim_params.tstop = 1000; % (ms) simulation end time
sim_params.record_dt = sim_params.dt; % (ms) simulation recording time step

% Synaptic Parameters
syn_params.E_exc = 0;
syn_params.E_inh = -80;
syn_params.gmax_exc = 0.1;
syn_params.gmax_inh = 0.4; % 0.4
syn_params.tau_exc = 2; % 2
syn_params.tau_inh = 7;

% Simulate
Data = simNetwork(net_params, sim_params, syn_params);
plotData(Data, net_params, 1, 1, 1, 1);
saveas(figure(1),strcat('Figures/3c/vm',int2str(N),'.png'));
saveas(figure(2),strcat('Figures/3c/raster',int2str(N),'.png'));
saveas(figure(3),strcat('Figures/3c/meanvm',int2str(N),'.png'));
saveas(figure(4),strcat('Figures/3c/firingrate',int2str(N),'.png'));

%% 3.D Large Network - Increased Inhibitory Synaptic Current
close all
% Network Parameters
N = 1000;
Ne = 0.8 * N;       Ni = 0.2 * N;
re = rand(Ne, 1);   ri = rand(Ni, 1);
net_params.N = N;
net_params.Izhi_params.a = [0.02*ones(Ne,1);    0.02+0.08*ri];
net_params.Izhi_params.b = [0.2*ones(Ne,1);     0.25-0.05*ri];
net_params.Izhi_params.c = [-65+15*re.^2;       -65*ones(Ni,1)];
net_params.Izhi_params.d = [8-6*re.^2;          2*ones(Ni,1)];
net_params.v0 = -65*ones(N,1);
net_params.I = [5*randn(Ne,1); 2*randn(Ni,1)];
net_params.cell_types = [ ones(Ne,1); zeros(Ni,1) ]; % 1 - excitatory, 0 - inhibitory
net_params.post_syn = mat2cell(randi(N,N),N,ones(N,1));
net_params.weights = mat2cell(rand(N,N),N,ones(N,1));

% Simulation Parameters
sim_params.dt = 0.1; % (ms) simulation time step
sim_params.tstop = 1000; % (ms) simulation end time
sim_params.record_dt = sim_params.dt; % (ms) simulation recording time step

% Synaptic Parameters
syn_params.E_exc = 0;
syn_params.E_inh = -80;
syn_params.gmax_exc = 0.1;
syn_params.gmax_inh = 0.4; % 0.4
syn_params.tau_exc = 2; % 2
syn_params.tau_inh = 10;

% Simulate
Data = simNetwork(net_params, sim_params, syn_params);
plotData(Data, net_params, 1, 1, 1, 1);
saveas(figure(1),strcat('Figures/3d/vm',int2str(N),'_t',int2str(syn_params.tau_inh),'.png'));
saveas(figure(2),strcat('Figures/3d/raster',int2str(N),'_t',int2str(syn_params.tau_inh),'.png'));
saveas(figure(3),strcat('Figures/3d/meanvm',int2str(N),'_t',int2str(syn_params.tau_inh),'.png'));
saveas(figure(4),strcat('Figures/3d/firingrate',int2str(N),'_t',int2str(syn_params.tau_inh),'.png'));

%% 4.A CPG
close all

% Simulation Parameters
sim_params.dt = 0.5; % (ms) simulation time step
sim_params.tstop = 1500; % (ms) simulation end time
sim_params.record_dt = sim_params.dt; % (ms) simulation recording time step


% currents for M, C, L, E, R
M = 0; % 4      x
C = 4.5; % 5      3.5
L = 3.5; % 0.5    0.35
E = 4.5; % 1.25      0.875
R = 5.5; % 5      3.5

% Network Parameters
N = 10;
net_params.N = N;
net_params.Izhi_params.a = 0.02*ones(N,1);
net_params.Izhi_params.b = 0.2*ones(N,1);
net_params.Izhi_params.c = -65*ones(N,1);
net_params.Izhi_params.d = 8*ones(N,1);
net_params.v0 = -65*ones(N,1);
net_params.I = cpg_current(sim_params.tstop,sim_params.dt,100, ...
                            M,M, ... % M
                            C,C, ... % C
                            L,L, ... % L
                            E,E, ... % E
                            R,R);    % R
net_params.cell_types = [1 1 0 0 0 0 1 1 1 1]'; % 1 - excitatory, 0 - inhibitory
LM = 1; RM = 2;
LC = 3; RC = 4; 
LL = 5; RL = 6;
LE = 7; RE = 8;
LR = 9; RR = 10;
net_params.post_syn = {[];[]; % M (L/R)
                        [RE RL RC RR LM];[LE LL LC LR RM]; % C (L/R)
                        [LC];[RC]; % L (L/R)
                        [LL LC LR LM];[RL RC RR RM]; % E (L/R)
                        [LE LL LC];[RE RL RC]; % R (L/R)
                        };
net_params.weights = {[];[]; % M (L/R)
                        [3.5 3.5 3.5 5.5 5.5]';[3.5 3.5 3.5 5.5 5.5]'; % C (L/R)
                        [3.5];[3.5]; % L (L/R)
                        [3.5 3.5 5.5 5.5]';[3.5 3.5 5.5 5.5]'; % E (L/R)
                        [1.5 0.55 0.7]';[1.5 0.55 0.7]'; % R (L/R)
                        };

% Synaptic Parameters
syn_params.E_exc = 0;
syn_params.E_inh = -80;
syn_params.gmax_exc = 0.1; 
syn_params.gmax_inh = 0.4; 
syn_params.tau_exc = 2;
syn_params.tau_inh = 7;

% Simulate
Data = simNetwork(net_params, sim_params, syn_params);
raster_ylabels = ["left M" "right M" "left C" "right C" "left L"...
                    "right L" "left E" "right E" "left R" "right R"];
plotData(Data, net_params, 1,1,1,1, raster_ylabels);
saveas(figure(1),strcat('Figures/4a/vm','.png'));
saveas(figure(2),strcat('Figures/4a/raster','.png'));
saveas(figure(3),strcat('Figures/4a/meanvm','.png'));
saveas(figure(4),strcat('Figures/4a/firingrate','.png'));

%% 4.B CPG Increased Oscillatory Frequency
close all

% Simulation Parameters
sim_params.dt = 0.5; % (ms) simulation time step
sim_params.tstop = 1500; % (ms) simulation end time
sim_params.record_dt = sim_params.dt; % (ms) simulation recording time step


% currents for M, C, L, E, R
M = 0; % 4      x
C = 6; % 5      3.5
L = 5; % 0.5    0.35
E = 6; % 1.25      0.875
R = 7; % 5      3.5

% Network Parameters
N = 10;
net_params.N = N;
net_params.Izhi_params.a = 0.02*ones(N,1);
net_params.Izhi_params.b = 0.2*ones(N,1);
net_params.Izhi_params.c = -65*ones(N,1);
net_params.Izhi_params.d = 8*ones(N,1);
net_params.v0 = -65*ones(N,1);
net_params.I = cpg_current(sim_params.tstop,sim_params.dt,100, ...
                            M,M, ... % M
                            C,C, ... % C
                            L,L, ... % L
                            E,E, ... % E
                            R,R);    % R
net_params.cell_types = [1 1 0 0 0 0 1 1 1 1]'; % 1 - excitatory, 0 - inhibitory
LM = 1; RM = 2;
LC = 3; RC = 4; 
LL = 5; RL = 6;
LE = 7; RE = 8;
LR = 9; RR = 10;
net_params.post_syn = {[];[]; % M (L/R)
                        [RE RL RC RR LM];[LE LL LC LR RM]; % C (L/R)
                        [LC];[RC]; % L (L/R)
                        [LL LC LR LM];[RL RC RR RM]; % E (L/R)
                        [LE LL LC];[RE RL RC]; % R (L/R)
                        };
net_params.weights = {[];[]; % M (L/R)
                        [3.5 3.5 3.5 5.5 5.5]';[3.5 3.5 3.5 5.5 5.5]'; % C (L/R)
                        [3.5];[3.5]; % L (L/R)
                        [3.5 3.5 5.5 5.5]';[3.5 3.5 5.5 5.5]'; % E (L/R)
                        [1.5 0.55 0.7]';[1.5 0.55 0.7]'; % R (L/R)
                        };

% Synaptic Parameters
syn_params.E_exc = 0;
syn_params.E_inh = -80;
syn_params.gmax_exc = 0.1; 
syn_params.gmax_inh = 0.4; 
syn_params.tau_exc = 2;
syn_params.tau_inh = 7;

% Simulate
Data = simNetwork(net_params, sim_params, syn_params);
raster_ylabels = ["left M" "right M" "left C" "right C" "left L"...
                    "right L" "left E" "right E" "left R" "right R"];
plotData(Data, net_params, 1,1,1,1, raster_ylabels);
saveas(figure(1),strcat('Figures/4b/vm','.png'));  % num2str(M),"_",num2str(C),"_",num2str(L),"_",num2str(E),"_",num2str(R)
saveas(figure(2),strcat('Figures/4b/raster','.png'));
saveas(figure(3),strcat('Figures/4b/meanvm','.png'));
saveas(figure(4),strcat('Figures/4b/firingrate','.png'));


