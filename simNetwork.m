% This function outputs a 'Network' structure, which can be input into
% simNetwork to run a network simulation, given an input structure of 
% network parameters, net_params, simulation parameters, sim_params, and
% synaptic conductance parameters, syn_params
% Inputs:
%   net_params - structure with fields:
%       N - total number of neurons
%       Izhi_params - structure with fields:
%           a - Nx1 vector of a Izhikevich neuron parameter values
%           b - Nx1 vector of b Izhikevich neuron parameter values
%           c - Nx1 vector of c Izhikevich neuron parameter values
%           d - Nx1 vector of d Izhikevich neuron parameter values
%       v0 - Nx1 vector of initial membrane potentials
%       I - Nx1 vector or Nxnum_time_tsteps array of current input values
%           for each cell. The latter version allows for time-varying current
%           inputs to each cell
%       cell_types - Nx1 vector with 1 or 0 for each cell, 1 corresponds to
%           excitatory cell, 0 corresponds to inhibitory cell
%       post_syn - Nx1 cell array. Each element i contains a 1xnsyn vector with
%       the indices of the post-synaptic cells to cell i
%       weights - Nx1 cell array. Each element i contains a 1xnsyn vector
%       with the weight of the connections defined in post_syn
%   sim_params - structure with fields:
%       dt - time step (ms)
%       tstop - simulation time stop (ms)
%       record_dt - recording time step (ms) (must be >= dt, usually fine to
%       set equal to dt for short simulations, e.g. couple seconds)
%   syn_params - structure with fields:
%       E_exc - reversal potential of excitatory synaptic conductance
%       E_inh - reveral potential of inhibitory synaptic conductance
%       gmax_exc - maximum conductance of excitatory synaptic conductance
%       gmax_inh - maximum conductance of inhibitory synaptic conductance
%       tau_exc - time constant of excitatory synaptic current
%       tau_inh - time constant of excitatory synaptic current
% Outputs:
%   Running without designating variable to catch output results in no
%   output and plotting data using plotData(Data), otherwise outputs:
%   Data - structure with fields:
%       t - 1xnum_time_steps time vetor
%       v - Nxnum_record_time_steps array with membrane potential
%           recordings from each neuron at each recording time 
%       spike_times - Nx1 cell array. Each element i is an 1xnspikes vector
%       of all the spike times for the ith cell    
function [varargout] = simNetwork(net_params,sim_params,syn_params)
if nargin < 1   
    N=2;
   net_params.N = N;
   net_params.Izhi_params.a = 0.02*ones(N,1);
   net_params.Izhi_params.b = 0.2*ones(N,1);
   net_params.Izhi_params.c = -65*ones(N,1);
   net_params.Izhi_params.d = 8*ones(N,1);
   net_params.v0 = -65*ones(N,1);
   net_params.I = [5; 4];
   net_params.cell_types = [ 1 ; 1 ]; % 1 - excitatory, 0 - inhibitory
   net_params.post_syn  = {[2];[1]};
   net_params.weights = {[1];[1]};     
end
if nargin < 2           
   sim_params.dt = 0.5; % (ms) simulation time step
   sim_params.tstop = 1000; % (ms) simulation end time
   sim_params.record_dt = sim_params.dt; % (ms) simulation recording time step
end
if nargin < 3
   syn_params.E_exc = 0;
   syn_params.E_inh = -80;
   syn_params.gmax_exc = 0.1; 
   syn_params.gmax_inh = 0.4;
   syn_params.tau_exc = 2; 
   syn_params.tau_inh = 7;
end
%% Initialize network 
N = net_params.N; 
dt = sim_params.dt;
tstop = sim_params.tstop;
record_dt = sim_params.record_dt;
if record_dt < dt
   record_dt = dt;
   fprintf('Recording time step set below dt, changing to dt\n'); 
end
fprintf('Simulating network of %g cells for %.2f ms\n',N,tstop);
t = 0:dt:tstop; % time vector (ms)
rt = 0:record_dt:tstop; % time vector for recording (ms)
v = net_params.v0.*ones(N,1); % initialize N x 1 cell array of instantaneous v values
u = net_params.Izhi_params.b.*v; % initialize N x 1 cell array of instantaneous u values
g_exc = zeros(N,1); % initialize excitatory synaptic conductances
g_inh = zeros(N,1); % initialize inhibitory synaptic conductances
Izhi_params = net_params.Izhi_params;
a = Izhi_params.a; b = Izhi_params.b; c = Izhi_params.c; d = Izhi_params.d;
I = net_params.I;
post_syn = net_params.post_syn;
weights = net_params.weights;
cell_types = net_params.cell_types;
E_exc = syn_params.E_exc;
E_inh = syn_params.E_inh;
% initialize recording arrays
spike_times = cell(N,1); % N x 1 cell array of spike times for each cell
v_rec = [v,zeros(N,length(rt)-1)]; % N x number of time steps array
rec_i = 2; % vm recording counter
for i = 2:length(t)    
   % Calculate synaptic currents based on previous time step
    I_exc = -1*g_exc.*(v-E_exc);
    I_inh = -1*g_inh.*(v-E_inh);
    % ODEs, forward Euler
    if size(I,2) == length(t) % AC current input
       Ii = I(:,i);
    else % DC current input for each cell
       Ii = I; 
    end
    % update the membrane potential
    v = v + dt*(0.04 * v.^2 + 5 * v + 140 - u + Ii + I_exc + I_inh);
    v(find(v>=30)) = 30;
    u = u + dt*(a.*(b.*v-u));  
    % update the synaptic conductances g_exc*(1-dt/tau)
    g_exc = g_exc * (1 - dt/syn_params.tau_exc); 
    g_inh = g_inh * (1 - dt/syn_params.tau_inh); 

    % save membrane potentials
    if t(i) == rt(rec_i)        
        v_rec(:,rec_i) = v;
        rec_i = rec_i+1; % increment next recording time
    end
    % Check cells that fired
    fired = find(v>=30); % indices of pre-synaptic cells that fired
    if any(fired)
%         fprintf('Spike at %.2f ms, cell %g\n',t(i),fired); % uncomment for troubleshooting
        num_spikes = length(fired);
        % loop through cells that fired and modify post-synaptic conductances
        for j = 1:num_spikes % for each pre-synaptic cell that fired
            pre_j = fired(j); % index of jth pre-synaptic cell that fired
            post_syn_j = post_syn{pre_j}; % all post-synaptic indices for jth pre-synaptic cell
            weights_j = weights{pre_j}; % weights of connections for jth pre-synaptic cell
            cell_type_j = cell_types(pre_j); % 1 - exc, 0 inh cell type
            spike_times{pre_j} = [spike_times{pre_j},t(i)]; % save spike time
            if cell_type_j % excitatory
                % increment by g_exc of cell's receiving APs by gmax_exc multiplied by weights
                g_exc(post_syn_j) = g_exc(post_syn_j) + syn_params.gmax_exc * weights_j;
            else % inhibitory
                % increment by g_inh of cell's receiving APs by gmax_inh multiplied by weights
                g_inh(post_syn_j) = g_inh(post_syn_j) + syn_params.gmax_inh * weights_j;
            end
        end
        % update v and u of cells that fired
        v(fired) = c(fired);
        u(fired) = u(fired) + d(fired);
    end
end
% output time vector, voltages, and spike times in Data structure 
Data.t = rt;
Data.v = v_rec;
Data.spike_times = spike_times;
if nargin == 0
   plotData(Data,net_params,1,1,1,1) 
end
if nargout == 1 % only returns Data if requesting an output
    varargout{1} = Data;
end
end