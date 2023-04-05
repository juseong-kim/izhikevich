% This function outputs figures to visualize your network activity, stored
% in the 'Data' structure and using information in the 'net_params'
% structure. You can select which plots to generate by setting the relevant
% argument to 1. If just Data and net_params are input, it plots all 4 figures
% by default. 
% Inputs:
%   Data - structure with fields:
%       t - 1xnum_time_steps time vetor
%       v - Nxnum_record_time_steps array with membrane potential
%           recordings from each neuron at each recording time 
%       spike_times - Nx1 cell array. Each element i is an 1xnspikes vector
%           of all the spike times for the ith cell   
%       net_params - structure with fields described in simNetwork, but only 
%               field used here:
%       cell_types - Nx1 vector with 1 or 0 for each cell, 1 corresponds to
%           excitatory cell, 0 corresponds to inhibitory cell
%   For rest of parameters below, set to 1 or 0 to generate plot specified
%       plot_vm - Plot of Vm vs. t for all cells.
%           When network is less than 10 cells (or a number you choose), 
%           it plots each Vm trace as a curve overlaid on the same time 
%           axes. For networks larger than 10 cells, it plots each cell's 
%           Vm at each time point as a horizontal row of pixels colored 
%           according to the membrane potential. Each row (y-axis) is a 
%           different cell, similar to a raster plot.
%       plot_raster - Raster plot for all cells. Each cell's spike times are 
%           plotted as a point, with y-coordinate corresponds to index of 
%           neuron, x-axis as time. 
%       plot_mean_vm - (top) Average Vm vs t for excitatory and inhibitory
%           cells in population, overlaid. (bottom) Fast-fourier transform
%           (FFT) of mean potentials in top to give frequency spectra for
%           both cell populations.
%       plot_firing_rate - Average firing rate per bin for excitatory and
%       inhibitory cells in population, overlaid. Change 'time_wind' to 
%       change time windows (bins) to count spikes in, i.e. 
%       firing rate(t) = (number of spikes in bin)/time_wind
function plotData(Data,net_params,plot_vm,plot_raster,plot_mean_vm,plot_firing_rate, raster_ylabels)
if nargin < 6
   plot_firing_rate = 1; % set default for firing rate plot
end
if nargin < 5
   plot_mean_vm = 1; % set default for mean Vm plot
end
if nargin < 4
    plot_raster = 1; % set default for raster plot
end
if nargin < 3
   plot_vm = 1; % set default for Vm plot
end
t = Data.t;
v = Data.v;
spike_times = Data.spike_times;
N = size(v,1); % number of neurons
%% Plot membrane potential
if plot_vm 
    figure(1); clf; 
    if N < 10 % change this to control how many cells' Vm(t) to plot on a single axis vs. as colormap
        plot(t,v); % all Vm's overlaid
        legend(num2str((1:N)','%g\n'))
        xlabel('time (ms)'); ylabel('Vm (mV)');  
        box off;
    else
        imagesc(v); % plot as rows of pixels with colors scaled by Vm
        colormap(jet); colorbar;
        xlabel('time-step'); % *x-axis is index of time vector, not time in ms*
        ylabel('Cell index');
        box off;
        ax = gca; ax.YDir = 'normal'; % 1-N going upward
    end
    title('Membrane potentials');
    if exist("raster_ylabels","var"); yticklabels(raster_ylabels); end
end
%% Plot Raster plot
if plot_raster
    figure(2); clf;
    for i = 1:N % plot point at spike time for each cell on ith row (y=i)
       plot(spike_times{i},i*ones(size(spike_times{i})),'k.'); hold on; 
    end
    xlabel('time (ms)'); ylabel('Cell index'); 
    if exist("raster_ylabels","var"); yticklabels(raster_ylabels); end
    ylim([0 N+0.5]); 
    xlim([0 t(end)]); 
    box off;
    title('Raster plot')
end
%% Plot Mean Vm of excitatory and inhibitory cells
% Can modify to plot cells corresponding to whatever property you specify
% by changing lines 87-88 to use a different condition, e.g.
% net_params.a==0.25 to average cells with a=0.25
if plot_mean_vm
   meanv_exc = mean(v(net_params.cell_types==1,:),1); 
   meanv_inh = mean(v(net_params.cell_types==0,:),1); 
   L = length(t);
   Ye = fft(detrend(meanv_exc)); Yi = fft(detrend(meanv_inh)); %fft of exc/inhib mean Vm, remove DC component
   P2e = abs(Ye/L); P2i = abs(Yi/L); % 2-sided spectrum
   P1e = P2e(1:floor(L/2)+1); P1e(2:end-1)=2*P1e(2:end-1); % 1-sided freq spectrum
   P1i = P2i(1:floor(L/2)+1); P1i(2:end-1)=2*P1i(2:end-1); % 1-sided freq spectrum
   Fe = P1e; Fi = P1i;
   freq = 1e3*(1/(t(2)-t(1)))*(0:(L/2))/L;
   figure(3); clf;
   subplot(2,1,1)
   plot(t,meanv_exc,t,meanv_inh);
   xlabel('time (ms)'); ylabel('mean Vm (mV)'); box off;
   legend('Excitatory cells','Inhibitory Cells')
   title('Mean potentials of excitatory and inhibitory cells');
   subplot(2,1,2)
   plot(freq,Fe,freq,Fi);
   xlabel('freq (Hz)'); ylabel('Amp'); 
   xlim([0 200]); box off;
   legend('Excitatory cells','Inhibitory Cells')      
   title('FFT of mean potentials of excitatory and inhibitory cells');
end
%% Plot firing rate over time in bins
% Set bin using 'time_wind'. Try out different bin sizes to see what
% happens to the firing rates over time
if plot_firing_rate
    time_wind = 100; % time window for binary, using spikes after LFMS pulse
    time_vec_wind = (0:time_wind:t(end));
    spiketimes_bin = cellfun(@(x) (histc(x,time_vec_wind)),spike_times,'UniformOutput',false); % convert to binary trains using time window above           
    for i = 1:N
       if isempty(spiketimes_bin{i})
           spiketimes_bin{i} = zeros(1,length(time_vec_wind)); 
       end
    end
    spiketimes_bin = cell2mat(spiketimes_bin); 
    firing_rates_exc = mean(spiketimes_bin(net_params.cell_types==1,:),1)/(time_wind*1e-3);
    firing_rates_inh = mean(spiketimes_bin(net_params.cell_types==0,:),1)/(time_wind*1e-3);
    figure(4); 
    plot(time_vec_wind,firing_rates_exc,time_vec_wind,firing_rates_inh)
    legend('Excitatory cells','Inhibitory Cells');
    xlabel('time (ms)'); 
    ylabel('Firing rate');
    title(sprintf('Firing rate over time (bin=%.1f ms)',time_wind))
    box off;
end
end