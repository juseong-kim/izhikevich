function [v, u] = izhikevich(t,v0,u0,I,a,b,c,d)
    % use default values if not provided
    % a = scale of u
    % b = sensitivity of u
    % c = after-spike reset value of v
    % d = after-spike reset of u
    dt = t(2) - t(1); % get time step
    v = zeros(size(t)); % v = membrane potential of neuron
    v(1) = v0; % set initial v
    u = zeros(size(t)); % u = membrane recovery variable
    u(1) = u0; % set initial u
    
    for i = 2:length(t)
        if v(i-1) >= 30
            v(i-1) = 30;
            v(i) = c;
            u(i) = u(i-1) + d;
        else
            dv_dt = 0.04 * v(i-1)^2 + 5 * v(i-1) + 140 - u(i-1) + I(i-1);
            du_dt = a * (b * v(i-1) - u(i-1));
        
            v(i) = v(i-1) + dt*dv_dt;
            u(i) = u(i-1) + dt*du_dt;
        end
    end
end