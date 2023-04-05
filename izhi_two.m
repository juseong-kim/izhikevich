function [v1, v2, I1, I2, u1, u2, g1, g2] = izhi_two(t,v0,u0,i0,g0,a,b,c,d,tau,gbar,vr, ei)
    % use default values if not provided
    % a = scale of u
    % b = sensitivity of u
    % c = after-spike reset value of v
    % d = after-spike reset of u
    % tau = time constant
    % gbar = maximal conductance
    % ddelta = array of size(t), 1 at time of pre-synaptic AP, otherwise 0
    
    dt = t(2) - t(1); % get time step
    v1 = zeros(size(t)); % membrane potential of neuron 1
    v1(1) = v0; % set initial v1
    v2 = zeros(size(t)); % membrane potential of neuron 2
    v2(1) = v0; % set initial v2
    u1 = zeros(size(t)); % u = membrane recovery variable
    u1(1) = u0; % set initial u
    u2 = zeros(size(t)); % u = membrane recovery variable
    u2(1) = u0; % set initial u
%     I1 = zeros(size(t));
%     I1(1) = i0;
    I1 = [zeros(1,20), 10*ones(1,size(t,2)-20)];   
    if ei
        I2 = zeros(size(t));
        I2(1) = i0;
    else
%         I2 = [zeros(1,20), 10*ones(1,size(t,2)-20)];
        I2 = zeros(size(t));
        I2(1) = i0;
    end
    g1 = zeros(size(t));
    g1(1) = g0;
    g2 = zeros(size(t));
    g2(1) = g0;
    ddelta1 = zeros(size(t));
    ddelta2 = zeros(size(t));
  
    for i = 2:length(t)
        dg_dt1 = -1 * g1(i-1) / tau + gbar * ddelta1(i-1);
        g1(i) = g1(i-1) + dt * dg_dt1;
        dg_dt2 = -1 * g2(i-1) / tau + gbar * ddelta2(i-1);
        g2(i) = g2(i-1) + dt * dg_dt2;

        dv_dt1 = 0.04 * v1(i-1)^2 + 5 * v1(i-1) + 140 - u1(i-1) + I1(i-1);
        du_dt1 = a * (b * v1(i-1) - u1(i-1));

        v1(i) = v1(i-1) + dt*dv_dt1;
        u1(i) = u1(i-1) + dt*du_dt1;
                
        dv_dt2 = 0.04 * v2(i-1)^2 + 5 * v2(i-1) + 140 - u2(i-1) + I2(i-1);
        du_dt2 = a * (b * v2(i-1) - u2(i-1));

        v2(i) = v2(i-1) + dt*dv_dt2;
        u2(i) = u2(i-1) + dt*du_dt2;
        if ei
            I2(i) = -1 * g2(i) * (v2(i) - vr);
        else 
            I2(i) = -1 * g2(i) * (v2(i) - vr);
        end

        if v1(i-1) >= 30
            v1(i-1) = 30;
            v1(i) = c;
            u1(i) = u1(i-1) + d;
            ddelta2(i) = 1;
        end
        if v2(i-1) >=30
            v2(i-1) = 30;
            v2(i) = c;
            u2(i) = u2(i-1) + d;
        end
    end
end