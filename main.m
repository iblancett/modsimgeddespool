% Plots the path of a pool 

table_width = 1.17; % m
table_length = 2.34; % m
ball_radius = 0.03; % m

dRdt = @(ti, Pi) derivcalc(ti, Pi);
options = odeset('Events', @ode_events);
timespan = [0 2000]; % time to run simulation (s)

T_master = [];
S_master = [];

% Start ball in the middle with an initial velocity of 1 m/s
S = [table_width/2, table_length/2, 1, 1];

bounces = 0;

while (bounces < 10)
    [t, S] = ode45(dRdt, timespan, S, options);
    T_master = [T_master; t];
    S_master = [S_master; S];
    x = S(end,1);
    y = S(end,2);
    vx = S(end,3);
    vy = S(end,4);
    
    tolerance = 1e-4;
    
    S = S(end,:);
    
    if (vx == 0 && vy == 0) % ball stopped, stop simulation
        break;
    end
    if (x - ball_radius < tolerance || table_width - ball_radius - x < tolerance) % hit left or right wall
        S(3) = -vx;
        bounces = bounces + 1;
    elseif (y - ball_radius < tolerance || table_length - ball_radius - y < tolerance) % hit top or bottom wall
        S(4) = -vy;
        bounces = bounces + 1;
    end
    
end

comet(S_master(:,1), S_master(:,2));