% Plots the path of a pool ball

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
    S = calculate_vectors_after_collision(S);
    if (S == false)
        break; % The balls stopped rolling
    end
    bounces = bounces + 1;
end

comet(S_master(:,1), S_master(:,2));