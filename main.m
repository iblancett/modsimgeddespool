% Plots the path of a pool ball
% By Kyle Combes and Isabel Blancett for Olin College ModSim 2016.

%% Define constants
% given
table_width = 1.17; % m
table_length = 2.34; % m
mass_cue = .17; % kg
mass_regular = .16; % kg
radius_cue = .055; % m
radius_regular = .05; % m

% calculations
volume_cue = 4/3 * pi * radius_cue^3; % m^3
volume_eight = 4/3 * pi * radius_regular^3; % m^3
rho_cue = mass_cue / volume_cue; % kg/m^3
rho_eight = mass_regular / volume_eight; % kg/m^3
a_cue = pi * radius_cue^2; % m^2
a_eight = pi * radius_regular^2; % m^2
c = .015; % unitless - range of 0.005 - 0.015

% Define initial ball positions and velocities
S = [table_width/2, table_length/2-0.5, 0, 1, ...
table_width/2+0.01, table_length/2+0.14, 0, 0, ...
table_width/2-0.05, table_length/2+0.25, 0, 0, ...
table_width/2+0.2, table_length/2, 0, 0];

% create vectors of ball properties for flow function
ball_count = length(S)/4;
m = [mass_cue mass_regular*ones(1,ball_count-1)];
rho = [rho_cue rho_eight*ones(1,ball_count-1)];
A = [a_cue a_eight*ones(1,ball_count-1)];
global radii
radii = [radius_cue radius_regular*ones(1,ball_count-1)];

%% Run the simulation

dRdt = @(ti, Pi) derivcalc(ti, Pi, m, rho, A, c);
options = odeset('Events', @ode_events);
startTime = 0; % seconds
endTime = 40; % seconds
timestep = 0.1; % seconds

T_master = [];
S_master = [];

bounces = 0; % Maximum number of bounces to run for (runaway stop)
global ball_dist_event_dirs

while (bounces < 40 && startTime < endTime-timestep)
    x = zeros(1,ball_count);
    y = zeros(1,ball_count);
    for i = 1:ball_count
        start_index = 4*(i-1)+1;
        x(i) = S(end, start_index);
        y(i) = S(end, start_index+1);
    end

    ball_dist_event_dirs = zeros(1,factorial(ball_count)/factorial(ball_count - 2));
    ball_dists = calc_ball_proximity(x, y, ball_count, radii, ball_dist_event_dirs);
    ball_dist_event_dirs = -ball_dists./abs(ball_dists);

    timespan = startTime:timestep:endTime;
    [t, S] = ode45(dRdt, timespan, S, options);
    T_master = [T_master; t];
    S_master = [S_master; S];
    startTime = t(end) + timestep;
    [S, balls_stopped] = calculate_vectors_after_collision(S(end,:), radii, m);
    if (balls_stopped)
        S_master(end, :) = S(1,:);
        break; % The balls stopped rolling
    end
    bounces = bounces + 1;
end

%% Display results

% Calculate total kinetic energy throughout simulation
rows = size(S_master,1);
K_tot = zeros(rows, 1);
for i = 1:ball_count
    start_index = 4*(i-1)+3;
    vx = S_master(:, start_index);
    vy = S_master(:, start_index+1);
    K_tot = K_tot + 0.5*m(i)*(vx.^2 + vy.^2);
end

table_dims = [0 0 table_width table_length];
axes = [0 table_length 0 table_length];
animate_results(T_master, S_master, axes, radii, table_dims, K_tot);