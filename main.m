% Plots the path of a pool ball

% change to get more balls
numballs = 2;

% given
table_width = 1.17; % m
table_length = 2.34; % m
ball_radius = 0.03; % m
mass_cue = .17; % kg
mass_eight = .16; % kg
radius_cue = .055; % m
global radius_eight
radius_eight = .05; % m

% calculations
volume_cue = 4/3 * pi * radius_cue^3; % m^3
volume_eight = 4/3 * pi * radius_eight^3; % m^3
rho_cue = mass_cue / volume_cue; % kg/m^3
rho_eight = mass_eight / volume_eight; % kg/m^3
a_cue = pi * radius_cue^2; % m^2
a_eight = pi * radius_eight^2; % m^2
c = .015; % unitless - range of 0.005 - 0.015


% create vectors of ball properties for flow function
m = zeros(1, numballs);
m(1) = mass_cue;
rho = zeros(1, numballs);
rho(1) = rho_cue;
a = zeros(1, numballs);
a(1) = a_cue;

for i = 1:numballs
m(i+1) = mass_eight;
rho(i+1) = rho_eight;
a(i+1) = a_eight;
end

dRdt = @(ti, Pi) derivcalc(ti, Pi, m, rho, a, c);
options = odeset('Events', @ode_events);
startTime = 0; % seconds
endTime = 30; % seconds
timestep = 0.1; % seconds

T_master = [];
S_master = [];

% Start ball in the middle with an initial velocity of 1 m/s
S = [1, table_length/2, 0, 1, ...
1.05, table_length/2 + 0.5, 0, 0];

bounces = 0;
ball_count = length(S)/4;
global ball_dist_event_dirs

while (bounces < 20 && startTime < endTime)
    x = zeros(1,ball_count);
    y = zeros(1,ball_count);
    for i = 1:ball_count
        start_index = 4*(i-1)+1;
        x(i) = S(end, start_index);
        y(i) = S(end, start_index+1);
    end

    ball_dist_event_dirs = zeros(1,factorial(ball_count)/(2*factorial(ball_count - 2)));
    ball_dists = calc_ball_proximity(x, y, ball_count, radius_eight, ball_dist_event_dirs);
    ball_dist_event_dirs = -ball_dists./abs(ball_dists);

    timespan = startTime:timestep:endTime;
    [t, S] = ode45(dRdt, timespan, S, options);
    T_master = [T_master; t];
    S_master = [S_master; S];
    startTime = t(end) + timestep;
    S = S(end,:);
    S = calculate_vectors_after_collision(S(end,:), radius_eight);
    if (S == false)
        break; % The balls stopped rolling
    end
    bounces = bounces + 1;
end
% comet(S_master(:,1), S_master(:,2));

table_dims = [0 0 table_width table_length];
axes = [0 table_length 0 table_length];
animate_results(T_master, S_master, axes, radius_cue, table_dims);