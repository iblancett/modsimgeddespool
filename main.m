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

for i = numballs
    m(i+1) = mass_eight;
    rho(i+1) = rho_eight;
    a(i+1) = a_eight;
end

dRdt = @(ti, Pi) derivcalc(ti, Pi, m, rho, a, c);
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
    
    if (vx < tolerance && vy < tolerance) % ball stopped, stop simulation
        break;
    elseif (x - ball_radius < tolerance || table_width - ball_radius - x < tolerance) % hit left or right wall
        S(3) = -vx;
        bounces = bounces + 1;
    elseif (y - ball_radius < tolerance || table_length - ball_radius - y < tolerance) % hit top or bottom wall
        S(4) = -vy;
        bounces = bounces + 1;
    end
    
end

comet(S_master(:,1), S_master(:,2));