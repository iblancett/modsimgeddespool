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
global r
r = zeros(1, numballs);
r(1) = radius_cue;

for i = 1:numballs
m(i+1) = mass_eight;
rho(i+1) = rho_eight;
a(i+1) = a_eight;
r(i+1) = radius_eight;
end

dRdt = @(ti, Pi) derivcalc(ti, Pi, m, rho, a, c);
options = odeset('Events', @ode_events);
timespan = 0:.02:100; % time to run simulation (s)

T_master = [];
S_master = [];

% Start ball in the middle with an initial velocity of 1 m/s
S = [table_width/2, table_length/2, 0.2, 1, ...
table_width/2, table_length/2 + 0.5, 0, 0];
    
figure(1); clf; hold on;
axis([0 table_width 0 table_length ]);
plot(S(5), S(6), 'ko');
bounces = 0;

while (bounces < 1000)
    [t, S] = ode45(dRdt, timespan, S, options);
    timespan = t(end)+.02:.02:100;
    T_master = [T_master; t];
    S_master = [S_master; S];
    S = calculate_vectors_after_collision(S(end,:), r);
    if (S == false)
        break; % The balls stopped rolling
    end
    bounces = bounces + 1;
end

comet(S_master(:,1), S_master(:,2));
%p1 := 