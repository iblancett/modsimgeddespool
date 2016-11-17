% Plots the path of a pool 

table_width = 1.17; % m
table_length = 2.34; % m
ball_radius = 0.03; % m

dPdt = @(ti, Pi) deriv_calc(ti, Pi);
options = odeset('Events', @ode_events);

% Start ball in the middle with an initial velocity of 1 m/s
S = [table_width/2, table_length/2, 1, 1];
