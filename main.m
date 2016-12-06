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
pocket_radius = 0.15; % m

% calculations
volume_cue = 4/3 * pi * radius_cue^3; % m^3
volume_eight = 4/3 * pi * radius_regular^3; % m^3
rho_cue = mass_cue / volume_cue; % kg/m^3
rho_eight = mass_regular / volume_eight; % kg/m^3
a_cue = pi * radius_cue^2; % m^2
a_eight = pi * radius_regular^2; % m^2
c = 0.35; % unitless - range of 0.005 - 0.015

% Define initial ball positions. Cue velocity set later.
S_init = [table_width/4-radius_regular, table_length/3,                 0, 0, ... % cue
    table_width/4-radius_regular,   table_length/2-3*radius_regular,    0, 0, ... % 8 ball
    table_width/8,                  table_length/2-4*radius_regular,    0, 0, ... % red
    table_width/4-3*radius_regular, table_length/2+3*radius_regular,  0, 0, ... % orange
    table_width/4-3*radius_regular+sqrt(2)*radius_regular, table_length/2+3*radius_regular+sqrt(2)*radius_regular,      0, 0]; % purple

% create vectors of ball properties for flow function
ball_count = length(S_init)/4;
m = [mass_cue mass_regular*ones(1,ball_count-1)];
rho = [rho_cue rho_eight*ones(1,ball_count-1)];
A = [a_cue a_eight*ones(1,ball_count-1)];
global radii
radii = [radius_cue radius_regular*ones(1,ball_count-1)];

%% Simulation parameters

theta_step = pi/24;
max_speed = 2; % maximum speed John Geddes can get the cue ball moving
sim_init_thetas = 1:theta_step:2*pi;
sim_init_speed = 2;%:1:max_speed;

%% Run the simulation

dRdt = @(ti, Pi) derivcalc(ti, Pi, m, rho, A, c);
options = odeset('Events', @ode_events);
start_time = 0; % seconds
end_time = 15; % seconds
timestep = 0.1; % seconds


global ball_dist_event_dirs

for t = 1:length(sim_init_thetas)
    for s = 1:length(sim_init_speed)

        S = S_init;
        T_master = [];
        S_master = [];
        
        % Set cue initial velocity
        S(3) = sim_init_speed(s)*cos(sim_init_thetas(t));
        S(4) = sim_init_speed(s)*sin(sim_init_thetas(t));
        
        result = 0;
        bounces = 0; % Maximum number of bounces to run for (runaway stop)
        resume_time = start_time;

        while (resume_time < end_time-timestep)
            x = zeros(1,ball_count);
            y = zeros(1,ball_count);
            for i = 1:ball_count
                start_index = 4*(i-1)+1;
                x(i) = S(end, start_index);
                y(i) = S(end, start_index+1);
            end

            ball_dist_event_dirs = zeros(ball_count, ball_count, 2);
            ball_dists = calc_ball_proximity(x, y, ball_count, radii, ball_dist_event_dirs);
            ball_dist_event_dirs = -ball_dists./abs(ball_dists);

            timespan = resume_time:timestep:end_time;
            [T, S] = ode45(dRdt, timespan, S, options);
            T_master = [T_master; T];
            S_master = [S_master; S];
            resume_time = T(end) + timestep;
            [S, balls_stopped] = calculate_vectors_after_collision(S(end,:), radii, m, pocket_radius);

            % Check if simulation has ended
            x = zeros(1,ball_count);
            y = zeros(1,ball_count);
            vx = zeros(1,ball_count);
            vy = zeros(1,ball_count);
            for i = 1:ball_count
                start_index = 4*(i-1)+1;
                x(i) = S(end, start_index);
                y(i) = S(end, start_index+1);
                vx(i) = S(end, start_index+2);
                vy(i) = S(end, start_index+3);
            end

            % Check which balls are pocketed
            % Cue first
            x_dist_center = abs(x-table_width/2);
            y_dist_center = abs(y-table_length/2);
            if (x_dist_center(1) == table_width/2 ...
                    && y_dist_center(1) == table_length/2) % Scratched
                result = 0;
                break;
            end
            % Now check eight ball
            if (x_dist_center(2) == table_width/2 ...
                    && y_dist_center(2) == table_length/2)
                % Eight ball pocketed, check if other balls are as well
                if (x_dist_center(3:ball_count) == table_width/2*ones(1,ball_count-2) ...
                        && y_dist_center(3:ball_count) == table_lenth/2*ones(1,ball_count-2))
                    % All other balls (aside from cue) pocketed already
                    result = 1;
                    break;
                else
                    % Not all other balls pocketed yet, 8 ball sunk too early
                    break;
                end
            end
            % Cue nor eight ball pocketed, check if balls are stopped
            if (balls_stopped)
                S_master(end, :) = S(1,:);
                break; % The balls stopped rolling
            end
            % Continue simulation
            bounces = bounces + 1;
        end
        sim(t,s).init_theta = sim_init_thetas(t);
        sim(t,s).init_speed = sim_init_speed(s);
        sim(t,s).result = result;
        sim(t,s).T = T_master;
        sim(t,s).S = S_master;
        sim(t,s).bounces = bounces;

    %% Display results

    % Calculate total kinetic energy throughout simulation
    rows = size(S_master,1);
    K_tot = zeros(rows, 1);
    speed = zeros(size(S_master,1), ball_count);
    for i = 1:ball_count
        start_index = 4*(i-1)+3;
        vx = S_master(:, start_index);
        vy = S_master(:, start_index+1);
        K_tot = K_tot + 0.5*m(i)*(vx.^2 + vy.^2);
        speed(:,i) = sqrt(vx.^2 + vy.^2);
    end

    table_dims = [0 0 table_width table_length];
    axes = [-pocket_radius table_length+pocket_radius -pocket_radius ...
        table_length+pocket_radius];
    figure(1); clf;
    animate_results(T_master, S_master, axes, radii, table_dims, K_tot, ...
        speed, pocket_radius, sim_init_thetas(t), sim_init_speed(s));
    end
end