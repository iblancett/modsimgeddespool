function [value, isterminal, direction] = ode_events(~, S)
    
    table_width = 1.17; % m
    table_length = 2.34; % m
    ball_radius = 0.03; % m
    ball_count = length(S)/4;
    
    x = zeros(1,ball_count);
    y = zeros(1,ball_count);
    vx = zeros(1,ball_count);
    vy = zeros(1,ball_count);
    
    %% Organize old values    
    for i = 1:ball_count
        start_index = 4*(i-1)+1;
        x(i) = S(start_index, end);
        y(i) = S(start_index+1, end);
        vx(i) = S(start_index+2, end);
        vy(i) = S(start_index+3, end);
    end
    
    % Calculate proximity to edges
    dist_le = x - ball_radius;
    dist_re = table_width - x - ball_radius;
    dist_te = table_length - y - ball_radius;
    dist_be = y - ball_radius;
    
    % Calculate proximity of balls to each other
    num_combinations = factorial(ball_count)/(2*factorial(ball_count - 2));
    ball_dists = zeros(1, num_combinations*2);
    for i = 1:ball_count-1
        for j = i+1:ball_count
            index = (i-1)*(ball_count-1)+1;
            ball_dists(2*index-1) = x(i)-x(j) - 2 * ball_radius;
            ball_dists(2*index) = y(i)-y(j) - 2 * ball_radius;
        end
    end
    
    value = [dist_le dist_re dist_te dist_be vx vy ball_dists];
    isterminal = ones(size(value));
    direction = [-1*ones(1,length(value)-2*ball_count) zeros(1,2*ball_count)];
    
end