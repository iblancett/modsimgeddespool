function [value, isterminal, direction] = ode_events(~, S)
    
    table_width = 1.17; % m
    table_length = 2.34; % m
    global radii ball_dist_event_dirs
    ball_count = length(S)/4;
    
    x = zeros(1,ball_count);
    y = zeros(1,ball_count);
    vx = zeros(1,ball_count);
    vy = zeros(1,ball_count);
    
    % Unpack old values
    for i = 1:ball_count
        start_index = 4*(i-1)+1;
        x(i) = S(start_index, end);
        y(i) = S(start_index+1, end);
        vx(i) = S(start_index+2, end);
        vy(i) = S(start_index+3, end);
    end
    
    ball_dists = calc_ball_proximity(x, y, ball_count, radii, ball_dist_event_dirs);
    
    % Calculate proximity to edges
    dist_le = x - radii;
    dist_re = table_width - x - radii;
    dist_te = table_length - y - radii;
    dist_be = y - radii;
    
    ball_combins = length(ball_dists);
    
    value = [dist_le dist_re dist_te dist_be vx vy ball_dists];
    isterminal = ones(size(value));
    direction = [-1*ones(1,length(value)-ball_combins) ball_dist_event_dirs];
    
end