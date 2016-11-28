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
        x(i) = S(end, start_index);
        y(i) = S(end, start_index+1);
        vx(i) = S(end, start_index+2);
        vy(i) = S(end, start_index+3);
    end
    
    x = S(1);
    y = S(2);
    vx = S(3);
    vy = S(4);
    
    % Calculate proximity to edges
    dist_le = x - ball_radius;
    dist_re = table_width - x - ball_radius;
    dist_te = table_length - y - ball_radius;
    dist_be = y - ball_radius;
    
    value = [dist_le dist_re dist_te dist_be vx vy];
    isterminal = ones(size(value));
    direction = -1*ones(size(value));
    
end