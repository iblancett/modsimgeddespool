function [value, isterminal, direction] = ode_events(~, S)
    
    table_width = 1.17; % m
    table_length = 2.34; % m
    ball_radius = 0.03; % m

    x = S(1);
    y = S(2);
    vx = S(3);
    vy = S(4);
    
    dist_to_left_edge = x - ball_radius;
    dist_to_right_edge = table_width - x - ball_radius;
    dist_to_top_edge = table_length - y - ball_radius;
    dist_to_bottom_edge = y - ball_radius;
    
    value = [dist_to_left_edge, dist_to_right_edge, ...
        dist_to_top_edge, dist_to_bottom_edge, ...
        vx, vy];
    isterminal = [1 1 1 1 1 1];
    direction = [0 0 0 0 0 0];
    
end