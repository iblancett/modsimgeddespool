function [value, isterminal, direction] = ode_events(~, S)
    
    table_width = 1.17; % m
    table_length = 2.34; % m
    speed_threshold = 0.1; % m/s (about 1 in/s)
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
    
    ball_dists_1d = []; %zeros(size(ball_dists,1) * size (ball_dists,2) * 2);
    ball_dists_event_dirs_1d = [];
    for i = 1:ball_count-1
        for j = 2:ball_count
            ball_dists_1d = [ball_dists_1d ball_dists(i,j,1) ball_dists(i,j,2)];
            ball_dists_event_dirs_1d = [ball_dists_event_dirs_1d ball_dist_event_dirs(i,j,1) ball_dist_event_dirs(i,j,2)];
        end
    end
    
    ball_combins = length(ball_dists_1d);
    
%     v_mag = sqrt(vx.^2 + vy.^2);
%     for i = 1:ball_count
%         if (v_mag(i) < speed_threshold && v_mag(i) > 0)
%             vx(i) = 0;
%             vy(i) = 0;
%         end
%     end
%     if (min(v_mag(v_mag>0)) < 0.2)
%         vx_thresh = vx ./ v_mag * speed_threshold;
%         vy_thresh = vy ./ v_mag * speed_threshold;
%         absvx = abs(vx);
%         absvy = abs(vy);
%         t1 = abs(absvx - abs(vx_thresh));
%         t2 = abs(absvy - abs(vy_thresh));
%         reduced_vx = t1.*vx./absvx;
%         reduced_vy = t2.*vy./absvy;
%         % Set NaNs back to 0 (originally 0)
%         reduced_vx(isnan(reduced_vx)) = 0;
%         reduced_vy(isnan(reduced_vy)) = 0;
%         newmag = norm([reduced_vx' reduced_vy']);
%     else
%         reduced_vx = vx;
%         reduced_vy = vy;
%     end
    value = [dist_le dist_re dist_te dist_be vx vy ball_dists_1d];
    isterminal = ones(size(value));
    direction = [-1*ones(1,length(value)-ball_combins) ball_dists_event_dirs_1d];
    
end