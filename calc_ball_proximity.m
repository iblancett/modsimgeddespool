function ball_dists = calc_ball_proximity(x, y, ball_count, radii, ball_dist_dirs)
    % Calculate proximity of balls to each other
%     num_combinations = factorial(ball_count)/(2*factorial(ball_count - 2));
    ball_dists = zeros(size(ball_dist_dirs));
    for i = 1:ball_count-1
        for j = i+1:ball_count
            dist_btwn_centers = radii(i) + radii(j);
            ball_dists(i,j,1) = x(i)-x(j) + dist_btwn_centers * ball_dist_dirs(i,j,1);
            ball_dists(i,j,2) = y(i)-y(j) + dist_btwn_centers * ball_dist_dirs(i,j,2);
        end
    end
end