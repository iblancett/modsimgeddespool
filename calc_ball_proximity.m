function ball_dists = calc_ball_proximity(x, y, ball_count, ball_radius, ball_dist_dirs)
    % Calculate proximity of balls to each other
    num_combinations = factorial(ball_count)/(2*factorial(ball_count - 2));
    diam = 2 * ball_radius;
    ball_dists = zeros(1, num_combinations*2);
    for i = 1:ball_count-1
        for j = i+1:ball_count
            index = (i-1)*(ball_count-1)+1;
            ball_dists(2*index-1) = x(i)-x(j) + diam * ball_dist_dirs(index);
            ball_dists(2*index) = y(i)-y(j) + diam * ball_dist_dirs(index);
        end
    end
end