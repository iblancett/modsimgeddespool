function ball_dists = calc_ball_proximity(x, y, ball_count, radii, ball_dist_dirs)
    % Calculate proximity of balls to each other
    num_combinations = factorial(ball_count)/(2*factorial(ball_count - 2));
    ball_dists = zeros(1, num_combinations*2);
    for i = 1:ball_count-1
        for j = i+1:ball_count
            dist_btwn_centers = radii(i) + radii(j);
            index = (i-1)*(ball_count-1)+1;
            ball_dists(2*index-1) = x(i)-x(j) + dist_btwn_centers * ball_dist_dirs(index);
            ball_dists(2*index) = y(i)-y(j) + dist_btwn_centers * ball_dist_dirs(index+1);
        end
    end
end