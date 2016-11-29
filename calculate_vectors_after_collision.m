function res = calculate_vectors_after_collision(S, r)
    % Calculates the new position and velocity vectors after a collision
    % or the balls have stopped rolling. Returns an updated version of the
    % same vector as was passed in or 'false' if all the balls have stopped
    % rolling.
    
    ball_count = length(S)/4;
    table_width = 1.17; % m
    table_length = 2.34; % m
    tolerance = 1e-4; % Tolerance when checking for collisions
    m = .165; % kg
    
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
    
    %% Check if all balls have stopped
    if (max([vx vy]) < tolerance)
        res = false;
        return;
    end
    
    %% Check for collisions with wall
    for i = 1:ball_count
        if (x(i) - r(i) < tolerance || table_width - r(i) - x(i) < tolerance) % hit left or right wall
            vx(i) = -vx(i);
        elseif (y(i) - r(i) < tolerance || table_length - r(i) - y(i) < tolerance) % hit top or bottom wall
            vy(i) = -vy(i); 
        end
    end
    
    %% Check for collisions with other balls
    P = [x' y'];
    V = [vx' vy']; % row for each ball
    for i = 2:ball_count
        for j = 1:ball_count
            if i>j && abs(x(i)-x(j)) < r(i) + r(j) + tolerance &&...
                    abs(y(i)-y(j)) < r(i) + r(j) + tolerance
                
                % Compares each pair of balls only once & checks that balls
                % are not being compared to themselves
                % Checks if current pair of balls are in proximity of each
                % other to be colliding
             
                N = P(i,:) - P(j,:); % subtract coordinates of balls
                N = norm(N); % normalizes N
                a1 = dot(V(i,:),N); % finds the velocity vector toward the other ball
                a2 = dot(V(j,:),N); 
                V(i,:) = V(i,:) - (2 .*(a1 - a2) ./ m .* N); % computes new velocity
                V(j,:) = V(j,:) - (2 .*(a1 - a2) ./ m .* N);
                vx = V(:,1)'; % transposes input into vector
                vy = V(:,2)';
            end
        end    
    end     
  
    
    
    %% Pack result
    res = size(S);
    for i = 1:ball_count
        start_index = 4*(i-1)+1;
        res(start_index) = x(i);
        res(start_index+1) = y(i);
        res(start_index+2) = vx(i);
        res(start_index+3) = vy(i);
    end
    
end