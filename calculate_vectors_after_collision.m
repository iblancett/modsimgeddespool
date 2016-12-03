function res = calculate_vectors_after_collision(S, ball_radius)
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
    if (max(abs([vx vy])) < tolerance)
        res = false;
        return;
    end
    
    %% Check for collisions with wall
    for i = 1:ball_count
        if (x(i) - ball_radius < tolerance || table_width - ball_radius - x(i) < tolerance) % hit left or right wall
            vx(i) = -vx(i);
        elseif (y(i) - ball_radius < tolerance || table_length - ball_radius - y(i) < tolerance) % hit top or bottom wall
            vy(i) = -vy(i); 
        end
    end
    
    %% Check for collisions with other balls
    P = [x' y'];
    V = [vx' vy']; % row for each ball
    for i = 1:ball_count
        for j = 2:ball_count
            if i<j && abs(x(i)-x(j)) < 2 * ball_radius + tolerance &&...
                    abs(y(i)-y(j)) < 2 * ball_radius + tolerance
                K_i = 0.5*m*((vx(1)+vy(1)).^2 + vx(2)+vy(2).^2);
                P_i = m*(V(i,:)+V(2,:));
                N = P(i,:) - P(j,:);
                magdP = norm(N);
                N = N./magdP;
                % Calc length of each component of the velocity
                % vectors along N
                a1 = dot(V(i,:),N);
                a2 = dot(V(j,:),N);
                optimizedP = (a1 - a2) / m;
                V(i,:) = V(i,:) - optimizedP .* N * m;
                V(j,:) = V(j,:) + optimizedP .* N * m;
                vx = V(:,1)';
                vy = V(:,2)';
                K_f = 0.5*m*((vx(1)+vy(1)).^2 + vx(2)+vy(2).^2);
                P_f = m*(V(i,:)+V(2,:));
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