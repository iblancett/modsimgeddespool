function res = calculate_vectors_after_collision(S, radii, m)
    % Calculates the new position and velocity vectors after a collision
    % or the balls have stopped rolling. Returns an updated version of the
    % same vector as was passed in or 'false' if all the balls have stopped
    % rolling.
    
    ball_count = length(S)/4;
    table_width = 1.17; % m
    table_length = 2.34; % m
    tolerance = 1e-4; % Tolerance when checking for collisions
    
    x = zeros(1,ball_count);
    y = zeros(1,ball_count);
    vx = zeros(1,ball_count);
    vy = zeros(1,ball_count);
    
    %% Unpack old values    
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
        radius = radii(i);
        if (x(i) - radius < tolerance || table_width - radius - x(i) < tolerance) % hit left or right wall
            vx(i) = -vx(i);
        elseif (y(i) - radius < tolerance || table_length - radius - y(i) < tolerance) % hit top or bottom wall
            vy(i) = -vy(i); 
        end
    end
    
    %% Check for collisions with other balls
    P = [x' y'];
    V = [vx' vy']; % row for each ball
    for i = 1:ball_count
        for j = i+1:ball_count
            dist_btwn_centers = radii(i) + radii(j);
            if abs(x(i)-x(j)) < dist_btwn_centers + tolerance &&...
                    abs(y(i)-y(j)) < dist_btwn_centers + tolerance
                K_i = 0.5*m*(V(:,1).^2+V(:,2).^2);
                P_i = m*(V(:,1).^2+V(:,2).^2);
                N = P(i,:) - P(j,:);
                magdP = norm(N);
                N = N./magdP;
                % Calc length of each component of the velocity
                % vectors along N
                a1 = dot(V(i,:),N);
                a2 = dot(V(j,:),N);
                optimizedP = (a1 - a2) * length(m) / sum(m);
                V(i,:) = V(i,:) - optimizedP .* N * m(i);
                V(j,:) = V(j,:) + optimizedP .* N * m(j);
                vx = V(:,1)';
                vy = V(:,2)';
                K_f = 0.5*m*(V(:,1).^2+V(:,2).^2);
                P_f = m*(V(:,1).^2+V(:,2).^2);
                if (K_f > K_i)
                    pause(0.001); % Just a place to put a breakpoint for debugging
                end
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