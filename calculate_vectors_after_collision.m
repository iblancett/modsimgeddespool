function res = calculate_vectors_after_collision(S, ball_radius)
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
    
    %% Organize old values    
    for i = 1:ball_count
        start_index = 4*(i-1)+1;
        x(i) = S(end, start_index);
        y(i) = S(end, start_index+1);
        vx(i) = S(end, start_index+2);
        vy(i) = S(end, start_index+3);
    end
    
    %% Check if all balls have stopped
    for i = 1:ball_count
        if (abs(vx(i)) > tolerance && abs(vy(i)) > tolerance)
            break; % At least one ball is still moving
        end
        res = false;
        return
    end
    
    %% Check for collisions with wall
    for i = 1:ball_count
        if (x(i) - ball_radius < tolerance || table_width - ball_radius - x(i) < tolerance) % hit left or right wall
            vx(i) = -vx(i);
        elseif (y(i) - ball_radius < tolerance || table_length - ball_radius - y(i) < tolerance) % hit top or bottom wall
            vy(i) = -vy(i);
        end
    end
    %{
    %% Check for collisions with other balls
    for i = 1:ball_count
        
        for n = 1:ball_count
            if (n == i)
                continue; % Don't check for collision with itself
            end
            
            xi = x(i);
            yi = y(i);
            xn = x(n);
            yn = y(n);
            
            % Calculate distance between balls
            d = sqrt( (xi - xn)^2 + (yi - yn)^2 );
            
            if ((d - 2*ball_radius) < tolerance)
                theta = atand( (yi-yn)/(xi-xn) );
            end
            
        end
    end
    %}
    
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