function [resS, balls_stopped] = calculate_vectors_after_collision(S, radii, m, pocket_radius)
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
    table_center_horz = table_width/2;
    table_center_vert = table_length/2;
    
    %% Unpack old values    
    for i = 1:ball_count
        start_index = 4*(i-1)+1;
        x(i) = S(end, start_index);
        y(i) = S(end, start_index+1);
        vx(i) = S(end, start_index+2);
        vy(i) = S(end, start_index+3);
    end
    
    %% Check for collisions with wall and for sinking into a pocket
    for i = 1:ball_count
        radius = radii(i);
        if (x(i) - radius < tolerance || table_width - radius - x(i) < tolerance) % hit left or right wall
            %% Hit left or right wall
            % Check for falling into center pocket along vertical walls
            if (abs(y(i) - table_center_vert) < pocket_radius)
                if (x(i) > table_center_horz)
                    x(i) = table_width; % Fell into right center pocket
                else
                    x(i) = 0; % Fell into left center pocket
                end
                y(i) = table_center_vert;
                vx(i) = 0; vy(i) = 0;
                radii(i) = radii(i)/2; % Just for visual cue
            % Check for falling into a corner pocket
            elseif (abs(y(i) - table_center_vert) < pocket_radius ...
                || abs(y(i) - table_center_vert) > table_center_vert - pocket_radius)
                [x(i), y(i)] = get_sunk_ball_pos(x(i), y(i));
                vx(i) = 0; vy(i) = 0;
                radii(i) = radii(i)/2; % Just for visual cue
            else
                %% Didn't fall in, just bouncing off wall
                vx(i) = -vx(i)/2; % 50% E lost in collision
            end
        elseif (y(i) - radius < tolerance || table_length - radius - y(i) < tolerance)
            %% Hit top or bottom wall
            % Check for falling into a corner pocket
            if (abs(x(i) - table_center_horz) > table_center_horz - pocket_radius)
                [x(i), y(i)] = get_sunk_ball_pos(x(i), y(i));
                vx(i) = 0; vy(i) = 0;
                radii(i) = radii(i)/2; % Just for visual cue
            else
                %% Didn't fall in, just bouncing off wall
                vy(i) = -vy(i)/2; % 50% E lost in collision
            end
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
    resS = size(S);
    for i = 1:ball_count
        start_index = 4*(i-1)+1;
        resS(start_index) = x(i);
        resS(start_index+1) = y(i);
        resS(start_index+2) = vx(i);
        resS(start_index+3) = vy(i);
    end
    
    %% Check if all balls have stopped
    if (max(abs([vx vy])) < tolerance)
        balls_stopped = true;
    else
        balls_stopped = false;
    end
    
    function [x, y] = get_sunk_ball_pos(x, y)
        if (x > table_center_horz)
            x = table_width; % Fell into right center pocket
        else
            x = 0; % Fell into left center pocket
        end
        if (y > table_center_vert)
            y = table_length;
        else
            y = 0;
        end
    end
end