function res = animate_results(Time, S, axes, radii, dims, K_tot, ...
    speeds, pocket_radius, init_theta, init_speed)
    
    felt_color = [10, 108, 3]./255;
    ball_colors = {'w', 'k', [224 102 102]/255, [246 178 107]/255, [103 78 167]/255};

    % Unpack positions
    ball_count = size(S,2)/4;
    X = zeros(size(S,1), ball_count);
    Y = zeros(size(S,1), ball_count);
    for t = 1:ball_count
        X(:,t) = S(:,4*t-3);
        Y(:,t) = S(:,4*t-2);
    end
    maxE = max(K_tot);
    maxV = max(max(speeds));
    dt = Time(2) - Time(1);
    
    % Calc values needed for plotting pockets
    piover2 = pi/2;
    pi3over2 = 3*pi/2;
    pi5over2 = 5*pi/2;
    table_width = dims(3);
    table_length = dims(4);
    pocket_xs = [0 0 0 table_width table_width table_width];
    pocket_ys = [table_length table_length/2 0 0 table_length/2 table_length];
    start_thetas = [0 piover2 piover2 pi pi3over2 pi3over2];
    end_thetas = [pi3over2 pi3over2 2*pi pi5over2 pi5over2 3*pi];
    
    figure(1);
    for t = 1:length(Time)
        clf;
        p = subplot(1,2,1); hold on;
        p.Color = felt_color;
        init_theta_deg = ceil(init_theta*180/pi);
        title(strcat('Init theta: ', num2str(init_theta_deg), '°, init speed: ', ...
            num2str(init_speed), 'm/s'));
        axis(axes);
        axis square
        
        % Draw table
        rectangle('Position', dims);
        
        % Draw pockets
        for i = 1:6
            draw_arc(pocket_xs(i), pocket_ys(i), pocket_radius, ...
                start_thetas(i), end_thetas(i), 'k');
        end
        
        % Draw balls
        for b = 1:ball_count
            x = X(t,b);
            y = Y(t,b);
            radius = radii(b);
            if (x == dims(1) || x == dims(3) || y == dims(2) || y == dims(4))
                radius = radius/2;
            end
            circles(x, y, radius, 'facecolor', ball_colors{b});
        end
        
        % Plot kinetic energy
        subplot(2,2,2); hold on;
        
        plot(Time, K_tot);
        %plot(Time, P_tot);
        plot([Time(t) Time(t)], [0 maxE]);
        title('Total Kinetic Energy Over Time');
        xlabel('Time (s)');
        ylabel('Energy (J)');
        %legend('Kinetic energy', 'Momentum');
        
        % Plot ball speeds
        subplot(2,2,4); hold on;
        title('Ball Velocities Over Time');
        xlabel('Time(s)');
        ylabel('Speed (m/s)');
        plot(Time, speeds);
        plot([Time(t) Time(t)], [0 maxV]);
        pause(dt/5);
    end
    
    function draw_circle(x,y,r,color)
        %x and y are the coordinates of the center of the circle
        %r is the radius of the circle
        %0.01 is the angle step, bigger values will draw the circle faster but
        %you might notice imperfections (not very smooth)
        ang=0:0.05:2*pi; 
        xp=r*cos(ang);
        yp=r*sin(ang);
        plot(x+xp,y+yp,color);
    end

    function draw_arc(x,y,r,startTheta,endTheta,color)
        %x and y are the coordinates of the center of the circle
        %r is the radius of the circle
        %0.01 is the angle step, bigger values will draw the circle faster but
        %you might notice imperfections (not very smooth)
        ang=startTheta:0.05:endTheta; 
        xp=r*cos(ang);
        yp=r*sin(ang);
        plot(x+xp,y+yp,color);
    end
end