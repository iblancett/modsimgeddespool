function res = animate_results(Time, S, axes, radii, dims, K_tot)
    
    % Unpack positions
    ball_count = size(S,2)/4;
    X = zeros(size(S,1), ball_count);
    Y = zeros(size(S,1), ball_count);
    for t = 1:ball_count
        X(:,t) = S(:,4*t-3);
        Y(:,t) = S(:,4*t-2);
    end
    maxE = max(K_tot)*1.1;
    dt = Time(2) - Time(1);
    
    for t = 1:length(Time)
        figure(1);
        clf; hold on;
        axis(axes);
        axis square
        rectangle('Position', dims);
        for b = 1:ball_count
            draw_circle(X(t,b), Y(t,b), radii(b));
        end
        figure(2); clf; hold on;
        
        plot(Time, K_tot);
        %plot(Time, P_tot);
        plot([Time(t) Time(t)], [0 maxE]);
        %legend('Kinetic energy', 'Momentum');
        pause(dt);
    end
    
    function draw_circle(x,y,r)
        %x and y are the coordinates of the center of the circle
        %r is the radius of the circle
        %0.01 is the angle step, bigger values will draw the circle faster but
        %you might notice imperfections (not very smooth)
        ang=0:0.05:2*pi; 
        xp=r*cos(ang);
        yp=r*sin(ang);
        plot(x+xp,y+yp);
    end
end