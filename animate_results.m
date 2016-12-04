function res = animate_results(Time, S, axes, radius, dims)
    
    % Unpack positions
    ball_count = size(S,2)/4;
    X = zeros(size(S,1), ball_count);
    Y = zeros(size(S,1), ball_count);
    for i = 1:ball_count
        X(:,i) = S(:,4*i-3);
        Y(:,i) = S(:,4*i-2);
    end
    
%     ball_colors = {'ko', 'go'};
    
    dt = Time(2) - Time(1);
    for i = 1:length(Time)
        clf; hold on;
        axis(axes);
        axis square
        rectangle('Position', dims);
%         for b = 1:ball_count
%             c = ball_colors(b);
            draw_circle(X(i,1), Y(i,1), radius);
            draw_circle(X(i,2), Y(i,2), radius);
%         end
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