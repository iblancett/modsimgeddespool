function res = derivcalc(T, S, m, rho, A, c)

numballs = length(S)/4;
speed_threshold = 0.03; % m/s (about 1 in/s)

x = zeros(1,numballs);
y = zeros(1,numballs);
vx = zeros(1,numballs);
vy = zeros(1,numballs);
for i = 1:numballs
    x(i) = S(i*4-3);
    y(i) = S(i*4-2);
    vx(i) = S(i*4-1);
    vy(i) = S(i*4);
    % Stop balls once they reach a certain threshold, otherwise
    % they'll asymptotically approach 0 m/s
    n = norm([vx(i) vy(i)]);
    if (0 < n && n < speed_threshold)
        vx(i) = 0;
        vy(i) = 0;
    end
end

% Force of drag
% Find:
% c = coefficient of drag
% rho = density of ball
% a = cross-section area
% m = mass of ball
fx = - 1/2 * c .* rho .* A .* vx .* (vx.^2 + vy.^2);
fy = - 1/2 * c .* rho .* A .* vy .* (vx.^2 + vy.^2);

% Get acceleration
ax = fx ./ m;
ay = fy ./ m;

T = zeros(length(S),1);

for k = 1:numballs
    T(k*4-3) = vx(k);
    T(k*4-2) = vy(k);
    T(k*4-1) = ax(k);
    T(k*4) = ay(k);
    if (T(end) > 1.7 && k ~= 1 && (max(ax(2:numballs) > 0) || max(ay(2:numballs)) > 0))
        temp = 5;
    end
end

res = T;

end