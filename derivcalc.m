function res = derivcalc(~, S, m, rho, A, c)

numballs = length(S)/4;

x = zeros(1,numballs);
y = zeros(1,numballs);
vx = zeros(1,numballs);
vy = zeros(1,numballs);
for i = 1:numballs
    x(i) = S(i*4-3);
    y(i) = S(i*4-2);
    vx(i) = S(i*4-1);
    vy(i) = S(i*4);
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
end

res = T;

end