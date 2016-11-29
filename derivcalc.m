function res = derivcalc(~, S, m, rho, a, c)

numballs = length(S)/4;

for i = 1:numballs
    x(i) = S(i*4-3);
    y(i) = S(i*4-2);
    vx(i) = S(i*4-1);
    vy(i) = S(i*4);
end

for j = 1:numballs
% Force of drag
% Find:
% c = coefficient of drag
% rho = density of ball
% a = cross-section area
% m = mass of ball

fx(j) = - 1/2 .* c .* rho(j) .* a(j) .* vx(j) .* sqrt(vx(j).^2 + vy(j).^2);
fy(j) = - 1/2 .* c .* rho(j) .* a(j) .* vy(j) .* sqrt(vx(j).^2 + vy(j).^2);

% Get acceleration
ax(j) = fx(j) ./ m(j);
ay(j) = fy(j) ./ m(j);
end

T = zeros(length(S),1);

for k = 1:numballs
    T(k*4-3) = vx(k);
    T(k*4-2) = vy(k);
    T(k*4-1) = ax(k);
    T(k*4) = ay(k);
end

res = T;

end