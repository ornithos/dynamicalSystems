% Spring dynamics (according to discretised ODE)
pend.f = @(x, dt) [x(1,:) + dt*x(2,:); x(2,:) - 9.81 * sin(x(1,:)) * dt];
pend.F = @(x, dt) [1, dt;-9.81*cos(x(1))*dt,  1];

dt     = 0.01;
f      = @(x) pend.f(x, dt);
F      = @(x) pend.F(x, dt);
h      = @(x) sin(x(1,:));
H      = @(x) [cos(x(1)), 0];
qc     = 1;
Q      = [qc * dt^3 / 3, qc * dt^2 / 2; qc * dt^2 / 2, qc * dt];
R      = 0.2;

dsPend  = ds.dynamicalSystem(2, 1, 'x0', [1;-0.6], 0.5, 'evolution', f, F, Q, 'emission', h, H, R, 'data', 400);


%%
dsPend  = dsPend.filter('ekf');
dsPend  = dsPend.smooth('ekf');
dsPend  = dsPend.save(':EKF');

dsPend  = dsPend.filter('ukf');
dsPend  = dsPend.smooth('ukf');
dsPend  = dsPend.save(':UKF');
ds.gui.posteriorGaussGUI(dsPend, ':EKF', ':UKF');

%%
dsPEnd  = dsPend.filterMix('EKF');

%% Sarkka stuff
% ... run pendulum_sim ...
m0      = [1.6;0]; % Slightly off
P0      = 0.1*eye(2);
DT      = 0.01;
g       = 9.81;
Q       = 0.01*[DT^3/3 DT^2/2; DT^2/2 DT];
R       = 0.1;

% f       = [m(1)+m(2)*DT; m(2)-g*sin(m(1))*DT];
% f       = @(x) [x(1,:) + x(2,:)*dt; x(2,:) - g*sin(x(1,:))*dt];
% F       = @(x) [1, dt; -g*cos(x(1))*dt, 1];
% F       = [1 DT; -g*cos(m(1))*DT 1];
dsPend  = ds.dynamicalSystem(2, 1, 'x0', m0, P0, 'evolution', f, F, Q, 'emission', h, H, R, 'data', Y, 'xtrue', X);


%%

u_inp  = repmat(0.05*[0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0.1 2], [2,20]);
B      = eye(2);
C      = [1,0];
pend.f = @(x, u, dt) [x(1,:) + dt*x(2,:); x(2,:) - 9.81 * sin(x(1,:)) * dt] + B*u;
pend.F = @(x, u, dt) [1, dt;-9.81*cos(x(1))*dt,  1];

dt     = 0.01;
f      = @(x,u) pend.f(x, u, dt);
F      = @(x,u) pend.F(x, u, dt);
h      = @(x,u) sin(x(1,:)) + C*u;
H      = @(x,u) [cos(x(1)), 0];
qc     = 1;
Q      = [qc * dt^3 / 3, qc * dt^2 / 2; qc * dt^2 / 2, qc * dt];
R      = 0.2;

dsPend  = ds.dynamicalSystem(2, 1, 'x0', [1;-0.6], 0.5, 'evolution', f, F, Q, 'emission', h, H, R, 'data', 400, 'control', u_inp, true, true);

%%

dsPend  = dsPend.filter('ekf');
dsPend  = dsPend.smooth('ekf');
dsPend  = dsPend.save(':EKF');

dsPend  = dsPend.filter('ukf');
dsPend  = dsPend.smooth('ukf');
dsPend  = dsPend.save(':UKF');
ds.gui.posteriorGaussGUI(dsPend, ':EKF', ':UKF');