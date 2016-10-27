% Spring dynamics (according to discretised ODE)
pend.f = @(x, dt) [x(1) + dt*x(2); x(2) - 9.81 * sin(x(1)) * dt];
pend.F = @(x, dt) [1, dt;-9.81*cos(x(1))*dt,  1];

dt     = 0.01;
f      = @(x) pend.f(x, dt);
F      = @(x) pend.F(x, dt);
h      = @(x) sin(x(1));
H      = @(x) [cos(x(1)), 0];
qc     = 1;
Q      = [qc * dt^3 / 3, qc * dt^2 / 2; qc * dt^2 / 2, qc * dt];
R      = 0.1;

dsPend  = dynamicalSystem(2, 1, 'x0', [-4;-4], 1e2, 'evolution', f, F, Q, 'emission', h, H, R, 'data', 400);


%%
dsPend  = dsPend.filterExtended;
dsPend  = dsPend.smoothExtended;
dsPend  = dsPend.save('original-params');
gui.posteriorGaussGUI(dsPend, 'original-params', []);