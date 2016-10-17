% initial "circular motion" trial
tmp = dynamicalSystem([], [cos(0.03) -sin(0.03); sin(0.03) cos(0.03)], eye(2), 0.5*eye(2), 0.2*eye(2), 50);
tmp = tmp.useInputParameters;
tmp = tmp.posteriorSmooth;
gui.posteriorGaussGUI(tmp);

%%
% Newtonian dynamics (4D state space to keep linear system with acc/v)
deltaT              = 0.5;
transition          = eye(4);
transition(1,3)     = deltaT;
transition(2,4)     = deltaT;

emission            = [1 0 0 0; 0 1 0 0];
Q                   = 0.01*eye(4);
Q(3,3)              = 0.2;
Q(4,4)              = 0.2;
% x0                  = struct('mu', zeros(4,1), 'sigma', eye(4)*1e2);
x0                  = [];
dsNewton  = dynamicalSystem(x0, transition, emission, Q, 5*eye(2), 150);
dsNewton  = dsNewton.useInputParameters;
dsNewton  = dsNewton.posteriorSmooth;
opts = struct('maxiter', 500, 'epsilon', 1e-3);
[dsNewton, llh] = dsNewton.parameterLearningEM(opts);

%%
dsNewtonPlot = dsNewton;
dsNewtonPlot.x = dsNewtonPlot.x(1:2,:);
H = dsNewtonPlot.H;
dsNewtonPlot.posterior.filter.mu = H * dsNewtonPlot.posterior.filter.mu;
for tt = 1:dsNewtonPlot.T
    dsNewtonPlot.posterior.filter.sigma{tt} = H * dsNewtonPlot.posterior.filter.sigma{tt} * H';
end
dsNewtonPlot.posterior.smooth.mu = H * dsNewtonPlot.posterior.smooth.mu;
for tt = 1:dsNewtonPlot.T
    dsNewtonPlot.posterior.smooth.sigma{tt} = H * dsNewtonPlot.posterior.smooth.sigma{tt} * H';
end
H = dsNewtonPlot.inH;
dsNewtonPlot.posterior.inFilter.mu = H * dsNewtonPlot.posterior.inFilter.mu;
for tt = 1:dsNewtonPlot.T
    dsNewtonPlot.posterior.inFilter.sigma{tt} = H * dsNewtonPlot.posterior.inFilter.sigma{tt} * H';
end
dsNewtonPlot.posterior.inSmooth.mu = H * dsNewtonPlot.posterior.inSmooth.mu;
for tt = 1:dsNewtonPlot.T
    dsNewtonPlot.posterior.inSmooth.sigma{tt} = H * dsNewtonPlot.posterior.inSmooth.sigma{tt} * H';
end
gui.posteriorGaussGUI(dsNewtonPlot);