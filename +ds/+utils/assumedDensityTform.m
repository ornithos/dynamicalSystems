function [m, P, C] = assumedDensityTform(pars, m, P, u, type, utpar)
    % [m, P, C] = assumedDensityTform(pars, m, P, u, type, utpar)
    % Transform Gaussian density with moments (m,P) through function
    % defined either by linearity par.A, or nonlinearity par.f with
    % gradient par.Df. Nonlinear transformation performed by
    %
    % type = 0: Kalman (Linear)
    % type = 1: EKF
    % type = 2: UKF
    %
    % It is assumed that there is additive noise in the transformation
    % which is zero mean white with covariance par.Q.
    % Control inputs are given 
    
    switch type
        case 0
            if isempty(u)
                m = pars.A * m + pars.bias;
            else
                m = pars.A * m + pars.B * u + pars.bias;
            end
            C         = P * pars.A';
            P         = pars.A * C + pars.Q;
        case 1
            if isempty(u)
                F         = pars.Df(m);
                m         = pars.f(m);
            else
                F         = pars.Df(m,u);
                m         = pars.f(m,u);
            end
            C         = P * F';
            P         = F * P * F' + pars.Q;
        case 2
            if isempty(u)
                f         = @(x) pars.f(x);
            else
                f         = @(x) pars.f(x,repmat(u, 1, 2*numel(m)+1));
            end
            [m, P, C] = ds.utils.unscentedTransform(f, m, P, ...
                            utpar.alpha, utpar.beta, utpar.kappa);
            P         = P + pars.Q;
        otherwise
            error('Unknown assumed density type. Try 0=Kalman,1=EKF,2=UKF');
    end
end