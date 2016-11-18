function [m, P, lam] = filterStep(obj, m, P, typePredict, parPredict, typeUpdate, parUpdate, utpar)
    %filterStep(obj, m, P, parPredict, parUpdate)
    %  internal filterStep function. Performs filter for specified mean and
    %  covariance using dynamics of the parPredict step and parUpdate step.
    %
    % INPUTS:
    %  obj         - a dynamicalSystems object
    %  m           - the mean at time step (t-1)
    %  P           - the covariance at time step (t-1)
    %  typePredict - the type of filter to use for predict/transition step
    %  parPredict  - the parameters of the predict/transition step.
    %  typeUpdate  - the type of filter to use for update/emission step
    %  parUpdate   - the parameters of the update/emission step.
    %  utpar       - unscented parameter struct (if applicable)
    %
    %  Note that types are enumerated: (0) = Linear; (1) = EKF; (2) = UKF.
    %  
    %  OUTPUTS
    %  m          - mean at time step (t)
    %  P          - covariance at time step (t)
    %  lam        - eigenvalues of matrix S
        
    % Prediction step
    if any(obj.hasControl)
        u_t = obj.u(:,tt);
    else
        u_t = [];
    end
    [m_minus, P_minus, ~] = ds.utils.assumedDensityTform(parPredict, m, P, u_t, typePredict, utpar);

    % Filter step
    [m_y, S, covxy]   = ds.utils.assumedDensityTform(parUpdate, m_minus, P_minus, u_t, typeUpdate, utpar);
    [Sinv, lam]       = utils.math.pinvAndEig(S, 1e-12);
    K                 = covxy * Sinv;
    deltaY            = obj.y(:,tt) - m_y;
    m                 = m_minus + K*deltaY;
    P                 = P_minus - K * S * K';
end