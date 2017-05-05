function filterStepIONLDS(obj, yDim, varargin)
    % plotStep2D(obj, posteriorType)
    %
    % INPUTS:
    % obj           - a dynamicalSystems object.
    % posteriorType - either 'filter' or 'smooth'
    %
    % DESCRIPTION:
    % Plot prediction and correction posteriors at each time step. Useful
    % demonstration of how kalman filter works. Note: coded assuming 2D
    % state space: if not will error.
    
    assert(isa(obj, 'ds.dynamicalSystem'), 'first argument must be a dynamicalSystem object');
    assert(isnumeric(yDim) && utils.is.scalarint(yDim, 1) && yDim <= obj.d.y, 'yDim must be a scalar in 1, ..., %d', obj.d.y);
    
    assert(isa(obj, 'ds.ionlds'), 'IONLDS stuff currently implemented and obj is not an IONLDS. Remove this to continue.');
    obj.ensureInference('filterStepIONLDS (plot)', 'filter');
    assert(isfield(obj.infer.smooth, 'mu'), 'need to run smoother..');
    
    optsDefault.showPrincipalAxis = false;
    optsDefault.showNonlinearity  = false;
    optsDefault.fileName          = '';
    optsDefault.maxT              = obj.d.T;
    
    opts         = utils.base.processVarargin(varargin, optsDefault);
    
    opts.createAVI = false;
    if ~isempty(opts.fileName)
        opts.createAVI = true; 
        v = VideoWriter(opts.fileName); %, 'Uncompressed AVI');
        v.FrameRate = 0.75;
        open(v);
    end
    
    figure('Position', [0 0 800 600]);
    hold on;

    usePosterior = obj.infer.filter;
    mu      = usePosterior.mu;
    sigma   = usePosterior.sigma;
    
    m       = obj.par.x0.mu;
    P       = obj.par.x0.sigma;
    u_t     = [];
    
    C    = obj.par.emiNLParams.C(yDim,:);
    bias = obj.par.emiNLParams.bias(yDim);
    
    xrng = quantile(C * obj.infer.smooth.mu + bias, [0.05, 0.95]);
    xrng(1)   = xrng(1) - abs(xrng(1))*0.05;
    xrng(2)   = xrng(2) + abs(xrng(2))*0.05;
    xlim(xrng);
    
    xvals = linspace(xrng(1), xrng(2), 100);
    yvals = NaN(1,100); 
    for ii = 1:100; yvals(:,ii) =  obj.par.emiNLParams.bspl.functionEval(xvals(ii), obj.par.emiNLParams.eta(yDim,:)); end  % xvals includes C*(x)+b
    yrng  = quantile(yvals, [0.05, 0.95]);
    yrng(1)   = yrng(1) - abs(yrng(1))*0.05;
    yrng(2)   = yrng(2) + abs(yrng(2))*0.05;
    ylim(yrng);
    
    utpar       = struct('alpha', 1, 'beta', 0, 'kappa', 0);
    pType       = (~obj.evoLinear)*2;
    eType       = (~obj.emiLinear)*2;
    parPredict  = ds.internal.getParams(obj, 1, pType);
    parEmission = ds.internal.getParams(obj, 2, eType);
    
    if opts.showNonlinearity
        plot(xvals, yvals, '-', 'Color', [177,177,177]./255);
    end
    
    for tt = 1:opts.maxT
        title(sprintf('Iteration %d', tt));
        % one step ahead
        if obj.hasControl(1); u_t = obj.u(:,tt); end
        [m_minus, P_minus] = ds.utils.assumedDensityTform(parPredict, m, P, u_t, pType, utpar);
        
        if tt > 1
            delete(pt2c);
            delete(pt2a); delete(pt2b);
            delete(pt2a2); delete(pt2b2);
            delete(pt3a); delete(pt3b);
            delete(pt4a);
        end
        
%         if tt == 25
%             stophere = 1;
%         end
%         
        % specific stuff to IONLDS
        
        Cm     = C * m_minus + bias;
        Cv     = C * P_minus * C';
        
        % draw marginal x transformed to y space
        yrng = ylim; 
        pt1b = line([C*m+bias, C*m+bias], [yrng(1), yrng(2)], 'LineStyle', ':', 'Color', 'b');    % (prev) posterior mean
        pt1c = line([Cm, Cm], [yrng(1), yrng(2)], 'LineStyle', '-.', 'Color', 'r');     % prior mean
        pt1a = plotGauss1D(Cm, Cv, 'r-', 0.3);  % mu, sigma, style, pctScreen           % prior variance
        
        if opts.createAVI
            writeVideo(v, getframe(gcf));
        else
            pause;
        end
        
        % draw predictive covariance for (x, y).
        if obj.hasControl(2); u_t = obj.u(:,tt); else u_t = []; end
        [m_y, P, Cxy] = ds.utils.assumedDensityTform(parEmission, m_minus, P_minus, u_t, eType, utpar);
        delete(pt1b);
        crosscov = [Cv, C*Cxy(:,yDim); Cxy(:,yDim)'*C', P(yDim,yDim)];
        [pt2a, pt2b]   = plotGaussLevel([Cm; m_y(yDim)], crosscov, 'm', ':');   % predictive covariance
        [pt2a2, pt2b2] = plotGaussLevel([Cm; m_y(yDim)], 4*crosscov, 'm', ':');   % predictive covariance
        if opts.createAVI
            writeVideo(v, getframe(gcf));
        else
            pause;
        end
        
        xrng = xlim;
        pt2c = plot([xrng(1),xrng(2)], [obj.y(yDim,tt),obj.y(yDim,tt)], 'k-');   % y
        if opts.createAVI
            writeVideo(v, getframe(gcf));
        else
            pause;
        end
        
        % correct with next y value
        delete(pt1a);
        delete(pt1c);
        yrng = ylim;
        pt3a = plotGauss1D(C*mu(:,tt)+bias, (C*sigma{tt}*C'), 'b:', 0.3);           % posterior variance
        pt3b = plot([C*mu(:,tt)+bias, C*mu(:,tt)+bias], [yrng(1), yrng(2)], 'b-');     % posterior mean
        
        pt4a = plotGauss1D(C*obj.infer.smooth.mu(:,tt)+bias, (C*obj.infer.smooth.sigma{tt}*C'), 'c:', 0.3);
        if opts.createAVI
            writeVideo(v, getframe(gcf));
        else
            pause;
        end
        
        m = mu(:,tt);
        P = sigma{tt};
        
    end
    
    if opts.createAVI
        close(v);
    end
end


function [pt1, pt2] = plotGaussLevel(mu, sigma, col, style, doVar)
    % plot mean and level curve (1 sd) of Gaussian.
    % Optionally return the pointers of the mean and var lines.
    
    if nargin < 5 || isempty(doVar)
        doVar = true;
    end
    
    origHoldOn = ishold;
    hold on
    pt1 = plot(mu(1), mu(2), [col, '+']);
    
    
    if doVar
        lc = utils.plot.gaussian2DLevelCurve(1, mu, sigma, 100);
        pt2 = plot(lc(:,1), lc(:,2), [col, style]);
    else
        pt2 = [];
    end
    
    if ~origHoldOn 
        hold off;
    end
end

function pt1 = plotGauss1D(mu, varnce, style, pctScreen)
    % plot pdf of 1D gaussian on abscissa
    % Optionally return the pointer of line.
    
    sigma   = sqrt(varnce);
    origHoldOn = ishold;
    hold on
    xrng    = xlim;
    yrng    = ylim;
    gridpts = linspace(xrng(1), xrng(2), 100);
    depgrid = mu + (-5:0.25:5)*sigma;
    depgrid(depgrid > xrng(2)) = [];
    depgrid(depgrid < xrng(1)) = [];
    gridpts = union(gridpts, depgrid);
    
    y   = normpdf(gridpts, mu, sigma);
    y   = (y ./ max(y)).*diff(yrng)*pctScreen;
    y   = y + yrng(1);
    
    pt1 = plot(gridpts, y, style);
    
    if ~origHoldOn 
        hold off;
    end
end