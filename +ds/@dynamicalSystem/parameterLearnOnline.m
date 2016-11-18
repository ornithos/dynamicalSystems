function [obj, llhHist] = parameterLearnOnline(obj, fType, opts, utpar)

    %% Process arguments
    if nargin < 4 || isempty(utpar)
        utpar = struct;
    end
    utparDefault = struct('alpha', 1, 'beta', 0, 'kappa', 0);
    utpar        = utils.struct.structCoalesce(utpar, utparDefault);

    if nargin < 3 || isempty(opts)
        opts = struct;
    end

    optsDefault     = struct('epsilon', 1e-3, 'maxiter', 300, 'type', 'full-sweep', ...
                        'typeOpt', 1, 'minT', 20,  'verbose', true, 'dbg', false);
    optsDefault     = utils.base.parse_argumentlist(obj.opts, optsDefault, false);      % bring in global opts
    opts            = utils.base.parse_argumentlist(opts, optsDefault, false);          % add user specified opts.

    emOpts          = utils.struct.structColFn(opts, {'epsilon', 'maxiter', 'verbose', 'dbg'}, 'keep');
    fOpts           = struct('bIgnoreHash', true);
    emOpts          = utils.struct.structCoalesce(emOpts, fOpts);
    emOpts.verbose  = false;

    if nargin < 2 || isempty(fType)
        if obj.evoLinear && obj.emiLinear
            fType = 'kal';
        else
            error(['filter type (arg 1) not specified for nonlinear model. ', ...
                'Please specify algorithm.']);
        end
    end

    bNumeric     = true;
    inpFtype     = fType;
    fType        = ds.utils.filterTypeLookup(fType, bNumeric) - 1;
    if fType>0 && obj.evoLinear && obj.emiLinear
        warning(['Model is linear: exact inference will be performed using ', ...
                'Kalman equations rather than %s type requested'], utils.filterTypeLookup(inpFtype));
    end

    fType1          = fType;
    fType2          = fType;
    if obj.evoLinear; fType1 = 0; end
    if obj.emiLinear; fType2 = 0; end

    obj.validationInference;  % ensure able to do inference    
    
    
    %% Parameter set-up
    startTime       = tic;
    parPredict      = ds.internal.getParams(obj, 1, fType1);
    parUpdate       = ds.internal.getParams(obj, 2, fType2);
    llhHist         = zeros(obj.d.T, 3);
    d               = obj.d.y;
    llh             = 0;
    
    filterMu        = zeros(obj.d.x, obj.d.T);
    filterSigma     = cell(obj.d.T, 1);
    filterG         = cell(obj.d.T, 1);
    
    % set up object on which to perform online estimation
    tmpobj                      = obj;
    tmpobj.infer.filter.mu      = filterMu;
    tmpobj.infer.filter.sigma   = filterSigma;
    tmpobj.infer.llh            = NaN;
    tmpobj.infer.fType          = inpFtype;
    tmpobj.infer.sType          = inpFtype;
    tmpobj.infer.filter.utpar   = utpar;
    tmpobj.d.T                  = 0;
    
    
    % initialise t-1 = 0 values to prior
%     [m, P, covttp1] = ds.utils.assumedDensityTform(parPredict, obj.par.x0.mu, obj.par.x0.sigma, zeros(obj.d.x,1), fType1, utpar);
    m               = obj.par.x0.mu;
    P               = obj.par.x0.sigma;
    
%     if strcmpi(opts.type, 'filter-only')
%         G               = covttp1 / P;
%         m               = obj.par.x0.mu + G * (obj.par.x0.mu - m_minus);
%         P               = fP_tt + G * (obj.par.x0.sigma - P_minus) * G';
%     end
    pb = utils.base.objProgressBar;
    pb.newProgressBar('EM Iterations: ', 25, true, true);
    for tt = 1:obj.d.T
        
        %% Filter
        % Predict step
        if any(obj.hasControl)
            u_t = obj.u(:,tt);
        else
            u_t = [];
        end
        [m_minus, P_minus, covttp1] = ds.utils.assumedDensityTform(parPredict, m, P, u_t, fType1, utpar);

        % Update step
        [m_y, S, covxy]   = ds.utils.assumedDensityTform(parUpdate, m_minus, P_minus, u_t, fType2, utpar);
        [Sinv, lam]       = utils.math.pinvAndEig(S, 1e-12);
        K                 = covxy * Sinv;
        deltaY            = obj.y(:,tt) - m_y;
        m                 = m_minus + K*deltaY;
        P                 = P_minus - K * S * K';
        
        llh               = llh - 0.5*d*log(2*pi) - 0.5*sum(log(lam)) - 0.5*deltaY'*Sinv*deltaY;
        filterG{tt}       = covttp1 / (P_minus);
        
        % save
        tmpobj.infer.filter.mu(:,tt)  = m;
        tmpobj.infer.filter.sigma{tt} = P;
        tmpobj.y                      = obj.y(:,1:tt);
        tmpobj.d.T                    = tt;
        
        tmpobjCopy                    = tmpobj;  % for debugging
        
        %% Smooth
        if tt >= opts.minT
            switch upper(opts.type)
                case 'FILTER-ONLY'
                    error('This procedure is currently mothballed');
                    tmpobj.infer.smooth.mu    = tmpobj.infer.filter.mu;
                    tmpobj.infer.smooth.sigma = tmpobj.infer.filter.sigma;
                    tmpobj.infer.smooth.G     = filterG;
                    tmpobj.infer.smooth.x0    = struct('mu', m, 'sigma', P, 'G', G);
                    tmpobj.infer.smooth.utpar = utpar;
                    tmpobj.infer.sType        = inpStype;

                    [tmpobj, llhEM, iters] = parameterLearningEM(tmpobj, emOpts); % <-- this will re-estimate the filter/smooth estimates!
                case 'FULL-SWEEP'
                    [tmpobj, llhEM, iters] = parameterLearningEM(tmpobj, emOpts);
                    pb.print(tt/obj.d.T);
                case 'ALT-MURPHY'
                    [tmpobj] = alex_learn_kalman(tmpobj, emOpts.maxiter);
                    llhEM    = tmpobj.logLikelihood;
                    iters    = emOpts.maxiter;
                otherwise
                    error('online learning type must be one of ''filter-only'', ''full-sweep''');
            end
            llhEM            = llhEM(end);
            tmpobj.infer.llh = llhEM;
            llhHist(tt, :)   = [llhEM, iters, toc(startTime)];
            
            m           = tmpobj.infer.filter.mu(:,tt);
            P           = tmpobj.infer.filter.sigma{tt};
            parPredict  = ds.internal.getParams(obj, 1, fType1);
            parUpdate   = ds.internal.getParams(obj, 2, fType2);
        end
    end
    pb.finish;
    llhHist = llhHist(1:tt,:);
end