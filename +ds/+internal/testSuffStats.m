function out = testSuffStats(obj, n, varargin)
    % monte carlo estimates of the sufficient statistics calculated in
    % suffStats.m
    
    assert(isa(obj, 'ds.dynamicalSystem'), 'object is not a dynamicalSystems object');
    if nargin < 2
        n = 50000;
    end
    assert(utils.is.scalarint(n, 1), 'n must be a positive integer');
    
    optsDefault.SIGMA  = true;
    optsDefault.B      = true;
    optsDefault.D      = true;
    
    optsDefault.center = false;
    optsDefault.var    = false;
    opts               = utils.base.processVarargin(varargin, optsDefault);
    
    pb          = utils.base.objProgressBar;
    updateEvery = 100;
    
    dx  = obj.d.x;
    dy  = obj.d.y;
    T   = obj.d.T;
    out = struct;
    
    %% Generate SIGMA
    if opts.SIGMA
        
    pb.newProgressBar('Generating SIGMA: ', 30, 'showElapsed', true); 
    SIGMA  = zeros(dx, dx); 
    tmpx   = zeros(dx, T); 
    for ii = 1:n
        for jj = 1:T 
            tmpx(:,jj) = obj.infer.smooth.mu(:,jj) + chol(obj.infer.smooth.sigma{jj})' * randn(dx, 1); 
        end
        if opts.center
            tmpx  = tmpx - obj.infer.smooth.mu(:,1:T);
        end
        
        SIGMA = SIGMA + (tmpx*tmpx'./T); 
        if mod(ii, updateEvery)==0
            pb.print(ii/n); 
        end
    end
    SIGMA = SIGMA./n;
    pb.finish;
    
    out.SIGMA = SIGMA;
    end
    %% Generate D
    if opts.D
        
    pb.newProgressBar('Generating D:     ', 30, 'showElapsed', true); 
    D    = zeros(dy, dy);
    Dvar = zeros(dy,dy);
    for ii = 1:n
        
        %tmpy = obj.impute_y('sample', true); 
        %%% COPIED FROM IMPUTE_Y TO AVOID OVERHEADS:
        % ---------------------------------------------------------------
        tmpy       = obj.y;
        missing    = isnan(tmpy);
        anyMissing = any(missing, 1);
        xhat       = obj.infer.smooth.mu;
        P          = obj.infer.smooth.sigma;
        % precompute cholesky if need samples
        cholR      = chol(obj.par.R);

        % impute missing values
        for tt = find(anyMissing)
%             u_t     = [];
%             if obj.hasControl(2); u_t = obj.u(:,tt); end
            tmpx     = xhat(:,tt) + chol(P{tt})' * randn(obj.d.x, 1);
            yhat_tt  = obj.par.H * tmpx + cholR' * randn(obj.d.y, 1);
            mask        = missing(:,tt);
            tmpy(mask,tt)  = yhat_tt(mask);
        end
        % ---------------------------------------------------------------
        
        x = (tmpy*tmpy'./T);
        
        if opts.var && ii > 1
            Dvar = Dvar + (((ii-1) * x - D).^2)./((ii-1)*ii);
        end
        
        D = D + x; 
        if mod(ii, updateEvery)==0
            pb.print(ii/n); 
        end
    end
    D    = D./n;
    Dvar = Dvar./(n+1);
    pb.finish;
    
    out.D     = D;
    if opts.var; out.Dvar = Dvar; end

    end
    %% Generate B
    if opts.B
    
    pb.newProgressBar('Generating B:     ', 30, 'showElapsed', true); 
    B      = zeros(dy, dx);
    Bvar   = zeros(dy, dx);
    cholR  = chol(obj.par.R);
    tmpx   = zeros(dx, T); 
    for ii = 1:n
        for jj = 1:T 
            tmpx(:,jj) = obj.infer.smooth.mu(:,jj) + chol(obj.infer.smooth.sigma{jj})' * randn(dx, 1); 
        end
        ysmp = obj.par.H * tmpx + cholR'*randn(dy, T);
        tmpy = obj.y; 
        tmpy(isnan(tmpy)) = ysmp(isnan(tmpy)); 
        x    = (tmpy*tmpx'./obj.d.T);
        
        if opts.var && ii > 1
            Bvar = Bvar + (((ii-1) * x - B).^2)./((ii-1)*ii);
        end
        
        B = B + x;
        if mod(ii, updateEvery)==0
            pb.print(ii/n); 
        end
    end
    B    = B./n;
    Bvar = Bvar./(n+1);
    pb.finish;
    
    out.B     = B;
    if opts.var; out.Bvar = Bvar; end
    end
    %%
    
    
    
    
end