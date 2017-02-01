function [a, q] = expLogJoint(obj, varargin)
    
    % expLogJoint
    % Method of a dynamical systems object. Calculates the expected log
    % joint **given the current posterior**, <log(p(x,z))>
    %
    % Usage:
    %  val           = obj.ExpLogJoint
    %          - Outputs the value of <log(p(x,z))>.
    %
    %  [val, decomp] = obj.ExpLogJoint
    %          - second output gives decomposition of value as normalising
    %            constants (1), (2), and quadratic forms of both distns as
    %            (3), (4).
    %
    %  val           = obj.ExpLogJoint('freeEnergy', true, ...)
    %          - Outputs the free energy <log(p(x,z))> + H(q).
    %
    %  val           = obj.ExpLogJoint(..., 'recalcPosterior', true)
    %          - Recalculates the posterior before calculating expLogJoint.
    %            Note that the new posterior will not be saved to the ds
    %            object, ie. it will not affect the current state of
    %            inference.
    %
    
    optsDefault = struct('freeEnergy', false, 'recalcPosterior', false);
    opts        = struct;
    narargin = numel(varargin);
    assert(mod(narargin,2)==0, 'arguments to expLogJoint must come in name-value pairs');
    for kk = 0:(narargin/2-1)
        assert(ischar(varargin{2*kk+1}), 'odd numbered arguments must be names (character strings)');
        opts.(varargin{2*kk+1}) = varargin{2*kk+2};
    end
    opts        = utils.base.parse_argumentlist(opts, optsDefault);
    ssopts      = struct('bIgnoreHash', ~opts.recalcPosterior);
    
    %% 
    tmpobj = obj.copy;
    s    = tmpobj.suffStats(ssopts);
    A    = tmpobj.par.A;
    H    = tmpobj.par.H;
    T    = tmpobj.d.T;
    q    = NaN(4,1);
    q(1)    = -0.5*T*log(det(2*pi*tmpobj.par.Q));
    q(2)    = -0.5*T*log(det(2*pi*tmpobj.par.R));
    
    if ~obj.hasControl(1)
        ctrlAdd = zeros(obj.d.x);
    else
        Bum     = obj.par.B * s.XU';
        BuAm_m  = obj.par.B * s.Xm1_U' * obj.par.A';
        ctrlAdd = obj.par.B*s.UU*obj.par.B' - Bum - Bum' + BuAm_m + BuAm_m';
    end
    
    q(3)    = -0.5*T*trace((tmpobj.par.Q)\(s.SIGMA - s.C*A' - A*s.C' + A*s.PHI*A' + ctrlAdd));
    
    if ~obj.hasControl(2)
        ctrlAdd = zeros(obj.d.y);
    else
        Cuy     = obj.par.C * s.YU';
        Cum_H   = obj.par.C * s.XU' * obj.par.H';
        ctrlAdd = obj.par.C*s.UU*obj.par.C' - Cuy - Cuy' + Cum_H + Cum_H';
    end
    
    q(4)    = -0.5*T*trace((tmpobj.par.R)\(s.D - s.B*H' - H*s.B' + H*s.SIGMA*H' + ctrlAdd));
    
    % x0
    parV0   = tmpobj.par.x0.sigma;
    parm0   = tmpobj.par.x0.mu;
    P0      = tmpobj.infer.smooth.x0.sigma;
    m0      = tmpobj.infer.smooth.x0.mu;
    q(5)    = -0.5*utils.math.logdet(2*pi*parV0) -0.5*trace((parV0)\(P0 + (m0-parm0)*(m0-parm0)'));
    
    %%
    if opts.freeEnergy
        % entropy calculation for free energy
        if obj.d.T*obj.d.x > 2000
            warning('Determinant of size %d x %d will be calculated...', obj.d.T*obj.d.x,obj.d.T*obj.d.x);
        end
        fullcov = ds.utils.fullJointCovariance(obj);
        detArg  = 2*pi*fullcov;
        L       = chol(detArg);
        q(6)    = 0.5*2*sum(log(diag(L))) + (obj.d.T+1)*obj.d.x./2;  % logdet = 2*sum(log(diag(chol(.))))
        % T + 1 since posterior from 0:T.
    end
    
    a       = sum(q);
end