function [a, D] = expLogJoint(obj)

    tmpobj = obj;
    opts = struct; %struct('bIgnoreHash', true);
    s    = tmpobj.suffStats(opts);
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
        ctrlAdd = obj.par.B*s.UU*obj.par.B' - Bum - Bum' - BuAm_m - BuAm_m';
    end
    
    q(3)    = -0.5*T*trace((tmpobj.par.Q)\(s.SIGMA - s.C*A' - A*s.C' + A*s.PHI*A' + ctrlAdd));
    
    if ~obj.hasControl(2)
        ctrlAdd = zeros(obj.d.y);
    else
        Cuy     = obj.par.C * s.YU';
        Cum_H   = obj.par.C * s.XU' * obj.par.H';
        ctrlAdd = obj.par.C*s.UU*obj.par.C' - Cuy - Cuy' - Cum_H - Cum_H';
    end
    
    q(4)    = -0.5*T*trace((tmpobj.par.R)\(s.D - s.B*H' - H*s.B' + H*s.SIGMA*H' + ctrlAdd));
    a       = sum(q);
    
    if nargout > 1
        D = 0;
    end
end