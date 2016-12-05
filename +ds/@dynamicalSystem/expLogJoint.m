function [a, q] = expLogJoint(obj)

    tmpobj = obj;
    opts = struct; %struct('bIgnoreHash', true);
    s    = tmpobj.suffStats(opts);
    A    = tmpobj.par.A;
    H    = tmpobj.par.H;
    T    = tmpobj.d.T;
    q    = NaN(4,1);
    q(1)    = -0.5*T*log(det(2*pi*tmpobj.par.Q));
    q(2)    = -0.5*T*log(det(2*pi*tmpobj.par.R));
    q(3)    = -0.5*T*trace(inv(tmpobj.par.Q)*(s.SIGMA - s.C*A' - A*s.C' + A*s.PHI*A'));
    q(4)    = -0.5*T*trace(inv(tmpobj.par.R)*(s.D - s.B*H' - H*s.B' + H*s.SIGMA*H'));
    a       = sum(q);
end