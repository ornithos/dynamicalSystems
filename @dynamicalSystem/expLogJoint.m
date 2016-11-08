function q = expLogJoint(obj)

    s    = obj.suffStats;
    A    = obj.par.A;
    H    = obj.par.H;
    T    = obj.d.T;
    q    = -0.5*T*log(det(2*pi*obj.par.Q));
    q    = q -0.5*T*log(det(2*pi*obj.par.R));
    q    = q -0.5*T*trace(inv(obj.par.Q)*(s.SIGMA - s.C*A' - A*s.C' + A*s.PHI*A'));
    q    = q -0.5*T*trace(inv(obj.par.R)*(s.D - s.B*H' - H*s.B' + H*s.SIGMA*H'));
end