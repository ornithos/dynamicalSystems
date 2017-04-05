function out = tmpllhcovfn(obj, x)
    y         = zeros(3,3);
    y(1,1)    = x(1);
    y(1:2,2)  = x(2:3);
    y(1:3,3)  = x(4:6);
%     y = zeros(2,2);
%     y(1,1) = x(1);
%     y(1:2,2) = x(2:3);
    obj.par.R = y'*y;
    out = -obj.expLogJoint;
end