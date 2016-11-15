function [f, g] = LDSoptimFn(theta, dsobj)
    dsobj = ds.utils.vecSetParams(dsobj, theta, true);  % final argument: project matrices onto PSD cone.
    dsobj = dsobj.filter('Linear',true);
    f     = -dsobj.infer.llh;
    
    if nargout > 1
        g = -ds.utils.vecGradient(dsobj);
    end
end