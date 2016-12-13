function out = testELJDiff2(obj, test, etay, etax)
    
%     tmpobj  = obj.copy;
    eta     = obj.par.emiNLParams.eta;
    obj.par.emiNLParams.eta(etay, etax) = obj.par.emiNLParams.eta(etay, etax) + test;
    
    [~,M2] = ds.utilIONLDS.utTransform_ymHx(obj);
    out     = -0.5*obj.d.T*trace(inv(obj.par.R)*M2);
    obj.par.emiNLParams.eta = eta;
end
    
    