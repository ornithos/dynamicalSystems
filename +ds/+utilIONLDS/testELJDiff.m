function out = testELJDiff(obj, test, etay, etax)
    
%     tmpobj  = obj.copy;
    eta     = obj.par.emiNLParams.eta;
    obj.par.emiNLParams.eta(etay, etax) = obj.par.emiNLParams.eta(etay, etax) + test;
    
    out     = obj.expLogJoint;
    obj.par.emiNLParams.eta = eta;
end
    
    