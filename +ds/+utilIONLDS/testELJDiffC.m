function out = testELJDiffC(obj, test, jj, ii)
    
%     tmpobj  = obj.copy;
    C     = obj.par.emiNLParams.C;
    obj.par.emiNLParams.C(jj, ii) = obj.par.emiNLParams.C(jj, ii) + test;
    
    out     = obj.expLogJoint;
    obj.par.emiNLParams.C = C;
end
    