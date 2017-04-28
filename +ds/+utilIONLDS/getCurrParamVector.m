function out =  getCurrParamVector(obj, etaMask, logSpace)
    % extract current emission function parameters from ionlds object -->
    % single vector for parameter learning.
    emiParams = obj.par.emiNLParams;
    cEta      = emiParams.eta(:, etaMask);
    cEta      = [cEta(:,1), diff(cEta(:,1:end), 1, 2)];
    if logSpace
        cEta = log(cEta); 
    end
    out       = [cEta(:); emiParams.C(:); emiParams.bias(:)];
end