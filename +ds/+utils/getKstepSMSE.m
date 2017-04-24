function smse = getKstepSMSE(obj, k, varargin)
    
    optsDefault = struct('average', true);
    opts     = utils.base.processVarargin(varargin, optsDefault);
    
    predvals = obj.getPredictedValues(k);
    
    % data may be cell if dsBatch
    if isa(obj, 'ds.dynamicalSystemBatch')
        Nmodels = obj.d.n;
        y       = obj.y;
%         d       = obj.ambientDimension;
    else
        Nmodels = 1;
        predvals = {predvals}; y = {obj.y};
%         d       = obj.d.y;
    end
    
    smse = zeros(obj.d.y, Nmodels);
    
    for nn = 1:Nmodels
        ix          = obj.d.T(nn) - size(predvals{nn},2)+1:obj.d.T(nn);
        smse(:, nn) = nanmean((y{nn}(:,ix) - predvals{nn}).^2,2)./nanvar(y{nn}, [], 2);   % smse by channel (STNDD BY VAR(*all* Y))
    end
    
    if opts.average
        smse = nanmean(smse, 1);    % still need nanmean in case all of one channel is obscured
    end
end