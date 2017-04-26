function smse = getKstepSMSE(obj, k, varargin)
    
    optsDefault = struct('average', true, 'standardize', 'allVar');
    opts     = utils.base.processVarargin(varargin, optsDefault);
    
    assert(ismember(opts.standardize, {'allVar', 'remVar'}), ['''standardize'' must be either ' ...
        ,'''allVar'' or ''remVar'', Current setting is %s'], opts.standardize)
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
        if strcmp(opts.standardize, 'allVar')
            smse(:, nn) = nanmean((y{nn}(:,ix) - predvals{nn}).^2,2)./nanvar(y{nn}, [], 2);   % smse by channel (STNDD BY VAR(*all* Y))
        elseif strcmp(opts.standardize, 'remVar')
            yvar        = nanvar(y{nn}(:,ix), [], 2);
            smse(:, nn) = nanmean((y{nn}(:,ix) - predvals{nn}).^2,2)./yvar;   % smse by channel (STNDD BY VAR(remaining Y))
        else
            error('invalid standardize choice');
        end
    end
    
    if opts.average
        smse = nanmean(smse, 1);    % still need nanmean in case all of one channel is obscured
    end
end