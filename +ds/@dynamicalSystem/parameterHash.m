function out = parameterHash(obj)
    hashes = cell(4,1);
    if obj.evoLinear
        hashes{1} = utils.base.DataHash(obj.par.A);
    else
        if obj.evoNLhasParams
            hashes{1} = utils.base.DataHash(obj.par.evoNLParams);
        else
            hashes{1} = '';
        end
    end

    if obj.emiLinear
        hashes{2} = utils.base.DataHash(obj.par.H);
    else
        if obj.emiNLhasParams
            hashes{2} = utils.base.DataHash(obj.par.emiNLParams);
        else
            hashes{2} = '';
        end
    end

    hashes{3}   = utils.base.DataHash(obj.par.Q);
    hashes{4}   = utils.base.DataHash(obj.par.R);
    out         = utils.base.DataHash(strjoin(hashes, '')); 
end