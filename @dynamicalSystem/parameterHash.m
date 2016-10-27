function out = parameterHash(obj)
    hashes = cell(4,1);
    if obj.evoLinear
        hashes{1} = utils.base.DataHash(obj.A);
    else
        if obj.evoNLhasParams
            hashes{1} = utils.base.DataHash(obj.evoNLParams);
        else
            hashes{2} = '';
        end
    end

    if obj.emiLinear
        hashes{2} = utils.base.DataHash(obj.H);
    else
        if obj.emiNLhasParams
            hashes{2} = utils.base.DataHash(obj.emiNLParams);
        else
            hashes{2} = '';
        end
    end

    hashes{3}   = utils.base.DataHash(obj.Q);
    hashes{4}   = utils.base.DataHash(obj.R);
    out         = utils.base.DataHash(strjoin(hashes, '')); 
end