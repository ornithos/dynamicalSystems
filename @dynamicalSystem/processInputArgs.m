function obj = processInputArgs(obj, args)
% dynamicalSystem(dim_x, dim_y, 'evolution', {A || f || []}, {Df, Q},
   %                 'emission', {H || h || []}, {Dh, R},
   %                 'data', {y || T})
   
        nargs = numel(args);
         
        % update options
        optsDefault     = struct('warnings', true, 'verbose', false);
        if isstruct(args{nargs})
            obj.opts        = utils.base.parse_argumentlist(opts, optsDefault);
            nargs           = nargs - 1;
        else
            obj.opts        = optsDefault;
        end
        
        assert(utils.base.isscalarint(args{1}), 'dimx is not an integer scalar');
        assert(utils.base.isscalarint(args{2}), 'dimy is not an integer scalar');
        obj.d.x = args{1};
        obj.d.y = args{2};
        
        trackDone = false(4,1);
        xtrue     = [];
        ii = 3;
        assert(ischar(args{ii}), 'Unexpected argument in position %d: expected string', ii);
        curSection = internalProcessChar(args{ii}, ii, obj.opts.warnings);
        ii = ii + 1;
        
        while ii <= nargs
            jj = ii;
            while jj < nargs && ~ischar(args{jj+1})
                jj = jj + 1;
            end
            
%             if ischar(args{jj})
%                 error('Insufficient arguments for section (%s)', curSection);
%             end
            
            switch curSection
                case 'evo'
                    obj = internalProcessEvo(obj, args(ii:jj));
                    trackDone(1) = true;
                case 'emi'
                    obj = internalProcessEmi(obj, args(ii:jj));
                    trackDone(2) = true;
                case 'x0 '
                    obj = internalProcessX0(obj, args(ii:jj));
                    trackDone(3) = true;                
                case 'dat'
                    dataStuff = args(ii:jj);
                    trackDone(4) = true;
                case 'xtr'
                    assert(ii == jj, 'xtrue can have at most one argument');
                    xtrue   = args{ii};
                    % optional: no trackDone
                otherwise
                    error('Unknown section marker ''%d''', curSection);
            end
            
            if jj + 1 >= nargs
                break
            end
            
            % update section to process
            ii         = jj + 1;
            curSection = internalProcessChar(args{ii}, ii, obj.opts.warnings);
            ii         = ii + 1;
            while ischar(args{ii})
                warning('Input: ignoring character string ''%s''', args{ii+1});
                ii         = ii + 1;
            end
        end
        
        obj      = internalProcessDat(obj, dataStuff);
        if ~isempty(xtrue)
            assert(all(size(xtrue) == [obj.d.x, obj.d.T]), 'x must be matrix %d x %d', ...
                    obj.d.x, obj.d.T);
            obj.x = xtrue;
        end
        
        allSects = {'evo', 'emi', 'x0 ', 'dat'};
        if ~all(trackDone)
            error('No arguments for sections: %s', strjoin(allSects(~trackDone), ','));
        end
end

function arg = internalProcessChar(arg, ii, doWarning)
    orig = arg;
    if numel(arg < 3); arg = [arg, '   ']; arg = arg(1:3); end
    assert(any(strcmpi(arg(1:3), {'evo', 'emi', 'x0 ', 'dat', 'xtr'})), ...
        ['string argument #%d (%s) is not one of', ...
        ' the valid markers {x0, evolution, emission, data}'], ii, orig);
    argsFullnames = {'x0', 'evolution', 'emission', 'data', 'xtrue'};
    if doWarning && ~any(strcmpi(orig, argsFullnames))
        chosen = find(arg, {'evo', 'emi', 'x0 ', 'dat', 'xtr'});
        warning('arg %d: interpreted %s as %s', orig, argsFullnames{chosen});
    end
end

function obj = internalProcessX0(obj, arg)
    nargs       = numel(arg);
    switch nargs
        case 1
            if all(size(arg{1})==1)
                obj.par.x0 = struct('mu', zeros(obj.d.x,1), 'sigma', eye(obj.d.x)*arg{1});
            elseif all(size(arg{1})==[obj.d.x, obj.d.x])
                obj.par.x0 = struct('mu', zeros(obj.d.x,1), 'sigma', arg{1});
            else
                error('Invalid dimensions of x0 prior covariance');
            end
        case 2
            assert(numel(setxor(size(arg{1}), [obj.d.x,1]))==0, 'x0mu incompatible with size of latent space');
            assert(all(size(arg{2})==1) || all(size(arg{2})==[obj.d.x, obj.d.x]), ...
                'x0sigma incompatible with size of latent space');
            obj.par.x0 = struct('mu', arg{1}, 'sigma', []);
            if all(size(arg{2})==1)
                obj.par.x0.sigma = eye(obj.d.x)*arg{2};
            else
                obj.par.x0.sigma = arg{2};
            end
        otherwise
            error('Unexpected number of arguments for x0 (%d). Expected 1 or 2.', nargs);
    end
end

function obj = internalProcessEvo(obj, arg)
    nargs         = numel(arg);
    obj.evoLinear = false;
    %{A || f || []}, {Df, Q}
    if isnumeric(arg{1})
        assert(all(size(arg{1})==[obj.d.x, obj.d.x]), 'matrix A is not conformable to dim(x)');
        obj.par.A = arg{1};
        obj.evoLinear = true;
    elseif isa(arg{1}, 'function_handle')
        obj.par.f = arg{1};
        if obj.evoNLhasParams
            f = @(x) obj.par.f(x, obj.evoNLParams);
        else
            f = obj.par.f;
        end
        internalTestFun(f, 'f', false, obj.d.x, [obj.d.x,1]);
    else
        error('Unknown argument type (%s) in Evolution argument %d', class(arg{1}), 1);
    end
    
    exhaustNum = false;
    exhaustFn = false;
    for ii = 2:nargs
        if isnumeric(arg{ii}) && ~exhaustNum
            assert(all(size(arg{ii})==[obj.d.x, obj.d.x]), 'matrix Q is not conformable to dim(x)');
            obj.par.Q      = arg{ii};
            exhaustNum = true;
        elseif isa(arg{ii}, 'function_handle') && ~exhaustFn
            obj.par.Df = arg{ii};
            if obj.evoNLhasParams
                f = @(x) obj.par.Df(x, obj.evoNLParams);
            else
                f = obj.par.Df;
            end
            internalTestFun(f, 'Df', true, obj.d.x, [obj.d.x, obj.d.x]);
            exhaustFn = true;
        elseif isstruct(arg{ii})
            obj.evoNLParams = arg{ii};
            obj.evoNLhasParams = true;
        else
            error('Don''t know what to do with argument %d in Evolution section', ii);
        end
    end
    
    if ~exhaustNum
        error('No covariance matrix detected for evolution distribution');
    end
end
            
function obj = internalProcessEmi(obj, arg)
    nargs         = numel(arg);
    obj.emiLinear = false;
    %{H || h || []}, {Dh, R},
    if isnumeric(arg{1})
        assert(all(size(arg{1})==[obj.d.y, obj.d.x]), 'matrix H is not conformable to dim(x), dim(y)');
        obj.par.H = arg{1};
        obj.emiLinear = true;
    elseif isa(arg{1}, 'function_handle')
        obj.par.h = arg{1};
        if obj.emiNLhasParams
            f = @(x) obj.par.h(x, obj.emiNLParams);
        else
            f = obj.par.h;
        end
        internalTestFun(f, 'h', false, obj.d.x, [obj.d.y,1]);
    else
        error('Unknown argument type (%s) in Emission argument %d', class(arg{1}), 1);
    end
    
    exhaustNum = false;
    exhaustFn = false;
    for ii = 2:nargs
        if isnumeric(arg{ii}) && ~exhaustNum
            assert(all(size(arg{ii})==[obj.d.y, obj.d.y]), 'matrix R is not conformable to dim(y)');
            obj.par.R      = arg{ii};
            exhaustNum = true;
        elseif isa(arg{ii}, 'function_handle') && ~exhaustFn
            obj.par.Dh = arg{ii};
            if obj.emiNLhasParams
                f = @(x) obj.par.Dh(x, obj.emiNLParams);
            else
                f = obj.par.Dh;
            end
            internalTestFun(f, 'Dh', true, obj.d.x, [obj.d.y, obj.d.x]);
        elseif isstruct(arg{ii})
            obj.par.emiNLParams = arg{ii};
            obj.emiNLhasParams = true;
        else
            error('Don''t know what to do with argument %d in Emission section', ii);
        end
    end
    if ~exhaustNum
        error('No covariance matrix detected for emission distribution');
    end
end

function obj = internalProcessDat(obj, arg)
    narg        = numel(arg);
    if narg > 2
        error('Too many data arguments. Should only be one: obs or number to generate');
    end
    
    [~, ord]    = sort(cellfun(@numel, arg), 'descend');
    arg         = arg(ord);
    if utils.base.isscalarint(arg{1})
        obj.d.T = arg{1};
        assert(narg == 1, 'Too many data arguments');
        obj   = obj.generateData;
    else
        assert(isnumeric(arg{1}), 'observation data must be a scalar matrix');
        obj.y = arg{1};
        obj.d.T = size(obj.y, 2);
        if obj.opts.warnings && size(obj.y,1) > obj.d.y
            warning('more dimensions in observations than timepoints. y is (d x T) matrix.');
        end
        if obj.opts.warnings && narg == 2
            if ~utils.base.isscalarint(arg{2})
                warning('extraneous argument in data section is unexpected type (%s). Ignoring..', class(arg{2}))
            else
                if arg{2} ~= obj.d.T
                    warning('additional data argument interpreted as dimension: different dim to supplied data y. Ignored.');
                end
            end
        end
    end
    
end

function internalTestFun(f, nm, isDif, inputDim, outputDim)
    try
        fRng = f(ones(inputDim,1));
    catch ME
        warning('Tried %s with input of ones(%d, 1). Output error message:\n', nm, inputDim);
        rethrow(ME);
    end
    assert(all(size(fRng)==outputDim), ['%s does not map into output space ', ...
        '(column vec in R^(%d x %d))'], nm, outputDim(1), outputDim(2));
    
    if ~isDif
        try
            fRng = f(ones(inputDim,10));
        catch ME
            warning('Tried %s with matrix of ones(%d, 10). %s must be columnwise vectorised:\n', nm, inputDim);
            rethrow(ME);
        end
        assert(all(size(fRng)==[outputDim(1),10]), ['%s does not process column-wise ', ...
            'inputs correctly (note D%s need not do so)'], nm, nm);
    end
end