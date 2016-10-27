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
        ii = 3;
        assert(ischar(args{ii}), 'Unexpected argument in position %d: expected string', ii);
        curSection = internalProcessChar(args{ii}, ii, obj.opts.warnings);
        ii = ii + 1;
        
        while ii <= nargs
            jj = ii;
            while jj < nargs && ~ischar(args{jj+1})
                jj = jj + 1;
            end
            
            if ischar(args{jj})
                error('Insufficient arguments for section (%s)', curSection);
            end
            
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
        end
        
        obj      = internalProcessDat(obj, dataStuff);
        allSects = {'evo', 'emi', 'x0 ', 'dat'};
        if ~all(trackDone)
            error('No arguments for sections: %s', strjoin(allSects(~trackDone), ','));
        end
end

function arg = internalProcessChar(arg, ii, doWarning)
    orig = arg;
    if numel(arg < 3); arg = [arg, '   ']; arg = arg(1:3); end
    assert(any(strcmpi(arg(1:3), {'evo', 'emi', 'x0 ', 'dat'})), ...
        ['string argument #%d (%s) is not one of', ...
        ' the valid markers {x0, evolution, emission, data}'], ii, orig);
    argsFullnames = {'x0', 'evolution', 'emission', 'data'};
    if doWarning && ~any(strcmpi(orig, argsFullnames))
        chosen = find(arg, {'evo', 'emi', 'x0 ', 'dat'});
        warning('arg %d: interpreted %s as %s', orig, argsFullnames{chosen});
    end
end

function obj = internalProcessX0(obj, arg)
    nargs       = numel(arg);
    switch nargs
        case 1
            if all(size(arg{1})==1)
                obj.x0 = struct('mu', zeros(obj.d.x,1), 'sigma', eye(obj.d.x)*arg{1});
            elseif all(size(arg{1})==[obj.d.x, obj.d.x])
                obj.x0 = struct('mu', zeros(obj.d.x,1), 'sigma', arg{1});
            else
                error('Invalid dimensions of x0 prior covariance');
            end
        case 2
            assert(numel(setxor(size(arg{1}), [obj.d.x,1]))==0, 'x0mu incompatible with size of latent space');
            assert(all(size(arg{2})==1) || all(size(arg{2})==[obj.d.x, obj.d.x]), ...
                'x0sigma incompatible with size of latent space');
            obj.x0 = struct('mu', arg{1}, 'sigma', []);
            if all(size(arg{2})==1)
                obj.x0.sigma = eye(obj.d.x)*arg{1};
            else
                obj.x0.sigma = arg{1};
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
        obj.A = arg{1};
        obj.evoLinear = true;
    elseif isa(arg{1}, 'function_handle')
        obj.f = arg{1};
        try
            if ~obj.evoNLhasParams
                fRng = obj.f(ones(obj.d.x,1));
            else
                fRng = obj.f(ones(obj.d.x,1), obj.evoNLParams);
            end
        catch ME
            warning('Tried f with input of ones(%d, 1). Output error message:\n', obj.d.x);
            rethrow(ME);
        end
        assert(all(size(fRng)==[obj.d.x,1]), 'f does not map back into (column vec) in R^%d', obj.d.x);
    else
        error('Unknown argument type (%s) in Evolution argument %d', class(arg{1}), 1);
    end
    
    exhaustNum = false;
    exhaustFn = false;
    for ii = 2:nargs
        if isnumeric(arg{ii}) && ~exhaustNum
            assert(all(size(arg{ii})==[obj.d.x, obj.d.x]), 'matrix Q is not conformable to dim(x)');
            obj.Q      = arg{ii};
            exhaustNum = true;
        elseif isa(arg{ii}, 'function_handle') && ~exhaustFn
            obj.Df = arg{ii};
            try
                if ~obj.evoNLhasParams
                    fRng = obj.Df(ones(obj.d.x,1));
                else
                    fRng = obj.Df(ones(obj.d.x,1), obj.evoNLParams);
                end
            catch ME
                warning('Tried Df with input of ones(%d, 1). Output error message:\n', obj.d.x);
                rethrow(ME);
            end
            assert(all(size(fRng)==[obj.d.x,obj.d.x]), 'Df does not map to R^(%d x %d) required by Hessian', obj.d.x, obj.d.x);
            exhaustFn = true;
        elseif isstruct(arg{ii})
            obj.evoNLParams = arg{ii};
            obj.evoNLhasParams = true;
        else
            error('Don''t know what to do with argument %d in Evolution section', ii);
        end
    end       
end
            
function obj = internalProcessEmi(obj, arg)
    nargs         = numel(arg);
    obj.emiLinear = false;
    %{H || h || []}, {Dh, R},
    if isnumeric(arg{1})
        assert(all(size(arg{1})==[obj.d.y, obj.d.x]), 'matrix H is not conformable to dim(x), dim(y)');
        obj.H = arg{1};
        obj.emiLinear = true;
    elseif isa(arg{1}, 'function_handle')
        obj.h = arg{1};
        try
            if ~obj.emiNLhasParams
                hRng = obj.h(ones(obj.d.x,1));
            else
                hRng = obj.h(ones(obj.d.x,1), obj.emiNLParams);
            end
        catch ME
            warning('Tried h with input of ones(%d, 1). Output error message:\n', obj.d.x);
            rethrow(ME);
        end
        assert(all(size(hRng)==[obj.d.y,1]), 'h does not map into output space (column vec in R^%d)', obj.d.y);
    else
        error('Unknown argument type (%s) in Emission argument %d', class(arg{1}), 1);
    end
    
    exhaustNum = false;
    exhaustFn = false;
    for ii = 2:nargs
        if isnumeric(arg{ii}) && ~exhaustNum
            assert(all(size(arg{ii})==[obj.d.y, obj.d.y]), 'matrix R is not conformable to dim(y)');
            obj.R      = arg{ii};
            exhaustNum = true;
        elseif isa(arg{ii}, 'function_handle') && ~exhaustFn
            obj.Dh = arg{ii};
            try
                if ~obj.emiNLhasParams
                    hRng = obj.Dh(ones(obj.d.x,1));
                else
                    hRng = obj.Dh(ones(obj.d.x,1), obj.emiNLParams);
                end
            catch ME
                warning('Tried Dh with input of ones(%d, 1). Output error message:\n', obj.d.x);
                rethrow(ME);
            end
            assert(all(size(hRng)==[obj.d.y,obj.d.x]), 'Df does not map to R^(%d x %d) required by Hessian', obj.d.y, obj.d.x);
            exhaustFn = true;
        elseif isstruct(arg{ii})
            obj.emiNLParams = arg{ii};
            obj.emiNLhasParams = true;
        else
            error('Don''t know what to do with argument %d in Emission section', ii);
        end
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
        if obj.opts.warnings
            if ~utils.base.isscalarint(arg{2})
                warning('extraneous argument in data section is unexpected type (%s). Ignoring..', class(arg{2}))
            else
                if arg{2} ~= obj.d.T
                    warning('additional data argument interpreted as dimension: different dim to supplied data y');
                end
            end
        end
    end
    
end