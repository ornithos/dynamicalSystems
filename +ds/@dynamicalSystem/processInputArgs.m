function processInputArgs(obj, args)
% dynamicalSystem(dim_x, dim_y, 'evolution', {A || f || []}, {Df, Q},
   %                 'emission', {H || h || []}, {Dh, R}, c, 
   %                 'data', {y || T},
   %                 'x0', {P0 || p0}, (m0),
   %                ('xtrue', x),
   %                ('control', u, B, (C))
   %                opts)
   % Note that opts is identified purely by whether the last argument is a
   % struct.
   
        nargs = numel(args);
        
        % Copy existing ds object by passing in as only argument.
        if nargs == 1 && isa(args{1}, 'ds.dynamicalSystem')
              % Copy all non-hidden properties.
              p = properties(args{1});
              for i = 1:length(p)
                  obj.(p{i}) = args{1}.(p{i});
              end
              return
        end
        
        if nargs < 4
            error('Insufficient number of inputs. Need dim_x, dim_y, evolution, emission, data and x0');
        end
        
        % update options
        optsDefault     = struct('warnings', true, 'verbose', true);
        if isstruct(args{nargs})
            obj.opts        = utils.base.parse_argumentlist(args{nargs}, optsDefault, true);
            nargs           = nargs - 1;
        else
            obj.opts        = optsDefault;
        end
        
        assert(utils.base.isscalarint(args{1}), 'dimx is not an integer scalar');
        assert(utils.base.isscalarint(args{2}), 'dimy is not an integer scalar');
        obj.d.x = args{1};
        obj.d.y = args{2};
        
        trackDone = false(4,1);
        isControl = false;
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
                    [obj, testsEvo] = internalProcessEvo(obj, args(ii:jj));
                    trackDone(1) = true;
                case 'emi'
                    [obj,testsEmi] = internalProcessEmi(obj, args(ii:jj));
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
                case 'con'
                    controlStuff = args(ii:jj);
                    isControl    = true;
                    % optional: no trackDone
                otherwise
                    error('Unknown section marker ''%s''', curSection);
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
        
        % check we have everything
        allSects = {'evo', 'emi', 'x0 ', 'dat'};
        if ~all(trackDone)
            error('No arguments for sections: %s', strjoin(allSects(~trackDone), ','));
        end
        
        %% This section processes arguments that have been seen already in
        % the main loop above, but for reasons (usually involving the
        % dimensionality of variables!), cannot be processed at the time.
        
        % process data
        [obj, bGen]  = internalProcessDat(obj, dataStuff);
        
        % process control
        if isControl
            obj = internalProcessControl(obj, controlStuff);
        end
        
        % test functions
        for ii = 1:size(testsEvo,1)
            t = testsEvo{ii,1};
            if obj.hasControl(1)
                t(ones(obj.d.u,1), obj.par.evoNLParams);
            else
                t([], obj.par.evoNLParams);
            end
        end
        for ii = 1:size(testsEmi,1)
            t = testsEmi{ii,1};
            if obj.hasControl(2)
                t(ones(obj.d.u,1), obj.par.emiNLParams);
            else
                t([], obj.par.emiNLParams);
            end
        end
        
        % data generation
        if bGen
            obj.generateData;
        end
        
        % do X0
        if ~isempty(xtrue)
            assert(all(size(xtrue) == [obj.d.x, obj.d.T]), 'x must be matrix %d x %d', ...
                    obj.d.x, obj.d.T);
            obj.x = xtrue;
        end
        
        % Check for known problems
        if any(any(obj.y)) && obj.hasControl(2)
            warning('caveat emptor: NaNs not handled yet in EM learning for emission control input. Use EM with EXTREME CAUTION (ie don''t).');
        end
end

function arg = internalProcessChar(arg, ii, doWarning)
    orig = arg;
    if numel(arg < 3); arg = [arg, '   ']; arg = arg(1:3); end
    assert(any(strcmpi(arg(1:3), {'evo', 'emi', 'x0 ', 'dat', 'xtr', 'con'})), ...
        ['string argument #%d (%s) is not one of', ...
        ' the valid markers {evolution, emission, x0, data, xtrue, control}'], ii, orig);
    argsFullnames = {'evolution', 'emission', 'x0', 'data', 'xtrue', 'control'};
    if doWarning && ~any(strcmpi(orig, argsFullnames))
        chosen = strcmpi(arg, {'evo', 'emi', 'x0 ', 'dat', 'xtr'});
        warning('arg %d: interpreted ''%s'' as ''%s''', ii, orig, argsFullnames{chosen});
    end
    arg = lower(arg);
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

function [obj, tests] = internalProcessEvo(obj, arg)
    nargs         = numel(arg);
    tests         = {};
    obj.evoLinear = false;
    % {A || f || []}, {Df, Q}
    if isnumeric(arg{1})
        if isempty(arg{1})
            % unknown transition matrix
            obj.evoLinear = true;
        else
            assert(all(size(arg{1})==[obj.d.x, obj.d.x]), 'matrix A is not conformable to dim(x)');
            obj.par.A = arg{1};
            obj.evoLinear = true;
        end
    elseif isa(arg{1}, 'function_handle')
        obj.par.f = arg{1};
        f = obj.par.f;
        tests = {@(u,pars) internalTestFun(f, 'f', u, false, obj.d.x, [obj.d.x,1],pars), 'evo'};
    else
        error('Unknown argument type (%s) in Evolution argument %d', class(arg{1}), 1);
    end
    
    exhaustNum = false;
    exhaustFn = false;
    for ii = 2:nargs
        if isnumeric(arg{ii}) && ~exhaustNum
            if numel(arg{ii}) == 1
                arg{ii} = arg{ii} * eye(obj.d.x);
            end
            assert(all(size(arg{ii})==[obj.d.x, obj.d.x]), 'matrix Q is not conformable to dim(x)');
            obj.par.Q      = arg{ii};
            exhaustNum = true;
        elseif isa(arg{ii}, 'function_handle') && ~exhaustFn
            obj.par.Df = arg{ii};
            f          = obj.par.Df;
            tests{end+1,1} = @(u,pars) internalTestFun(f, 'Df', u, true, obj.d.x, [obj.d.x, obj.d.x], pars);  %#ok
            tests{end,2}   = 'evo';
            exhaustFn = true;
        elseif isstruct(arg{ii})
            obj.par.evoNLParams = arg{ii};
            obj.evoNLhasParams = true;
        else
            error('Don''t know what to do with argument %d in Evolution section', ii);
        end
    end
    
    if ~exhaustNum
        error('No covariance matrix detected for evolution distribution');
    end
end
            
function [obj,tests] = internalProcessEmi(obj, arg)
    nargs         = numel(arg);
    obj.emiLinear = false;
    tests         = {};
    % {H || h || []}, {Dh, R},
    if isnumeric(arg{1})
        if isempty(arg{1})
            % unknown emission matrix
            obj.emiLinear = true;
        else
            assert(all(size(arg{1})==[obj.d.y, obj.d.x]), 'matrix H is not conformable to dim(x), dim(y)');
            obj.par.H = arg{1};
            obj.emiLinear = true;
        end  
    elseif isa(arg{1}, 'function_handle')
        obj.par.h = arg{1};
        f = obj.par.h;
        tests{end+1,1} = @(u, pars) internalTestFun(f, 'h', u, false, obj.d.x, [obj.d.y,1], pars);
        tests{end, 2}  = 'emi';
    else
        error('Unknown argument type (%s) in Emission argument %d', class(arg{1}), 1);
    end
    % bias?
    if isnumeric(arg{end})
        if isempty(arg{end})
            % unknown emission bias
            obj.par.c = zeros(obj.d.y,1);
            arg(end) = []; nargs = nargs - 1;
        elseif all(size(arg{end})==[obj.d.y, 1])
            obj.par.c = arg{end};
            arg(end) = []; nargs = nargs -1;
        end  
    end
    exhaustNum = false;
    exhaustFn = false;
    for ii = 2:nargs
        if isnumeric(arg{ii}) && ~exhaustNum
            if numel(arg{ii}) == 1
                arg{ii} = arg{ii} * eye(obj.d.y);
            end
            assert(all(size(arg{ii})==[obj.d.y, obj.d.y]), 'matrix R is not conformable to dim(y)');
            obj.par.R      = arg{ii};
            exhaustNum     = true;
        elseif isa(arg{ii}, 'function_handle') && ~exhaustFn
            obj.par.Dh = arg{ii};
            f = obj.par.Dh;
            tests{end+1,1} = @(u, pars) internalTestFun(f, 'Dh', u, true, obj.d.x, [obj.d.y, obj.d.x], pars); %#ok
            tests{end,2}   = 'emi';
        elseif isstruct(arg{ii})
            obj.par.emiNLParams = arg{ii};
            obj.emiNLhasParams = true;
        else
            if isnumeric(arg{ii}) && isempty(obj.d.c)
                warning('Cannot parse argument %d in emission. If it is intended to be bias, ensure it is conformable to y');
            end
            error('Don''t know what to do with argument %d in Emission section', ii);
        end
    end
    if ~exhaustNum
        error('No covariance matrix detected for emission distribution');
    end
end

function [obj, doGenerate] = internalProcessDat(obj, arg)
    narg        = numel(arg);
    doGenerate  = false;
    if narg > 2
        error('Too many data arguments. Should only be one: obs or number to generate');
    end
    
    [~, ord]    = sort(cellfun(@numel, arg), 'descend');
    arg         = arg(ord);
    if utils.base.isscalarint(arg{1})
        obj.d.T = arg{1};
        assert(narg == 1, 'Too many data arguments');
        doGenerate = true;
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

function obj = internalProcessControl(obj, arg)
    nargs       = numel(arg);
    if nargs < 3 || nargs > 3
        error('control must have 3 arguments: (u,B,C)');
    end
    
    sizes = zeros(2,3);
    
    % control inputs
    u = arg{1};
    if size(u,1)==obj.d.T && size(u,2) ~= obj.d.T
        warning('in future please specify ''u'' as (k x T) matrix.');
        u = u'; 
    end
    assert(size(u,2)==obj.d.T, 'u must have same length as emissions');
    obj.d.u = size(u,1);
    obj.u   = u;
    
    for ii = 2:3
        % scalar input
        if (islogical(arg{ii}) || isempty(arg{ii})) && isscalar(arg{ii})
            tmp = arg{ii};
            if isempty(tmp); tmp = false; end
            obj.hasControl(ii-1) = tmp;
            if tmp && ii == 2 
                if obj.evoLinear
                    obj.par.B = zeros(obj.d.x,obj.d.u);
                end
            elseif tmp && ii == 3
                if obj.emiLinear
                    obj.par.C = zeros(obj.d.y,obj.d.u);
                end
            end
            
        % matrix input
        elseif isnumeric(arg{ii}) && ismatrix(arg{ii})
            obj.hasControl(ii-1) = true;
            if ii == 2
                B = arg{ii};
                if size(B,2)==obj.d.x && size(B,1) ~= obj.d.x
                    warning('B is not currently conformable to x, but B transpose is. I''ve corrected.');
                    B = B'; 
                end
                obj.par.B = B;
                assert(all(size(obj.par.B)==[obj.d.x,obj.d.u]), 'B is not conformable to both x and u');
            elseif ii == 3
                C = arg{ii};
                if size(C,2)==obj.d.y && size(C,1) ~= obj.d.y
                    warning('C is not currently conformable to x, but C transpose is. I''ve corrected.');
                    C = C'; 
                end
                obj.par.C = C;
                assert(all(size(obj.par.C)==[obj.d.y,obj.d.u]), 'C is not conformable to both y and u');
            end
        else
            error('control input %d must be numeric matrix or scalar logical', ii);
        end
    end
    
    obj.par.uu        = u * u';
   
end



function internalTestFun(f, nm, u, isDif, inputDim, outputDim, pars)
    if isempty(fieldnames(pars)); haspars = false;
    else, haspars = true; end

    try 
        % switch (hascontrol, haspars)
        if isempty(u)
            if ~haspars
                fRng = f(ones(inputDim,1));
            else
                fRng = f(ones(inputDim,1), pars);
            end
        else
            if ~haspars
                fRng = f(ones(inputDim,1), u);
            else
                fRng = f(ones(inputDim,1), u, pars);
            end
        end
    catch ME
        extramsg = '';
        if strcmp(ME.identifier, 'MATLAB:minrhs') && isempty(u)
            extramsg = ['If the function contains control inputs, ensure the option is switched on in constructor call.'];
        end
        warning('Tried %s with input of ones(%d, 1). %sOutput error message:\n', nm, inputDim, extramsg);
        rethrow(ME);
    end
    assert(all(size(fRng)==outputDim), ['%s does not map into output space ', ...
        '(column vec in R^(%d x %d))'], nm, outputDim(1), outputDim(2));
    
    if ~isDif
        try
            % switch (hascontrol, haspars)
            if isempty(u)
                if ~haspars
                    fRng = f(ones(inputDim,10));
                else
                    fRng = f(ones(inputDim,10), pars);
                end
            else
                if ~haspars
                    fRng = f(ones(inputDim,10), u);
                else
                    fRng = f(ones(inputDim,10), u, pars);
                end
            end
        catch ME
            warning('Tried %s with matrix of ones(%d, 10). %s must be columnwise vectorised:\n', nm, inputDim);
            rethrow(ME);
        end
        assert(all(size(fRng)==[outputDim(1),10]), ['%s does not process column-wise ', ...
            'inputs correctly (note D%s need not do so)'], nm, nm);
    end
end