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
        if nargs == 1 && isa(args{1}, 'ds.dynamicalSystemBatch')
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
                    x0stuff      = args(ii:jj);
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
        
        % process x0
        obj = internalProcessX0(obj, x0stuff);
        
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
        
        % do X (if known)
        if ~isempty(xtrue)
            
            if isnumeric(xtrue); xtrue = {xtrue}; end
            assert(iscell(xtrue), 'xtrue must be a cell array');
            assert(obj.d.n == numel(xtrue), '%d xtrue series given for %d observation series', obj.d.n);
            cellDimChecks(curMu, obj.d.n, obj.d.x, obj.d.T, 'xtrue')  %inp, nElements, nrows, ncols, name
            
            obj.x = xtrue;
        end
        
        % Check for known problems
        if any(cellfun(@(x) any(any(isnan(x))), obj.y)) && obj.hasControl(2)
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
    nargs            = numel(arg);
    obj.par.x0       = struct('mu',[],'sigma',[]);
    
    assert(nargs > 0 && nargs <= 2, 'x0 must have either 1 or 2 arguments');
    if nargs == 1
        obj.par.x0.mu = repmat({zeros(obj.d.x,1)}, obj.d.n, 1);
        fnum = 2;
    else
        fnum = 1;
    end
    
    if fnum == 1
        curMu = arg{1};
        % convert numeric --> cell if given
        if isnumeric(curMu)
            assert(ismatrix(curMu), 'x0 mu must be a matrix or cell of matrices');
            if numel(curMu) == obj.d.x
                curMu = repmat({curMu}, obj.d.n,1);
            elseif all(size(curMu) == [obj.d.x, obj.d.n])
                curMu = mat2cell(curMu, obj.d.x, ones(1,obj.d.n))';
            else
                error('Do not know what to do with x0 mu. Expecting cell of %d x0 means, or single x0 mu', obj.d.n);
            end
        end
        assert(iscell(curMu), 'x0 mu must be a matrix or cell of matrices');
        cellDimChecks(curMu, obj.d.n, obj.d.x, 1, 'x0 mu')  %inp, nElements, nrows, ncols, name
        obj.par.x0.mu = curMu;
        fnum = fnum + 1;
    end
    
    if fnum == 2
        curSigma = arg{2};
        % convert numeric --> cell if given
        if isnumeric(curSigma)
            assert(ismatrix(curSigma), 'x0 sigma must be a matrix or cell of matrices');
            if numel(curSigma) == 1
                curSigma = repmat({curSigma*eye(obj.d.x)}, obj.d.n,1);
            elseif all(size(curSigma) == repmat(obj.d.x,2,1))
                curSigma = repmat({curSigma}, obj.d.n,1);
            else
                error('Do not know what to do with x0 sigma. Expecting cell of %d x0 sigma, or single x0 sigma', obj.d.n);
            end
        end
        assert(iscell(curSigma), 'x0 sigma must be a matrix or cell of matrices');
        cellDimChecks(curSigma, obj.d.n, obj.d.x, obj.d.x, 'x0 sigma')  %inp, nElements, nrows, ncols, name
        obj.par.x0.sigma = curSigma;
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

    if min(size(arg{1})) == 1 && all(all(arrayfun(@utils.base.isscalarint, arg{1})))
        if numel(arg{1}) > 1 && obj.opts.warnings
            warning(['Interpreting ''data'' vector as list of lengths for each dynamicalSystem', ...
                ' not observations. To force as observations, please place in cell array {[...]}.']);
        end
        obj.d.n = numel(arg{1});
        obj.d.T = arg{1};
        doGenerate = true;
    else
        if isnumeric(arg{1})
            assert(ismatrix(arg{1}), 'Tensors not supported in ''data''. If multiple time series, use cell.');
            obj.y = {arg{1}};
            obj.d.T = size(arg{1}, 2);
            obj.d.n = 1;
            if obj.opts.warnings && size(obj.y,1) > obj.d.y
                warning('more dimensions in observations than timepoints. y is (d x T) matrix.');
            end
        else
            assert(iscell(arg{1}), 'observation data must be a cell of matrices');
            obj.d.n   = numel(arg{1});
            filledNaN = false;
            for ii = 1:obj.d.n
                cY     = arg{1}{ii};
                assert(isnumeric(cY) && ismatrix(cY), '''data'' field: cell element %d not a numeric matrix.', ii);
                assert(size(cY,1) <= obj.d.y, '''data'' field: cell element %d of greater dimension (%d) than specified (%d)', ii, size(cY,1), obj.d.y);
                if obj.d.y > size(cY,1)
                    filledNaN = true;
                    arg{1}{ii} = [arg{1}{ii}; NaN(obj.d.y - size(cY,1), size(cY,2))];
                end
            end
            if filledNaN && obj.opts.warnings
                warning('Length of some data series is < dimy. Additional dimensions filled with NaN');
            end
            obj.y   = arg{1};
            obj.d.T = cellfun(@(x) size(x,2), obj.y);
        end
    end
    
end

function obj = internalProcessControl(obj, arg)
    nargs       = numel(arg);
    assert(nargs == 3, 'control must have 3 arguments: (u,B,C)');
    
    % control inputs
    u = arg{1};
    if isnumeric(u); u = {u}; end
    assert(iscell(u), 'control series u must be a cell array');
    assert(obj.d.n == numel(u), '%d control series given for %d observation series', obj.d.n);
    
    obj.d.u = size(u{1},1);
    
    % initial test sizes -- user may have supplied U in wrong orientation.
    assert(size(u{1},2) == obj.d.T(1), 'First control series is of length %d; expecting %d. Bear in mind time is the 2nd dimension.', size(u{1},2), obj.d.T(1))
    uSizes  = cell2mat(cellfun(@(x) size(x), u, 'Un', 0));
    assert(all(uSizes(:,2)==obj.d.T), 'control series %s different length to observations', utils.base.numjoin(find(uSizes(:,2)~=obj.d.T), ',', 0)); %#ok
    assert(all(uSizes(:,1)==uSizes(1,1)), 'control series are not of uniform dimension.');
    
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

function cellDimChecks(inp, nElements, eachNRows, eachNCols, inpType)
    % dimension checks
    assert(numel(inp) == nElements', '%d %s series given: expecting %d.', numel(inp), inpType, nElements);
    xSizes  = cell2mat(cellfun(@(x) size(x), inp, 'Un', 0));
    assert(all(xSizes(:,2)==eachNCols), '%s series: %s different length to observations', ...
            inpType, utils.base.numjoin(find(xSizes(:,2)~=eachNCols), ',', 0)); %#ok
    assert(all(xSizes(:,1)==eachNRows), '%s series are not all of dimension dimx (%s)', ...
            inpType, utils.base.numjoin(find(xSizes(:,2)~=eachNRows), ',', 0));  %#ok
end