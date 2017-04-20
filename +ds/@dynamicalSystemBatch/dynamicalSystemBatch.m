classdef dynamicalSystemBatch < handle
    % dynamicalSystemBatch: create batch version of dynamicalSystems object
    %
    % dynamicalSystemBatch(dimx, dimy, 'evolution', {A || f || []}, {Df, Q}
    %                 'emission', {H || h || []}, {Dh, R},
    %                 'data', {y || T}, 'x0' {x0mu, x0cov}, 
    %                 ('xtrue', x),
    %                 ('control', u, (C), (D)),
    %                 opts)
    %
    % dynamicalSystemBatch({dsObj1, dsObj2, ...});
    %
    % Each argument is the same as for dynamicalSystem, but is optionally a
    % cell containing multiple copies of time series with assumed same
    % parameters.

    properties
        dsArray = cell(0);
    end
    properties (GetAccess = public)
        par = []
    end

    methods
        function obj = dynamicalSystemBatch(varargin)
            % CONSTRUCTOR
            isCellInp     = cellfun(@iscell, varargin);
            cellLens      = cellfun(@numel, varargin(isCellInp));
            assert(all(cellLens - mean(cellLens) == 0), 'cells must contain same number of elements');
            cellLens      = min([cellLens(:);1]);
            obj.dsArray   = cell(cellLens,1);
            obj.par       = ds.internal.dynamicalSystemBatchPars(obj);
            
            for ii = 1:cellLens
                objInput      = varargin;
                for jj = 1:sum(isCellInp)
                    objInput{jj} = objInput{jj}{ii};
                end
                obj.dsArray{ii} = ds.dynamicalSystem(objInput{:});  % do superclass constructor
            end
        end

        
        %% DYNAMICAL SYSTEM METHODS
        function filter(obj, varargin)
            dsarr  = obj.allNonEmptyDS;
            for ii = 1:numel(dsarr)
                dsarr{ii}.filter(varargin{:});
            end
        end
        
        function smooth(obj, varargin)
            dsarr  = obj.allNonEmptyDS;
            for ii = 1:numel(dsarr)
                dsarr{ii}.smooth(varargin{:});
            end
        end
        
        function ssid(obj, varargin)
            dsarr  = obj.allNonEmptyDS;
            for ii = 1:numel(dsarr)
                dsarr{ii}.ssid(varargin{:});
                ssidrescale       = max(abs(eig(dsarr{ii}.par.A)));
                ssidrescale       = max(ssidrescale, 1e-9);    % ensure stable trans matrix (hack!)
                dsarr{ii}.par.A   = dsarr{ii}.par.A ./ ssidrescale;
                dsarr{ii}.par.H   = dsarr{ii}.par.H .* ssidrescale;
            end
        end
        
        function s = suffStats(obj, varargin)
            
            dsarr  = obj.allNonEmptyDS;
            narr   = numel(dsarr);
            assert(narr > 0, 'No dynamicalSystems in this batch!');
            s      = dsarr{1}.suffStats(varargin{:});
            fnms   = fieldnames(s);
            
            for ii = 2:narr
                stmp = dsarr{ii}.suffStats(varargin{:});
                for kk = 1:numel(fnms)
                    s.(fnms{kk}) = s.(fnms{kk}) + stmp.(fnms{kk});
                end
            end
        end
        
        function llh = logLikelihood(obj, txt, useExisting)  %#ok
            if nargin < 3; useExisting = false; end
            llh = 0;
            dsarr  = obj.allNonEmptyDS;
            for ii = 1:numel(dsarr)
                if useExisting
                    llh = llh + dsarr{ii}.infer.llh;
                else
                    llh = llh + dsarr{ii}.logLikelihood;
                end
            end
        end
        
        function llh = inferLlh(obj)
            llh = 0;
            dsarr  = obj.allNonEmptyDS;
            for ii = 1:numel(dsarr)
                llh = llh + dsarr{ii}.infer.llh;
            end
        end
        
        % prototype
        [llh, niters] = parameterLearningEM(obj, opts);
        % ***********************************************************
        % --> NEED TO IMPLEMENT PARAMETERLEARNINGMSTEP (INCL x0!!) FOR BTCH
        % ***********************************************************
      
        
        %% MISC FUNCTIONS
        
        function out = n(obj)
            out = numel(obj.dsArray);
        end

        function idx = firstNonEmptyIdx(obj)
            idx = find(~cellfun(@isempty, obj.dsArray), 1, 'first');
        end

        function arr = allNonEmptyDS(obj)
            idx = ~cellfun(@isempty, obj.dsArray);
            arr = obj.dsArray(idx);
        end
        
        function garbageCollect(obj)
            elEmpty      = cellfun(@isempty, obj.dsArray);
            obj.dsArray(elEmpty) = [];
        end

        function ok  = validationInference(obj, varargin)
            dsarr  = obj.allNonEmptyDS;
            ok     = true;
            for ii = 1:numel(dsarr)
                ok = ok && dsarr{ii}.validationInference(varargin{:});
            end
        end
        
        % Make a copy of a handle object.
        function new = copy(obj)
            % Instantiate new object of the same class.
            jj           = obj.firstNonEmptyIdx;
            new          = ds.dynamicalSystemBatch(obj.dsArray{jj}.copy);
            if jj ~= 1
                new.dsArray = [repmat({[]}, jj-1, 1); new.dsArray];
            end
            for ii = jj+1:obj.n
                cur                = obj.dsArray{ii};
                if isempty(cur)
                    new.dsArray{ii} = [];
                else
                    new.dsArray{ii} = obj.dsArray{ix}.copy;
                end
            end
        end


        % index access to object
        function sref = subsref(obj,s)
            switch s(1).type
                case '.'
                    noOutput = {'filter', 'smooth', 'ssid', 'garbageCollect'};
                    if ismember(s(1).subs, noOutput)
                        builtin('subsref', obj, s);
                    else
                        sref = builtin('subsref',obj,s);
                    end
%                 case '()'
%                     sref = builtin('subsref',obj,s);
                case '{}'
                    % in case of multiple subscripts > pass remaining to ds
                    ll = length(s);

                    curS = s(1);
                    curS.type = '{}';
                    assert(numel(curS.subs) == 1, 'dynamicalSystemBatch object is 1D');

                    sref = builtin('subsref', obj.dsArray, curS);

                    if ll > 1
                        curS = s(2:end);
                        sref = builtin('subsref', sref, curS);                     
                    end

                otherwise
                    error('bad suscript type: %s (-- I only accept {})', s(1).type);
            end
        end

    end
end