classdef dynamicalSystemBatchPars
    % exists soleley for use with dynamicalSystemBatch - allows overloading
    % of parameter struct fields so as to update each dynamicalSystem
    % within a dynamicalSystemBatch collection.
    %
    % REENGINEERED dynamicalSystemBatch ==> THIS IS NOW DEPRECATED!
    %
    
    properties (GetAccess = private)
        dsBatchObj
        A
        Q
        H
        R
    end
    
    methods
        function obj = dynamicalSystemBatchPars(dsBatchHandle)
            obj.dsBatchObj = dsBatchHandle;
        end
        
        function obj = set.A(obj, val)
            s        = struct('type', '.', 'subs', {'dsBatchObj'});
            superObj = builtin('subsref',obj, s);
            dsarr    = superObj.allNonEmptyDS;
            for ii = 1:numel(dsarr)
                dsarr{ii}.par.A = val;
            end
        end
        
        function obj = set.Q(obj, val)
            s        = struct('type', '.', 'subs', {'dsBatchObj'});
            superObj = builtin('subsref',obj, s);
            dsarr    = superObj.allNonEmptyDS;
            for ii = 1:numel(dsarr)
                dsarr{ii}.par.Q = val;
            end
        end
        
        function obj = set.H(obj, val)
            s        = struct('type', '.', 'subs', {'dsBatchObj'});
            superObj = builtin('subsref',obj, s);
            dsarr    = superObj.allNonEmptyDS;
            for ii = 1:numel(dsarr)
                dsarr{ii}.par.H = val;
            end
        end
        
        function obj = set.R(obj, val)
            s        = struct('type', '.', 'subs', {'dsBatchObj'});
            superObj = builtin('subsref',obj, s);
            dsarr    = superObj.allNonEmptyDS;
            for ii = 1:numel(dsarr)
                dsarr{ii}.par.H = val;
            end
        end
            
    end
end