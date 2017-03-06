function [val, D] = logKstep(obj, fType, utpar, k, opts)
    if nargin < 5 || isempty(opts); opts = struct; end
    if nargin < 4 || isempty(k); k = 1; end
    if nargin < 3 || isempty(utpar); utpar = struct; end
    if nargin < 2 || isempty(fType); fType = obj.infer.sType; end
    
    assert(isa(obj, 'ds.dynamicalSystem'), 'first argument must be a dynamical systems object');
    assert(k==1, 'only implemented k = 1');
    
    optsDefaults = struct('T', 1);
    opts         = utils.struct.structCoalesce(opts, optsDefaults);
    
    optsFilter   = struct('bDoValidation', true, 'bIgnoreHash', false, ...
                    'bCollectGradient', true, 'T', opts.T);

    tmpobj       = obj.copy;
    D            = tmpobj.filter(fType, true, utpar, optsFilter);
    val          = tmpobj.infer.llh;
end