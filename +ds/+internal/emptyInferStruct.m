function s = emptyInferStruct
    s = struct;
    s.filter = struct('mu', [], 'sigma', [], 'utpar', []');
    s.smooth = struct('mu', [], 'sigma', [], 'utpar', [], 'G', [], 'x0', []);
    s.smooth.x0 = struct('mu', [], 'sigma', [], 'G', []);
    s.llh = 0;
    s.fType = '';
    s.sType = '';
    s.fpHash = '';
end