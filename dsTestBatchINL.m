% Setup globals
    cols          = zeros(3,3);
    cnums         = [2 4 5];
    for kk = 1:3; cols(kk,:) = utils.plot.varyColor2(cnums(kk)); end
    litecols      = utils.plot.colortint(cols, 0.8);
    cols          = utils.plot.colortint(litecols, 1.4);
    pforward      = 150;

    
doPlot  = true;

% get patient covariates
demogLkp  = tblDemog;
pxLkp     = ionlds.utils.getPatientLookup;
demogLkp.ID = pxLkp(:,2);
demogLkp  = sortrows(demogLkp, 'ID');

% setup spline
%[~,bspl] = ionlds.explore.bsplineNonlin('knots', 6, 'doPlot', false);
knots  = [-1000.0000  -70.6237   -8.3125   -1.5334    2.7607    8.7608   70.7608  1000.0000];
bspl   = utils.spline.bspline(3, knots);
h      = @(z, pars) pars.bspl.functionEval(pars.C * z + pars.bias, pars.eta);

% setup shorthand for different nonlinear optimisation schemes
nlOpts = struct;
nlOpts.CoodAsc   = struct('chgEtas', logical([0 1 1 1 1 1 1 1 1 0]), 'display', 'none', ...
                 'optimType', 'fminunc', 'bfgsSpline', false);
nlOpts.AllConstr = struct('chgEtas', logical([0 1 1 1 1 1 1 1 1 0]), 'display', 'none', ...
                 'optimType', 'fmincon', 'bfgsSpline', true, 'small_iter', 200);
nlOpts.fixEta    = struct('optimType', 'fminunc', 'fixEta', true, 'bfgsSpline', false);
nlOpts.fixC      = struct('optimType', 'fminunc', 'fixC', true, 'bfgsSpline', false);
nlOpts.fixAll    = struct('optimType', 'fminunc', 'fixC', true, 'fixEta', true, 'bfgsSpline', false);

% longest / shortest effect types
% (half life --> base):    2^(-1/n)
% (base --> half-life):    -1/(log2(a))
hlfLifeRng       = [10, 70];
hlfLifeRng       = 2.^(-1./hlfLifeRng);


%%
A          = dsBatch.par.A;
os         = dsBatch.d.y;
ss         = 4;   %dsBatch.d.x;
Q          = dsBatch.par.Q;
R          = dsBatch.par.R;

% initialise nonlinearity
% ---> making it as close to linear as possible. However, now bias has to
% take into account that 0 corresponds to median value  :S
etaUB      = [140, 80, 95]';
% eta        = repmat([0,0,(1 + exp(-0.5*linspace(-5,5,6))).^(-1),1,1],os,1);
eta        = repmat([0,0,(knots(2:end-1)-knots(2))./(diff(knots([2,end-1]))),1,1], os, 1);
eta        = bsxfun(@times, eta, etaUB);
hpar       = struct('eta', eta, 'C', dsBatch.par.H, 'bias', dsBatch.par.c - ones(3,1)*knots(2), 'bspl', bspl);
opts       = struct('warnings', true);

x0mu       = dsBatch.par.x0.mu;

% ____ Get IONLDS object __________________________________________
dsINLBat   = ds.ionldsBatch(ss, os, 'x0', x0mu, eye(ss)*1e-6, 'evolution', A, Q, ...
                'emission', h, hpar, R, 'data', Ys, 'control', Us, ...
                dsBatch.par.B, false, opts);

dsINLBat.smooth('ukf', [], 'forceFilter', true);
dsINLBat.save('initial-inference');

%     dsINL.useSavedParameters('initial-inference', false);  % restart
%%

% ___ Parameter Estimation _________________________________________
nlOpts                = utils.struct.addFieldToSubStructs(nlOpts, 'monotoneWarning', false);
nlOpts                = utils.struct.addFieldToSubStructs(nlOpts, 'etaUB', etaUB);
doNL                  = nlOpts.CoodAsc;

%     doNL.dbgVisualisation = true;
emOpts  = struct('maxiter', 10, 'epsilon', 5e-4, 'multistep', 4, 'ssid', false, ...
          'verbose', true, 'stableVerbose', false, 'annealingSchedule', Inf, 'annealingIter', 100, ...
          'annealingMin', 1e-3, 'diagA', true, 'diagQ', false, 'diagAconstraints', hlfLifeRng, ...
          'strictNegativeCheck', false, 'filterType', 'ukf', 'fixBias2', false, 'nonlinOptimOpts', doNL);

% Learn spline coefficients (with batch pars fixed)

fixAll  = struct('fixA', true, 'fixB', true, 'fixQ', true, 'fixH', true, 'fixR', true);
emOptsFix = utils.struct.structCoalesce(fixAll, emOpts);
% emOptsFix.nonlinOptimOpts.fixC = true;
emOptsFix.maxiter = 1;
% fprintf('Learning spline...\n');
% dsINLBat.parameterLearningEM_bspl(emOptsFix);

% Now learn all emission parameters (with transition fixed).
emOptsFix.fixR = false;
% emOptsFix.nonlinOptimOpts.fixC = false;
fprintf('Learning Emission...\n');
dsINLBat.parameterLearningEM_bspl(emOptsFix);
dsINLBat.save('emission-adapt1');

% Learning begins in anger..
llhhist = dsINLBat.parameterLearningEM_bspl(emOpts);

dsINLBat.smooth('ukf', [], 'forceFilter', true);
dsINLBat.getFittedValues;
dsINLBat.save('em-initial');

ds.plot.predictions(dsINLBat, [1,10,20,Inf], pump.target);


%%
emOptsFix.fixR = true; emOptsFix.nonlinOptimOpts.fixC = true; emOptsFix.nonlinOptimOpts.fixEta = true; emOptsFix.fixB = false;
emOptsFix.maxiter = 5;
llhhist = dsINLBat.parameterLearningEM_bspl(emOptsFix);
    %% Visualise difference between batch parameters, and those with individual biases.
    % The individual bias SMSE is better, but also qualitatively very
    % encouraging for use with nonlinearity. Nevertheless, one imagines
    % emission parameters may be required to be very different.
    figh1 = figure;
    figh2 = figure;
    
for nn = 1:40
    ds.plot.predictions(dsBatchBiasCov, [1,10,20,Inf], pump.target, 'nn', nn, 'titlePrefix', sprintf('Patient %d',nn), 'verbose', false, 'figure', figh1); 
    ds.plot.predictions(dsBatch, [1,10,20,Inf], pump.target, 'nn', nn, 'titlePrefix', sprintf('Patient %d',nn), 'verbose', false, 'figure', figh2);
    pause
end
