## dynamicalSystems
A MATLAB package for inference and learning in (discrete time) Linear Dynamical System models (includes Gaussian State Space Models or Kalman Filters). The implementation is from the machine learning, or generative perspective rather than the engineering, LS/MVUE perspective. Many Machine Learning textbooks give a good treatment of this, such as Bishop, Barber and Murphy<sup>1</sup>. My hope is that this is a more friendly implementation of this family of models than some. 

### Features
* Linear Dynamical System (sampling, filtering, smoothing and learning).
* Nonlinear Dynamical System (sampling, filtering, smoothing via EKF or UKF).
* Control inputs may be specified in both transition and/or emission state, as required.
* IONLDS (Georgatzis et al. <sup>2</sup>) learning implemented (requires optimisation toolbox).
* Loglikelihood, expected log joint, and free energy available for debugging at any stage.
* GUI available for comparing filter/smooth/learning in the 2D case.
* Missing values -- both partial and full -- are handled in Linear case (Inference and Learning).
* (NOT YET AVAILABLE) Missing values not fully derived and engineered in the nonlinear case.
* (NOT YET AVAILABLE) Multiple time series can be learned simultaneously.
* (NOT YET AVAILABLE) Particle Filter.

#### Contents:
* [Examples](#examples)
* [Creating a model](#creating-a-model)
* [Sampling](#sampling)
* [Inference](#inference)
* [Learning](#learning)
* [Saving](#saving)
* [Other Methods](#other-methods)

-------------------------------------------------------------
## Using the package
No installation is necessary. The package can be used by cloning the contents into `/your-folder-here/dynamicalSystems` and calling
```matlab
addpath /your-folder-here/dynamicalSystems
```
It can be added permanently to the path by adding the statement to your `startup.m` file. 

###### Dependencies
`dynamicalSystems` depends on certain files in the `matlabUtils` package: this must also be downloaded and the top level folder added to path in the same way as above. Currently the MATLAB optimisation toolbox is required for learning in the IONLDS model. However, this can easily be replaced by a free toolbox if required. (The [minFunc](https://www.cs.ubc.ca/~schmidtm/Software/minFunc.html) project is an excellent choice.)

### Examples

#### Circular motion with noise
```matlab
dsCircle = ds.dynamicalSystem(2, 2, 'x0', 1e-5, 'evolution', [cos(0.03) -sin(0.03); sin(0.03) cos(0.03)], eye(2), ...
           'emission', 0.5*eye(2), 0.2*eye(2), 'data', 50, struct('verbose', false));   % sample from circular path
dsCircle.smooth;                  % perform RTS smoothing
dsCircle.save('original');        % save smoothed inference
dsCircle.save('original-copy');   %  "     "          "
ds.gui.posteriorGaussGUI(dsCircle, 'original', 'original-copy');   % hack to use the GUI for easy viewing
```

#### Newton dynamics
```matlab
% Newtonian dynamics (4D state space to keep linear system with acc/veloc)
deltaT              = 0.5;
transition          = eye(4);
transition(1,3)     = deltaT;
transition(2,4)     = deltaT;

emission            = [1 0 0 0; 0 1 0 0];
Q                   = 0.01*eye(4);
Q(3,3)              = 0.2;
Q(4,4)              = 0.2;
dsNewton  = ds.dynamicalSystem(4, 2, 'x0', 1e-5, 'evolution', transition, Q, 'emission', emission, 5*eye(2), 'data', 100);
dsNewton.smooth;
dsNewton.save('original-params');
```

#### Regression with drift...
Can't do this yet since `A` matrix fixed. If this is allowed to vary and equal the inputs, then we're done.


-------------------------------------------------------------
### Creating a model
The input parser is fairly flexible and may be a little intimidating initially due to the number of possible options. My hope is that it is relatively simple in actuality, and examples are given below. The syntax is built to be fairly robust and forgiving. The full specification is:

```matlab
myDsObject = dynamicalSystem(dimx, dimy, 'evolution', {A || f || []}, (Df), (Q), ...
              'emission', {H || h || []}, (Dh), (R), ...
              'data', {y || T}, 'x0' (x0mu), (x0cov), ...
              ('xtrue', x), ...
              ('control', u, (C), (D)), ...
              opts)
```
Let's break this down. `dimx` and `dimy` are the dimensions of the latent and output states respectively, and must be specified as the first two arguments. All following arguments may be specified in any order, but we will take them in the order listed. 

#### Evolution/transition specification (* required)
This is the (hidden) transition from time step `t` to `t+1`. Let's assume first that we have a linear transition:

```matlab
dynamicalSystem(..., 'evolution', A, Q, ...
```
would specify the known linear transition matrix `A` and the covariance matrix `Q`. If either/both of these are not known, they may be unspecified by passing in an empty vector as required:

```matlab
dynamicalSystem(..., 'evolution', [], [], ...
```

If a nonlinear function is used in the transition instead, pass in the corresponding function `f` instead of the transition matrix. Further, if using the EKF algorithm, the Jacobian must also be specified (`Df`), which may be passed in second, or ignored if using the UKF algorithm. The function(s) given will be tested for conformability as part of the input checks. The covariance matrix `Q` may be passed in if known, or left empty otherwise. If all these parameters were known for instance, the call would be:

```matlab
dynamicalSystem(..., 'evolution', f, Df, Q, ...
```

#### Emission specification (* required)
This is entirely analogous to the previous specfication: e.g.
```matlab
dynamicalSystem(..., 'emission', H, R, ...
```
```matlab
dynamicalSystem(..., 'emission', h, Dh, R, ...
```

#### Data (* required)
There are two expected possibilities when creating a dynamicalSystem object. Either the user has data they wish to model, or they have a given model they wish to sample from. In the first case, the observed data may be given as follows:
```matlab
dynamicalSystem(..., 'data', y, ...
```
where `y` is a `d x T` matrix - each observation is a column, meaning `T` observations corresponds to `T` columns. Otherwise the user is wishing to sample from the model (note that in this case the evolution and emission parameters must be specified), and the argument will simply be the number of time points to sample:
```matlab
dynamicalSystem(..., 'data', T, ...
```
Missing values may be specified as `NaN` where relevant (note the nonlinear case has not been worked through properly -- caveat emptor!).

#### Prior (* required)
```matlab
dynamicalSystem(..., 'x0', m0, P0, ...
```
`x0` is the mean of the prior distribution of the latent space (if not given, this is assumed zero), and `P0` is the covariance matrix. `P0` must be specified, as it can have a considerable impact on the inference. If unknown, use a relatively uninformative prior *with respect to the expected scale of the latent space*. For brevity a scalar may be given instead of a full covariance matrix, which will be interpreted to mean a spherical covariance matrix with the given variance in each dimension. This is also true for the evolution and emission parameters.

#### Actual latent state
```matlab
dynamicalSystem(..., 'xtrue', x, ...
```
This is an optional input to specify the actual latent state if known. This will typically be if the user has already sampled the data and wishes to compare models. It can be used to sample different observation sequences from the same hidden state, and also for visualising the inference performed vs. the actual latent state.

#### Control inputs
```matlab
dynamicalSystem(..., 'control', u, C, D, ...
```
If a dynamical system has control inputs, they may be specified in the vector `u` (see graphical model). The matrix `C` coresponds to the relationship of the control inputs to the latent space, and the matrix `D` gives the relationship to the emission. Frequently only one of these is used, but both may be specified. If the matrices are not known in advance, a boolean value must be given for each of these, so the parser knows which relationship exists. For example,
```matlab
dynamicalSystem(..., 'control', u, true, false, ...
```
is interpreted to mean that the latent state is dependent on the control input, but the emission is not (directly) dependent. This may of course be reversed.

#### Other options
```matlab
dynamicalSystem(..., 'opts', s, ...
```
where `s` is a struct containing a subset of the fields `{'warnings', 'verbose'}`. These are scalar boolean fields that dictate the verbosity of the dynamicalSystem object. No other options are currently accepted.


-------------------------------------------------------------
## Methods
Let's assume we have an object
```matlab
myDsObject = dynamicalSystem(....)
```

### Sampling
Note that in the intial examples, if no data is specified, the constructor automatically generates a sample. This can be requested manually, and existing samples overwritten by the following command:
```matlab
myDsObject.generateData;
```
Samples may also be added to the workspace via `y = ds.generateData;`.

### Inference
Inference may be performed on an existing dynamicalSystems object, say `myDsObject` using the following commands.
```matlab
myDsObject.filter;
myDsObject.smooth;
```
If you're unfamiliar with the concepts of filtering and smoothing, please consult one of the references given at the beginning. These methods take a couple of (optional) arguments.
```matlab
myDsObject.filter(bDoLlh, fType, opts, utpar);
```
It is easiest to calculate the log likelihood during the filtering step, so the first argument `bDoLlh` is a boolean scalar which determines if this is collected during the filtering step. Default is false, since there is a little computational overhead in doing so. Next the filter type (or `ftype`) can be given - this is currently one of `{'Kalman', 'EKF', 'Unscented'}` (various shorthands are accepted) which should be given if the system is nonlinear. `'Unscented'` is typically preferred to `'EKF'` for these problems, but it is useful to compare the results. The particle filter is unfortunately not yet implemented. `'opts'` are really designed for internal use, and `utpar` is the struct of  unscented transform parameters if performing UKF. The struct contains the fields `alpha`, `beta`, `kappa` which default to 1, 0, 0 respectively. The smoothing arguments are the same, except without the first (`bDoLlh` ).

The results of filtering or smoothing can be found in:
```matlab
myDsObject.infer.filter.mu; myDsObject.infer.filter.sigma;
myDsObject.infer.smooth.mu; myDsObject.infer.smooth.sigma;
```
and if collected, the log likelihood can be found in `myDsObject.infer.llh`. If this wasn't collected, a call can be made to the method
```matlab
myDsObject.logLikelihood;
```
which returns the quantity directly.

### Learning
The well known subspace identification algorithm, or spectral method for estimating the model parameters is available via the call:
```matlab
myDsObject.ssid(L)
```
where `L` is the size of the Hankel matrix. This needs to be at least as great as the dimension of the latent state. Larger values may give better results, but must be traded off against the length of the chain.

While the spectral method is known to give a good estimate, it is only justified in the limit of infinite data, and a local search will typically perform better. Therefore it is usually used to initialise the EM algorithm for parameter estimation. To this end, I recommend the following code instead:
```matlab
opts           = struct('maxiter', 1000, 'epsilon', 1e-3, 'ssid', true, 'sampleStability', 10);
history        = myDsObject.parameterLearningEM(opts);
```
This performs the EM algorithm for Maximum Likelihood Estimation, initialised from the subspace identification method (note `'ssid', true`). A maximum of 1000 iterations is specified with a convergence tolerance (in the log likelihood) of 10^(-3). The `sampleStability` option refers to the stable learning of the transition matrix in (Siddiqi et al. 2007<sup>3</sup>) to ensure the latent state doesn't blow up exponentially. The number specified is how frequently the algorithm will perform the check - it will increase sampling time automatically if the spectral radius starts knocking up against the upper limit. If unsure, I recommend leaving this as-is. Other options include varying levels of verbosity, diagonal covariance matrices, and deterministic annealing. See the `optsDefault` structure in the `parameterLearningEM` code for more details.

### Saving
It is common that several different types of inference or learning may be performed. This may be comparing EKF vs UKF, changing the unscented parameters, changing the learning parameters, or trying different initialisations. In all these cases, the current parameters and latent distributions may be saved and recalled as in the following example:

```matlab
myDsObject.filter;
myDsObject.smooth;
myDsObject.ssid(5)
myDsObject.save('initial-params');

opts           = struct('maxiter', 1000, 'epsilon', 1e-3, 'sampleStability', 10);
myDsObject.parameterLearningEM(opts);
myDsObject.save('EM-parameters1');

myDsObject.useSavedParameters('initial-params');
```
which recalls the original estimation of the parameters once the EM learning is complete. The different parameters learned can then be compared (use `getSaved` to return the relevant structs directly). The full list of save points can be viewed using `myDsObject.savedList`. 

### Other methods

```matlab
myDsObject.copy;
```
The dynamicalSystems object is a handle object, so assigning the object to a second variable will not create a copy. Use the `.copy` method instead.

```matlab
yhat = myDsObject.getFittedValues;
```
Return the fitted values for the given time series (ie. the smoothed latent mean transformed to the observation space.

```matlab
yhat = myDsObject.getPredictedValues(nlookahead);
```
As above, but now gives the predicted value of `nlookahead` steps ahead for each time point. Special cases include `0` returning the filtered observation estimates and `-Inf` returning the smoothed fitted values.


```matlab
ds.gui.posteriorGaussGUI(myDsObject, 'initial-params', 'EM-parameters1');
```
Visualise the filtered and smoothed estimates at the specified saved points against each other. If the observation is higher than 2 dimensional, the gui will plot only the first 2 dimensions.

-------------------------------------------------------------
###### Footnotes
<sup>1</sup> Chris Bishop - Pattern Recognition and Machine Learning, 2006, Springer, Chapter 13; David Barber - Bayesian Reasoning and Machine Learning, Cambridge University Press, 2012, Chapter 24; Kevin Murphy - Machine Learning: A Probabilistic Perspective, MIT Press, 2012, Chapter 18. All these are available online in pdf format.

<sup>2</sup> Georgatzis et al., 2016 [Input-Output Non-Linear Dynamical Systems applied to Physiological Condition Monitoring](https://arxiv.org/abs/1608.00242)

<sup>3</sup> Siddiqi et al., 2007, NIPS, A Constraint Generation Approach to Learning Stable Linear Dynamical Systems.
