## dynamicalSystems
A MATLAB package for inference and learning in (discrete time) Linear Dynamical System models (includes Gaussian State Space Models or Kalman Filters). The implementation is from the machine learning, or generative perspective rather than the engineering, MLVUE perspective. Many Machine Learning textbooks give a good treatment of this, such as Bishop, Barber and Murphy<sup>1</sup>. My hope is that this is a more friendly implementation of this family of models than some. 

### Features
* Linear Dynamical System (sampling, filtering, smoothing and learning).
* Nonlinear Dynamical System (sampling, filtering, smoothing via EKF or UKF).
* Control inputs may be specified in both transition and/or emission state, as required.
* IONLDS (Georgatzis et al. <sup>2</sup>) learning implemented (requires optimisation toolbox).
* Loglikelihood, expected log joint, and free energy available for debugging at any stage.
* GUI available for comparing filter/smooth/learning in the 2D case.

-------------------------------------------------------------
## Using the package
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
### Specifying the object
The input parser is fairly flexible and thus may be a little intimidating initially. This demonstrates the flexibility of the models permitted by the dynamicalSystems object. It is built to be fairly robust and forgiving: examples are given below. A full specification is:

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

#### Emissionn specification (* required)
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
#### Prior (* required)
```matlab
dynamicalSystem(..., 'x0', m0, P0, ...
```
`x0` is the mean of the prior distribution of the latent space (if not given, this is assumed zero), and `P0` is the covariance matrix. `P0` must be specified, as it can have a considerable impact on the inference. If unknown, use an ormative prior *with respect to the expected scale of the latent space*. For brevity a scalar may be given instead of a full covariance matrix, which will be interpreted to mean a spherical covariance matrix with the given variance in each dimension. This is also true for the evolution and emission parameters.

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
### Methods
Let's assume we have an object
```matlab
ds = dynamicalSystem(....)
```

#### Sampling
Note that in the intial examples, if no data is specified, the constructor automatically generates a sample. This can be requested manually, and existing samples overwritten by the following command:
```matlab
ds.generateData;
```
Samples may also be added to the workspace via `y = ds.generateData;`.

#### Inference

#### Learning

#### Save



###### Footnotes
<sup>1</sup> Chris Bishop - Pattern Recognition and Machine Learning, 2006, Springer, Chapter 13; David Barber - Bayesian Reasoning and Machine Learning, Cambridge University Press, 2012, Chapter 24; Kevin Murphy - Machine Learning: A Probabilistic Perspective, MIT Press, 2012, Chapter 18. All these are available online in pdf format.

<sup>2</sup> Georgatzis et al., 2016 [Input-Output Non-Linear Dynamical Systems applied to Physiological Condition Monitoring](https://arxiv.org/abs/1608.00242)
