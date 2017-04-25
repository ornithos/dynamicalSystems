function par = getParams(obj, stage, type, fns, nn)
    % par = getParams(obj, stage, type)
    % get parameters from dynamicalSystems object, either for transition
    % or emission stage.
    %
    % INPUTS:
    % obj     - a dynamicalSystems object
    % stage   - (1) = transition; (2) = emission.
    % type    - (0) = linear;     (2) = nonlinear.
    % fns     - struct containing output of obj.functionInterfaces (used to
    %           save time for repeated calling. OPTIONAL.)
    % nn      - series number (if dsBatch object) - will error o.w.
    %
    % OUTPUTS:
    % A struct with the following fields:
    % Q       - output covariance noise.
    %   ----- either one of the following two depending on type ------
    % {A, B}  - plant/emission matrix A; control matrix B if applicable.
    % {f, Df} - plant/emission function f and deriv Df.
    
    par = struct('control', false);
    if nargin < 4 || isempty(fns)
        returnAsStruct = true;
        fns    = obj.functionInterfaces(returnAsStruct);
    end
    givenSeriesNN = nargin == 5 && ~isempty(nn);
    
    if stage == 1
            par.Q = obj.par.Q;
            if type == 0
                par.A  = obj.par.A;
                if obj.hasControl(1) && ~isempty(obj.par.B)
                    par.B = obj.par.B;
                    par.control = true;
                end
                par.bias = zeros(obj.d.x, 1); %obj.par.b;
            else
                par.f  = fns.f;
                par.Df = fns.Df;
            end
            
        elseif stage == 2
            par.Q = obj.par.R; % ok - Q refers to the covariance matrix regardless of stage.
            if type == 0
                
                % EVOLUTION AND CONTROL
                par.A = obj.par.H; % ok - A refers to the transition/emission matrix regardless of stage.
                if obj.hasControl(2)
                    par.B = obj.par.C;   % ok - B refers to the linear control matrix regardless of stage.
                end
                
                % BIAS
                if isempty(obj.par.c)
                    par.bias = zeros(obj.d.y, 1);
                else
                    if givenSeriesNN && iscell(obj.par.c)
                        par.bias = obj.par.c{nn};
                    else
                        par.bias = obj.par.c;
                    end
                end
                
            else
                par.f  = fns.h;
                par.Df = fns.Dh;
            end
        else
            error('Unknown stage requested. Try 1=prediction or 2=update');
    end
end