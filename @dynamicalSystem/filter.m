function obj = filter(obj, type, varargin)
    anyNonLinear = ~(obj.evoLinear && obj.emiLinear);
    
    if anyNonLinear 
        if nargin < 2
            type = 'Extended';
            if obj.opts.warning
                warning('Non-linear relations present. In future, please specify the filter required. EKF has been selected by default.');
            end
        end
    else
        if isempty(type)
            type = 'Kalman';
        end
    end
    
    type = utils.filterTypeLookup(type);
    if ~anyNonLinear && type ~= 'Kalman' && obj.opts.warning
        warning('Filter type ignored since exact inference possible');
        type = 'Kalman';
    end
   
    % check validity of specification
    isValidType(type, anyNonLinear);
    
    % Signpost
    switch lower(type)
        case 'kalman'
            obj = obj.filterKalman(varargin{:});
        case 'ekf'
            obj = obj.filterExtended(varargin{:});
        case 'unscented'
            obj = obj.filterUnscented(varargin{:});
        case 'particle'
            obj = obj.filterParticle(varargin{:});
        otherwise
            error('Unknown filter type encountered: (%s)', type);
    end
     
end

    
function isValidType(type, anyNonLinear)
    nlTypes = {'EKF', 'Unscented', 'Particle'};
    nlTypesStr = strjoin(nlTypes, ', ');
    if anyNonLinear
        if ~ismember(type, nlTypes)
            error('Invalid non-linear filter type specified. Should be one of (%s)', nlTypesStr);
        end
    else
        if ~ismember(type, {'Kalman'})
            error('Invalid linear filter type specified. Should be Kalman./');
        end
    end
end