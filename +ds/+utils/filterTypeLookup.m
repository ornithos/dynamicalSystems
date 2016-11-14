function out = filterTypeLookup(type, bNumeric)
    %'Kalman', 'EKF', 'Unscented', 'Particle'
    assert(ischar(type), 'specified filter type must be character');
    if nargin < 2 || isempty(bNumeric)
        bNumeric = false;
    end
    assert(isscalar(bNumeric) && islogical(bNumeric), 'bNumeric must be true/false');
    
    type      = [type, '   '];
    shorttype = lower(type(1:3));

    switch lower(shorttype)
        case {'kal', 'lin'}
            out = 1;
        case {'ext', 'ekf', 'fir', 'tay'}
            out = 2;
        case {'uns', 'ukf'}
            out = 3;
        case {'pf ', 'par'}
            out = 4;
        otherwise
            error('Unknown type given. Try ''linear'',''EKF'',or ''UKF''');
    end
    
    if out == 4
        error(['Particle Filters have not yet been implemented. Please use ', ...
            'either the EKF or UKF approximations']);
    end
    
    if ~bNumeric
        types = {'Kalman', 'EKF', 'Unscented', 'Particle'};
        out = types{out};
    end
end