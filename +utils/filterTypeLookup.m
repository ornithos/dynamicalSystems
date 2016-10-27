function out = filterTypeLookup(type)
    %'Kalman', 'EKF', 'Unscented', 'Particle'
    assert(ischar(type), 'specified filter type must be character');
    type      = [type, '   '];
    shorttype = lower(type(1:3));

    switch lower(shorttype)
        case {'kal', 'lin'}
            out = 'Kalman';
        case {'ext', 'ekf', 'fir', 'tay'}
            out = 'EKF';
        case {'uns', 'ukf'}
            out = 'Unscented';
        case {'pf ', 'par'}
            out = 'Particle';
        otherwise
            error('Unknown filter type specified'}
    end
end