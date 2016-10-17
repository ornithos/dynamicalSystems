function obj = validationInference(obj)
    
    assert(~isempty(obj.A), 'No A found in object. Please intialise A');
    assert(~isempty(obj.H), 'No H found in object. Please intialise H');
    assert(~isempty(obj.Q), 'No Q found in object. Please intialise Q');
    assert(~isempty(obj.R), 'No R found in object. Please intialise R');
    assert(obj.isnummat(obj.A) && all(size(obj.A)==[obj.d.x, obj.d.x]), 'A is not a conformable matrix to x');
    assert(obj.isnummat(obj.H) && all(size(obj.H)==[obj.d.y, obj.d.x]), 'H is not a conformable matrix to x and/or y');
    assert(obj.isnummat(obj.Q) && all(size(obj.Q)==[obj.d.x, obj.d.x]), 'Q is not a conformable matrix to x');
    assert(obj.isnummat(obj.R) && all(size(obj.R)==[obj.d.y, obj.d.y]), 'R is not a conformable matrix to y');
    
    
end