function ok = validationInference(obj, doError)
    
    if nargin < 2
        doError = true;
    end
    ok = true(8,1);
    
    if obj.evoLinear
        ok(1) = testCondition(~isempty(obj.A), 'No A found in object. Please intialise A', doError);
    else
        ok(1) = testCondition(~isempty(obj.f), 'No transition function (f) found. Cannot perform inference', doError);
        if obj.evoNLhasParams
            ok(2) = testCondition(~isempty(obj.evoNLParams), 'Parameter object for transition function (f) not found.', doError);
        end
    end
    
    if obj.emiLinear
        ok(3) = testCondition(~isempty(obj.H), 'No H found in object. Please intialise H', doError);
    else
        ok(3) = testCondition(~isempty(obj.h), 'No emission function (h) found. Cannot perform inference', doError);
        if obj.emiNLhasParams
            ok(4) = testCondition(~isempty(obj.emiNLParams), 'Parameter object for emission function (h) not found.', doError);
        end
    end
    
    ok(5) = testCondition(~isempty(obj.Q), 'No Q found in object. Please intialise Q', doError);
    ok(6) = testCondition(~isempty(obj.R), 'No R found in object. Please intialise R', doError);
%     assert(obj.isnummat(obj.A) && all(size(obj.A)==[obj.d.x, obj.d.x]), 'A is not a conformable matrix to x');
%     assert(obj.isnummat(obj.H) && all(size(obj.H)==[obj.d.y, obj.d.x]), 'H is not a conformable matrix to x and/or y');
%     assert(obj.isnummat(obj.Q) && all(size(obj.Q)==[obj.d.x, obj.d.x]), 'Q is not a conformable matrix to x');
%     assert(obj.isnummat(obj.R) && all(size(obj.R)==[obj.d.y, obj.d.y]), 'R is not a conformable matrix to y');
    ok = all(ok);
end

function cond = testCondition(cond, errorMsg, doError)
    if doError && ~cond
        error(errorMsg);
    end 
end