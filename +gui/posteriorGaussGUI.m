function posteriorGaussGUI(dsObject, save1, save2)
    % posteriorGaussGUI(dsObject, save1, save2)
    % compare (2D) posteriors between Ground truth, Filter and Smoother for
    % 2 different saved results
    %
    % INPUTS:
    % dsObject - a dynamicalSystems object
    % save1    - a save point. If numeric, the relevant index is retrieved
    %            from the object stack. If character, the save point with
    %            matching description is retrieved. If empty, the current
    %            parameters are used.
    % save2    - the comparison point (acceptable values as above).
    
    assert(nargin == 3, 'Valid call is posteriorGaussGUI(dsObject, save1, save2)');
    addpath(utils.system.matlabPath('dynamicalSystems','guide'));
    
    [sp1, descr1] = getSavePoint(dsObject, save1);
    [sp2, descr2] = getSavePoint(dsObject, save2);
    sp1.descr = descr1;
    sp2.descr = descr2;
    guidePosteriorGaussGUI(dsObject, sp1, sp2);
end

        
function [out, descr] = getSavePoint(dso, sp)
    if isempty(sp)
        try
            dso      = dso.save(['current-', char(floor(75*rand(1, 4)) + 48)]);
        catch ME
            fprintf('Error trying to retrieve current status as save point:\n');
            rethrow(ME);
        end
        [out, descr]     = dso.stackTop;
    else
        sp       = dso.stackFind(sp);
        out      = dso.stack{sp, 1};
        descr    = dso.stack{sp, 2};
    end
end