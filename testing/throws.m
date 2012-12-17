function caught = throws(in, task, desired_exception, varargin)
    
        caught = 0;
        try
           eval(task);
        catch e
           caught = strcmp(e.identifier,desired_exception);
           
           if (~caught)
               fprintf('Warning: we got a different error - %s', e.message);
           end
           
        end