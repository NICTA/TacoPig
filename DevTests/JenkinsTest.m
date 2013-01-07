classdef JenkinsTest
    
    properties
        name
        err
        fail
        skip
        time
        classname
        message
        description
    end
    
    methods
        
        function this = JenkinsTest(classname,name)
            this.name = name;
            this.err = 0;
            this.fail = 0;
            this.skip = 0;
            this.time = 0;
            this.classname = classname;
            this.message = '';
            this.description = '';
            
        end
        
        function str = print(this)
            str = sprintf('        <testcase classname="%s" name="%s" time="%d"',...
                             this.classname, this.name, this.time);
            if (this.err)
                str = sprintf('%s>\\n            <error message=%s type="MATLAB:assert:failed">%s</failure>\\n        </testcase>\\n', str, this.message, this.description);
            elseif(this.fail)
                str = sprintf('%s>\\n            <failure message=%s type="MATLAB:assert:failed">%s</failure>\\n        </testcase>\\n', str, this.message, this.description);
            elseif (this.skip)
                str = sprintf('%s>\\n            <skipped/>\\n        </testcase>\\n', str);
            else
                str = strcat(str, '/>\n');
            end
        end
        
    end
    
end