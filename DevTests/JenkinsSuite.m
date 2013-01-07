classdef JenkinsSuite <handle
        
    properties
        name;
        tests;
    end
    methods
        function this = JenkinsSuite(name)
            %this = this@handle();
            this.name = name;
            this.tests = {};
        end
        function addtest(this,newtest)
            this.tests{end+1} = newtest;
        end
        function count = errs(this)
            count = 0;
            for i=1:length(this.tests)
                count = count + this.tests{i}.err;
            end
        end
        function count = fails(this)
            count = 0;
            for i=1:length(this.tests)
                count = count + this.tests{i}.fail;
            end
        end
        function count = skip(this)
            count = 0;
            for i=1:length(this.tests)
                count = count + this.tests{i}.skip;
            end
        end
        function count = testcount(this)
            count = length(this.tests);
        end
        function count = time(this)
            count = 0;
            for i=1:length(this.tests)
                count = count + this.tests{i}.time;
            end
        end
        
        function outstr = print(this)
            
            errs = 0;
            fails = 0;
            skip = 0;
            testcount = length(this.tests);
            time = 0;
            for i=1:testcount
                current = this.tests{i};
                errs = errs + current.err;
                fails = fails + current.fail;
                skip = skip + current.skip;
                time = time + current.time;
            end
            outstr = sprintf('    <testsuite errors="%d" failures="%d" name="%s" tests="%d" time="%f">\\n',...
                                    errs, fails, this.name, testcount, time);
            for i=1:testcount
               outstr = strcat(outstr, this.tests{i}.print());
            end
            
            outstr = strcat(outstr, '    </testsuite>\n');
        end
    end
    
end