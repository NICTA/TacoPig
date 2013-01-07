classdef JenkinsReport <handle
    
    properties(Constant)
        header = '<?xml version="1.0" encoding="utf-8"?>\n';
    end
    properties
        name;
        suites;
    end
    methods
        function this = JenkinsReport(name)
            this.name = name;
            this.suites = {};
        end
        function addsuite(this,newsuite)
            this.suites{end+1} = newsuite;
        end
        function outstr = print(this)
            outstr = this.header;
            errs = 0;
            fails = 0;
            skip = 0;
            tests = 0;
            time = 0;
            for i=1:length(this.suites)
                current = this.suites{i};
                errs = errs + current.errs();
                fails = fails + current.fails();
                skip = skip + current.skip();
                tests = tests + current.testcount();
                time = time + current.time();
            end
            suiteshdr = sprintf('<testsuites errors="%d" failures="%d" name="%s" tests="%d" time="%f">\\n',...
                                    errs, fails, this.name, tests, time);
            outstr = strcat(outstr, suiteshdr);
            for i=1:length(this.suites)
               outstr = strcat(outstr, this.suites{i}.print());
            end
            outstr = strcat(outstr, '</testsuites>\n');
        end
    end
    
end