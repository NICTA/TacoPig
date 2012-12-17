 function msg = geterr(e)
     str = e.message;
     msg = [];
     active = false; % cuts off the start and the HREF
     spacecount = 0;
     
     
                
     newline = sprintf('\n');
     
     for i=1:length(str)
         
         if (str(i) == '<')
             active = false;
         elseif (str(i) == '>')
             active = true;
         elseif (active)

             if (str(i)==' ')
                 spacecount = spacecount+1;
                 if (spacecount ==2)
                    msg = msg(6:end);
                    msg = sprintf('%-11s ',msg);
                 else
                     msg(end+1) = ' ';
                 end
             else
                 if (str(i) ~= newline)
                    msg(end+1) = str(i); 
                 else
                    msg = [msg, ' :: '];
                 end
             end
         end
     end
     msg(end+1)  = newline;
     