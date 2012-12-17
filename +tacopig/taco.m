% All the tacopig classes inherit from the handle class so that they can be
% passed by reference. However, the neccessary mechanisms for this should not 
% appear in the documentation 
classdef taco < handle

    % We will wrap all the methods of handle but make them hidden
     methods (Hidden)

         function varargout = addlistener(varargin)
             % Weird matlab case for passing variable output arguments
             [varargout{1:nargout}] = addlistener@handle(varargin{:});
         end

         function varargout = delete(varargin)      
             [varargout{1:nargout}] = delete@handle(varargin{:});
         end

         function varargout = eq(varargin)           
             [varargout{1:nargout}] = eq@handle(varargin{:});
         end

         function varargout = findobj(varargin)      
             [varargout{1:nargout}] = findobj@handle(varargin{:});      
         end

         function varargout = findprop(varargin)     
             [varargout{1:nargout}] = findprop@handle(varargin{:});
         end

         function varargout = ge(varargin)           
             [varargout{1:nargout}] = ge@handle(varargin{:});
         end

         function varargout = gt(varargin)           
             [varargout{1:nargout}] = gt@handle(varargin{:});
         end

         function varargout = le(varargin)           
             [varargout{1:nargout}] = le@handle(varargin{:});
         end

         function varargout = lt(varargin)           
             [varargout{1:nargout}] = lt@handle(varargin{:});
         end

         function varargout = ne(varargin)           
             [varargout{1:nargout}] = ne@handle(varargin{:});
         end

         function varargout = notify(varargin)   
             [varargout{1:nargout}] = notify@handle(varargin{:});
         end
         
         % Sealed
         %function varargout = isvalid(varargin)
         %   [varargout{1:nargout}] = isvalid@handle(varargin{:});
         %end
     
     end
     
     
     
end