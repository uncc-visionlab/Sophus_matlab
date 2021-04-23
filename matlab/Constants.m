classdef Constants < handle
    % Sophus constants    
    properties (Constant)
        M_PI    = 3.14159265358979323846264338328      %/* pi */
        epsilon = 1e-10
    end
    methods (Static)
        function retval = pi()
            retval = Constants.M_PI;
        end
        %function retval = epsilon()
        %    retval = epsilon;
        %end
    end
end

