classdef Constants < handle
    % Sophus constants    
    properties (Constant)
        M_PI    = 3.14159265358979323846264338328      %/* pi */
        % could use MATLAB's eps('double') or eps('single') when made
        % data type agnostic Eigen often uses sqrt(eps('double')) for
        % epsilon in Numerical differentiation (line 72)
        % https://github.com/libigl/eigen/blob/master/unsupported/Eigen/src/NumericalDiff/NumericalDiff.h
        % that works out to about 1.5e-8
        epsilon = (eps('double')^(1/3)) 
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

