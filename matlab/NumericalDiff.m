classdef NumericalDiff < Constants
    % Numerically calculate the Jacobian of the function f()
    % with respect to the input arguments x
    properties (Constant)
        CENTRAL_DIFFERENCE = 1
        FORWARD_DIFFERENCE = 2
    end
    properties
        f
        inputs
        values
        method
    end
    
    methods (Static)
        function J = demo()
            sqr = @(n) n.^2;
            num_diff = NumericalDiff(sqr, 1, 1);
            for x = 0:.1:2
                J = num_diff.df(x)
            end
        end
    end
    
    methods
        function self = NumericalDiff(fhandle, dim_input, dim_value, method)
            % fhandle = function handle of function to differentiate
            % dim_input = the number of values provided as input to the
            % function to differentiate
            % dim_output = the number of values provided as output from the
            % function to differentiate
            % method = 'Forward' (default) or 'Central'
            self.f = fhandle;
            self.inputs = dim_input;
            self.values = dim_value;
            if (nargin == 4 && strcmpi(method,'central') == true)
                self.method = NumericalDiff.CENTRAL_DIFFERENCE;
            else
                self.method = NumericalDiff.FORWARD_DIFFERENCE;
            end
        end
        
        function J = df(self, xin)
            %nfev=0;
            epsilon = Constants.epsilon; % sqrt(eps('double'));
            val1 = zeros(self.inputs,1);
            val2 = zeros(self.inputs,1);
            J = zeros(self.values, self.inputs);
            x = xin;
            % initialization
            switch(self.method)
                case NumericalDiff.FORWARD_DIFFERENCE
                    % compute f(x)
                    val1 = self.f(x);
                    %nfev = nfev + 1;
                case NumericalDiff.CENTRAL_DIFFERENCE
                    % do nothing
                otherwise
                    printf(1,'NumericalDiff: Error no such method!\n');
            end
            % Function Body
            for dim=1:self.inputs
                h = epsilon*abs(x(dim));
                if (h==0)
                    h = epsilon;
                end
                switch(self.method)
                    case NumericalDiff.FORWARD_DIFFERENCE
                        x(dim) = x(dim) + h;
                        val2 = self.f(x);
                        %nfev = nfev + 1;
                        x(dim) = xin(dim);
                        J(:,dim) = (val2-val1)/h;
                    case NumericalDiff.CENTRAL_DIFFERENCE
                        x(dim) = x(dim) + h;
                        val2 = self.f(x);
                        %nfev = nfev + 1;
                        x(dim) = x(dim) - 2*h;
                        val1 = self.f(x);
                        %nfev = nfev + 1;
                        x(dim) = xin(dim);
                        J(:,dim) = (val2-val1)/(2*h);
                    otherwise
                        printf(1,'NumericalDiff: Error no such method!\n');
                end
            end
        end
    end
end
    
