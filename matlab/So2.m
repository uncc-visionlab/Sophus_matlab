classdef So2 < RotationMatrix
    %  SO2 base type - implements SO2 class but is storage agnostic.
    %
    %  SO(2) is the group of rotations in 2d. As a matrix group, it is the set of
    %  matrices which are orthogonal such that ``R * R' = I`` (with ``R'`` being the
    %  transpose of ``R``) and have a positive determinant. In particular, the
    %  determinant is 1. Let ``theta`` be the rotation angle, the rotation matrix
    %  can be written in close form:
    %
    %   | cos(theta) -sin(theta) |
    %   | sin(theta)  cos(theta) |
    %
    %  As a matter of fact, the first column of those matrices is isomorph to the
    %  set of unit complex numbers U(1). Thus, internally, SO2 is represented as
    %  complex number with length 1.
    %
    %  SO(2) is a compact and commutative group. First it is compact since the set
    %  of rotation matrices is a closed and bounded set. Second it is commutative
    %  since ``R(alpha) * R(beta) = R(beta) * R(alpha``,  simply because ``alpha +
    %  beta = beta + alpha`` with ``alpha`` and ``beta`` being rotation angles
    %  (about the same axis).
    %
    %  Class invairant: The 2-norm of ``unit_complex`` must be close to 1.
    %  Technically speaking, it must hold that:
    %
    %    ``|unit_complex().squaredNorm() - 1| <= Constants<Scalar>::epsilon()``.

    properties (Constant)
        % Degrees of freedom of manifold, number of dimensions in tangent space (one
        % since we only have in-plane rotations).
        DoF = 1;
        % Number of internal parameters used (complex numbers are a tuples).
        num_parameters = 2;
        % Group transformations are 2x2 matrices.
        N = 2;
    end

    properties
        z
    end

    methods (Static)
        %////////////////////////////////////////////////////////////////////////////
        %// public static functions
        %////////////////////////////////////////////////////////////////////////////


        function retval =  exp(theta)
            % Group exponential
            %
            % This functions takes in an element of tangent space (= rotation angle
            % ``theta``) and returns the corresponding element of the group SO(2).
            %
            % To be more specific, this function computes ``expmat(hat(omega))``
            % with ``expmat(.)`` being the matrix exponential and ``hat(.)`` being the
            % hat()-operator of SO(2).
            %
            retval = So2(complex(cos(theta) + 1i*sin(theta)));
        end
        function retval = Dx_exp_x(theta)
            % /// Returns derivative of exp(x) wrt. x.
            % ///
            retval = [-sin(theta); cos(theta)];
        end
        
        function retval = Dx_exp_x_at_0()
            % /// Returns derivative of exp(x) wrt. x_i at x=0.
            % ///
            retval = [0; 1];
        end

        function retval = Dxi_exp_x_matrix_at_0(i)
            % /// Returns derivative of exp(x).matrix() wrt. ``x_i at x=0``.
            % ///
            retval = generator();
        end
        
        function retval = generator()
            % Returns the infinitesimal generators of SO3.
            %
            % The infinitesimal generators of SO(2) is:
            %
            %   |  0  1 |
            %   | -1  0 |
            %
            retval = hat(1);
        end
        function retval =  hat(theta)
            % hat-operator
            %
            % It takes in the scalar representation ``theta`` (= rotation angle) and
            % returns the corresponding matrix representation of Lie algebra element.
            %
            % Formally, the ``hat()`` operator of SO(2) is defined as
            %
            %   ``hat(.): R^2 -> R^{2x2},  hat(theta) = theta * G``
            %
            % with ``G`` being the infinitesimal generator of SO(2).
            %
            % The corresponding inverse is the ``vee``-operator, see below.
            %
            retval = [[0, -theta]; [theta, 0]];
        end
        function retval = fitToSO2(R)
            % /// Returns closed SO2 given arbitrary 2x2 matrix.
            % ///
            retval = So2.makeRotationMatrix(R);
        end
        function retval = lieBracket(a, b)
            % Lie bracket
            %
            % It returns the Lie bracket of SO(2). Since SO(2) is a commutative group,
            % the Lie bracket is simply ``0``.
            %
            retval = 0;
        end
        function retval = log(z)
            % Logarithmic map
            %
            % Computes the logarithm, the inverse of the group exponential which maps
            % element of the group (rotation matrices) to elements of the tangent space
            % (rotation angles).
            %
            % To be specific, this function computes ``vee(logmat(.))`` with
            % ``logmat(.)`` being the matrix logarithm and ``vee(.)`` the vee-operator
            % of SO(2).
            %
            retval = atan2(imag(z), real(z));
        end
        function retval = sampleUniform()
            % Draw uniform sample from SO(2) manifold.
            % 
            uniform = (rand(1)-0.5)*2*Constants.M_PI;
            retval = So2(uniform);
        end
        function retval = vee(omega)
            %  vee-operator
            %
            %  It takes the 2x2-matrix representation ``Omega`` and maps it to the
            %  corresponding scalar representation of Lie algebra.
            %
            %  This is the inverse of the hat-operator, see above.
            %
            %  Precondition: ``Omega`` must have the following structure:
            %
            %                 |  0 -a |
            %                 |  a  0 |
            %
            if (norm(diag(omega),1) > Constants.epsilon)
                fprintf(1,"Omega: \n% %f", omega);
            end
            if (abs(omega(2,1) + omega(1,2)) < Constants.epsilon)
                fprintf(1,"Omega: \n% %f", omega);
            end
            retval = omega(2, 1);
        end
    end

    methods
        function self = So2(varargin)
            % internally represented by a unit complex number z
            if (nargin == 1 && isequal(class(varargin{1}),'So2'))
                self = So2(varargin{1}.unit_complex());
            elseif (nargin == 1 && numel(varargin{1}) == 1 && imag(varargin{1}) > 0)
                % input is a complex number
                self.z = complex(varargin{1});
            elseif (nargin == 1 && numel(varargin{1}) == 1)
                % input is an angle theta in radians
                self = So2.exp(varargin{1});
            elseif (nargin == 1 && numel(varargin{1}) == 4 && size(varargin{1},1) == 2 && size(varargin{1},2) == 2)
                % input is a 2x2 rotation matrix
                self.real = complex(0.5*(varargin{1}(1,1) + varargin{1}(2,2)) + ...
                    1i*(0.5*(varargin{1}(2,1) + varargin{1}(1,2))));
                % TODO
                % SOPHUS_ENSURE(isOrthogonal(R), "R is not orthogonal:\n %", R);
                % SOPHUS_ENSURE(R.determinant() > Scalar(0), "det(R) is not positive: %",
                % R.determinant());
                if (abs(det(varargin{1}) - 1) <= Constants.epsilon)
                    fprintf(1,"det(R) should be (close to) 1.\n")
                end
            elseif (nargin == 2)
                % input is a pair of arguments arg1 + i*arg2
                self.z = complex(varargin{1} + 1i*varargin{2});
            else
                fprintf(1,"Error could not construct So2 object.\n");
            end
        end

        function normalize(self)
            % normalize the complex number to unit magnitude
            length = abs(self.z);
            if (length < Constants.epsilon)
                fprintf(1,"Complex number should not be close to zero!\n");
            end
            self.z = self.z/length;
        end

        function retval = Adj(self)
            % It simply ``1``, since ``SO(2)`` is a commutative group.
            retval = 1;
        end

        function retval = inverse(self)
            % Returns the ith generator of internal U(1) representation.
            retval = So2(conj(self.z));
        end

        function retval = exp_(self)
            retval = So2.exp(self.z);
        end

        function retval = log_(self)
            retval = So2.log(self.z);
        end

        function retval = print(self)
            retval = sprintf("So2: %f + i%f",real(self.z),imag(self.z))
        end

        function retval = matrix(self)
            % Returns 2x2 matrix representation of the instance.
            % 
            % For SO(2), the matrix representation is an orthogonal matrix ``R`` with
            % ``det(R)=1``, thus the so-called "rotation matrix".
            % 
            retval = [real(self.z), -imag(self.z); imag(self.z),  real(self.z)];
        end

        function retval = mtimes(self, other)
            % overload function for the * operator in MATLAB
            if (isequal(class(other),'So2') == true)
                % Group multiplication, which is rotation concatenation.
                % We can assume that the squared-norm is close to 1 since we deal with a
                % unit complex number. Due to numerical precision issues, there might
                % be a small drift after pose concatenation. Hence, we need to renormalizes
                % the complex number here.
                % Since squared-norm is close to 1, we do not need to calculate the costly
                % square-root, but can use an approximation around 1 (see
                % http://stackoverflow.com/a/12934750 for details).
                retval = So2(self.unit_complex()*other.unit_complex());
                retval.normalize();
            elseif numel(other) == 2
                % Group action on 2-points.
                %
                % This function rotates a 2 dimensional point ``p`` by the SO2 element
                %  ``bar_R_foo`` (= rotation matrix): ``p_bar = bar_R_foo * p_foo``.
                %
                retval = other;
                retval(1) = real(self.z)*other(1) - imag(self.z)*other(2);
                retval(2) = imag(self.z)*other(1) + real(self.z)*other(2);
            end
        end
        
        % TODO:
        % /// Group action on homogeneous 2-points.
        % ///
        % /// This function rotates a homogeneous 2 dimensional point ``p`` by the SO2
        % /// element ``bar_R_foo`` (= rotation matrix): ``p_bar = bar_R_foo * p_foo``.
        % ///
        % template <typename HPointDerived,
        % typename = typename std::enable_if<
        %     IsFixedSizeVector<HPointDerived, 3>::value>::type>
        % SOPHUS_FUNC HomogeneousPointProduct<HPointDerived> operator*(
        % Eigen::MatrixBase<HPointDerived> const& p) const {
        % Scalar const& real = unit_complex().x();
        % Scalar const& imag = unit_complex().y();
        % return HomogeneousPointProduct<HPointDerived>(
        % real * p[0] - imag * p[1], imag * p[0] + real * p[1], p[2]);
        % }
        %
        % /// Group action on lines.
        % ///
        % /// This function rotates a parametrized line ``l(t) = o + t * d`` by the SO2
        % /// element:
        % ///
        % /// Both direction ``d`` and origin ``o`` are rotated as a 2 dimensional point
        % ///
        % SOPHUS_FUNC Line operator*(Line const& l) const {
        % return Line((*this) * l.origin(), (*this) * l.direction());
        % }
        
        function retval = Dx_this_mul_exp_x_at_0(self)
            % Returns derivative of  this * SO2::exp(x)  wrt. x at x=0.
            % 
            retval = [-imag(z); real(z)];
        end

        function retval = params(self)
            % Returns internal parameters of SO(2).
            %
            % It returns (c[0], c[1]), with c being the unit complex number.
            %
            retval = unit_complex();
        end
        
        function setComplex(self, z)
            %  Takes in complex number / tuple and normalizes it.
            %
            %  Precondition: The complex number must not be close to zero.
            %
            self.z = z;
            self.normalize();
        end

        function retval = unit_complex(self)
            retval = self.z;
        end
    end
end
