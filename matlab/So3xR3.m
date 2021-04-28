classdef So3xR3 < Constants
    %  SO3 x R3 base type - implements SO3 x R3 class.
    %
    %

    properties (Constant)
        % Degrees of freedom of manifold, number of dimensions in tangent space
        % (three for translation, three for rotation).
        DoF = 6;
        % Number of internal parameters used (4-tuple for quaternion, three for
        % translation).
        num_parameters = 7;
        % Group transformations are 4x4 matrices.
        N = 4;
    end

    properties
        so3
        translation
    end

    methods (Static)
        %////////////////////////////////////////////////////////////////////////////
        %// public static functions
        %////////////////////////////////////////////////////////////////////////////


        function J = Dx_exp_x(omega)
            %  Returns derivative of exp(x) wrt. x.
            %
            J_So3 = So3.Dx_exp_x(omega(4:6));
            J_R3 = eye(3);
            J = [zeros(4,3), J_So3;
                J_R3, zeros(3,3)];
        end        
        function J = Dx_exp_x_at_0()
            %  Returns derivative of exp(x) wrt. x_i at x=0.
            %
            %  clang-format off
            J_So3 = So3.Dx_exp_x_at_0();
            J_R3 = eye(3);
            J = [zeros(4,3), J_So3;
                J_R3, zeros(3,3)];
            %  clang-format on
        end

        function retval = Dxi_exp_x_matrix_at_0(i)
            %  Returns derivative of exp(x).matrix() wrt. ``x_i at x=0``.
            %
            retval = generator(i);
        end

        function retval = pseudo_exp(omega)
            %  Group exponential of only the SO3 component
            %
            %  This functions takes in an element of tangent space (= rotation vector
            %  ``omega``) and returns the corresponding element of the group SO(3).
            %
            %  To be more specific, this function computes ``expmat(hat(omega))``
            %  with ``expmat(.)`` being the matrix exponential and ``hat(.)`` being the
            %  hat()-operator of SO(3).
            %
            [so3_group, ~] = So3.expAndTheta(omega(4:6));
            tangent = [omega(1:3); so3_group.log()];
            retval = So3xR3(tangent);
        end

        function retval = fitToSO3xR3(T)
            %  Returns closest SO3xR3 given arbitrary 4x4 matrix.
            % 
            retval = So3xR3( So3.fitToSO3(T(1:3,1:3)), T(1:3,4));
        end
        function retval = generator(i)
            %  Returns the ith infinitesimal generators of SE(3).
            % 
            %  The infinitesimal generators of SE(3) are:
            % 
            %  ```
            %          |  0  0  0  1 |
            %    G_0 = |  0  0  0  0 |
            %          |  0  0  0  0 |
            %          |  0  0  0  0 |
            % 
            %          |  0  0  0  0 |
            %    G_1 = |  0  0  0  1 |
            %          |  0  0  0  0 |
            %          |  0  0  0  0 |
            % 
            %          |  0  0  0  0 |
            %    G_2 = |  0  0  0  0 |
            %          |  0  0  0  1 |
            %          |  0  0  0  0 |
            % 
            %          |  0  0  0  0 |
            %    G_3 = |  0  0 -1  0 |
            %          |  0  1  0  0 |
            %          |  0  0  0  0 |
            % 
            %          |  0  0  1  0 |
            %    G_4 = |  0  0  0  0 |
            %          | -1  0  0  0 |
            %          |  0  0  0  0 |
            % 
            %          |  0 -1  0  0 |
            %    G_5 = |  1  0  0  0 |
            %          |  0  0  0  0 |
            %          |  0  0  0  0 |
            %  ```
            % 
            %  Precondition: ``i`` must be in [0, 5].
            % 

            if (i < 0 || i > 5)
                fprintf(1,"i should be an integer in range [0,5].\n");
            end
            e = zeros(1,N);
            e(i+1) = 1;
            retval = So3xR3.pseudo_hat(e);
        end
        function retval = pseudo_hat(a)
            %  hat-operator
            % 
            %  It takes in the 6-vector representation (= twist) and returns the
            %  corresponding matrix representation of Lie algebra element.
            % 
            %  Formally, the hat()-operator of SE(3) is defined as
            % 
            %    ``hat(.): R^6 -> R^{4x4},  hat(a) = sum_i a_i * G_i``  (for i=0,...,5)
            % 
            %  with ``G_i`` being the ith infinitesimal generator of SE(3).
            % 
            %  The corresponding inverse is the vee()-operator, see below.
            % 
            retval = [So3.hat(a(4:6)), [a(1); a(2); a(3)]; 0 0 0 0];
        end
        function retval = rotX(x)
            %  Construct x-axis rotation.
            %
            retval = So3.exp([x; (0); (0)]);
        end
        function retval =  rotY(y)
            %  Construct y-axis rotation.
            %
            retval = So3.exp([(0); y; (0)]);
        end
        function retval =  rotZ(z)
            %  Construct z-axis rotation.
            %
            retval = So3.exp([(0); (0); z]);
        end
    end

    methods
        function self = So3xR3(varargin)
            % internally represented by a unit complex number z
            if (nargin == 1 && isequal(class(varargin{1}),'So3xR3'))
                self.translation = varargin{1}.getTranslation();
                self.so3 = varargin{1}.getSo3();
            elseif (nargin == 1 && numel(varargin{1}) == 6)
                % input is an tangent vector
                self.translation = [varargin{1}(1); varargin{1}(2); varargin{1}(3)];
                self.so3 = So3.exp(varargin{1}(4:6));
            elseif (nargin == 2 && isequal(class(varargin{1}),'So3') && size(varargin{2},1) == 3)
                % constructor from So3() object and translation vector
                self.so3 = varargin{1};
                self.translation = [varargin{2}(1); varargin{2}(2); varargin{2}(3)];
            elseif (nargin == 2 && numel(varargin{1}) == 9 && size(varargin{1},1) == 3 && size(varargin{1},2) == 3 && ...
                    numel(varargin{2}) == 3 && size(varargin{2},1) == 3 && size(varargin{2},2) == 1)
                % input is a 3x3 rotation matrix and 3x1 translation vector
                self.so3 = So3(varargin{1}(1:3,1:3));
                self.translation = [varargin{2}(1); varargin{2}(2); varargin{2}(3)];
            else
                fprintf(1,"Error could not construct Se2 object.\n");
            end
        end

        function retval = pseudo_Adj(self)
            % Adjoint transformation
            %
            %  This function return the adjoint transformation ``Ad`` of the group
            %  element ``A`` such that for all ``x`` it holds that
            %  ``hat(Ad_A * x) = A * hat(x) A^{-1}``. See hat-operator below.
            %
            R = self.so3().matrix();
            res = [R, So3.hat(self.translation) * R;
                zeros(3,3), R];
        end

        function retval = angleX(self)
            %  Extract rotation angle about canonical X-axis
            % 
            retval = self.so3().angleX();
        end

        function retval = angleY(self)
            %  Extract rotation angle about canonical Y-axis
            % 
            retval = self.so3().angleY(self)
        end

        function retval = angleZ(self)
            %  Extract rotation angle about canonical Z-axis
            % 
            retval = self.so3().angleZ();
        end
        
        function J = Dx_this_mul_exp_x_at_0(self)
            %  Returns derivative of  this * SO3::exp(x)  wrt. x at x=0.
            %
            J_so3 = self.so3.Dx_this_mul_exp_x_at_0();
            J_R3 = eye(3);
            J = [zeros(4,3), J_so3;
                J_R3, zeros(3,3)];
        end

        function retval = inverse(self)
            %  Returns group inverse.
            invR = self.so3.inverse();
            retval = So3xR3(invR, invR*(self.translation * -1));
        end

        function retval = pseudo_log(self)
            %  Logarithmic map
            %
            %  Computes the logarithm, the inverse of the group exponential which maps
            %  element of the group (rotation matrices) to elements of the tangent space
            %  (rotation-vector).
            %
            %  To be specific, this function computes ``vee(logmat(.))`` with
            %  ``logmat(.)`` being the matrix logarithm and ``vee(.)`` the vee-operator
            %  of SO(3).
            %
            retval = [self.translation; self.so3.log()];
        end

        function normalize(self)
            %  It re-normalizes the SO3 element.
            % 
            %  Note: Because of the class invariant of SO3, there is typically no need to
            %  call this function directly.
            % 
            self.so3.normalize();
        end
        
        function retval = toString(self)
            retval = sprintf("So3xR3: (%3g, %3g, %3g, %s)", self.translation(1), ...
                self.translation(2), self.translation(3), self.so3.toString());
        end

        function homogeneous_matrix = matrix(self)
            % Returns 4x4 matrix representation of the instance.
            %
            % It has the following form:
            %
            %    | R t |
            %    | o 1 |
            %
            % where ``R`` is a 3x3 rotation matrix, ``t`` a translation 3-vector and
            % ``o`` a 3-column vector of zeros.
            %
            homogeneous_matrix = [ self.matrix3x4(); (0), (0), (0), (1)];
        end

        function matrix = matrix3x4(self)
            %  Returns the significant first three rows of the matrix above.
            % 
            matrix = [self.rotationMatrix(), self.translation];
        end

        function retval = mtimes(self, other)
            % overload function for the * operator in MATLAB
            if (isequal(class(other),'So3xR3') == true)
                retval = So3xR3(self.so3 * other.so3, self.getTranslation() + self.so3 * other.getTranslation());
            elseif (isscalar(other) == true)
                retval = So3xR3.pseudo_exp([self.translation * other; self.so3.log() * other]);
            elseif numel(other) == 3
                %  Group action on 3-points.
                %
                %  This function rotates and translates a three dimensional point ``p`` by the
                %  SO(3) element ``bar_T_foo = (bar_R_foo, t_bar)`` (= rigid body
                %  transformation):
                %
                %    ``p_bar = bar_R_foo * p_foo + t_bar``.
                %
                retval = self.so3 * other + self.translation;
            elseif numel(other) == 4
                %  Group action on homogeneous 3-points. See above for more details.
                %
                retval = [self.so3 * other(1:3) + other(4) * self.translation(); other(4)];
            elseif false
                %  Group action on lines.
                % 
                %  This function rotates and translates a parametrized line
                %  ``l(t) = o + t * d`` by the SE(3) element:
                % 
                %  Origin is transformed using SE(3) action
                %  Direction is transformed using rotation part
                % 
                % SOPHUS_FUNC Line operator*(Line const& l) const {
                % return Line((*this) * l.origin(), so3() * l.direction());
                % }
            else
                fprintf(1,"Could not multiply So3xR3 object.");
            end
        end

        function retval = rotationMatrix(self)
            %  Returns rotation matrix.
            %
            retval = self.so3.matrix();
        end

        function setTranslation(self, trans)
            %  Sets the translation component using trans
            self.translation = trans;
        end
        
        function setQuaternion(self, quat)
            %  Takes in complex number / tuple and normalizes it.
            %
            %  Precondition: The complex number must not be close to zero.
            %
            self.so3.setQuaternion(quat);
        end

        function setRotationMatrix(self, R)
            %  Sets ``so3`` using ``rotation_matrix``.
            %
            %  Precondition: ``R`` must be orthogonal and ``det(R)=1``.
            %
            % TODO:
            % SOPHUS_ENSURE(isOrthogonal(R), "R is not orthogonal:\n %", R);
            % SOPHUS_ENSURE(R.determinant() > Scalar(0), "det(R) is not positive: %",
            % R.determinant());
            self.so3.setQuaternion(quaternion(rotm2quat(eye(3))));
        end

        function retval = getSo3(self)
            %  Accessor of SO3 group.
            %
            retval = self.so3;
        end

        function retval = getTranslation(self)
            %  Accessor of translation vector
            %
            retval = self.translation;
        end
    end
end


