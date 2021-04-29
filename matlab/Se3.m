classdef Se3 < Constants
    % SE3 base type - implements SE3 class but is storage agnostic.
    %
    % SE(3) is the group of rotations  and translation in 3d. It is the
    % semi-direct product of SO(3) and the 3d Euclidean vector space.  The class
    % is represented using a composition of SO3  for rotation and a one 3-vector
    % for translation.
    %
    % SE(3) is neither compact, nor a commutative group.
    %
    % See SO3 for more details of the rotation representation in 3d.
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
        % ////////////////////////////////////////////////////////////////////////////
        % // public static functions
        % ////////////////////////////////////////////////////////////////////////////

        function J = Dx_exp_x(upsilon_omega)
            %  Returns derivative of exp(x) wrt. x.
            % 
            upsilon = upsilon_omega(1:3);
            omega = upsilon_omega(4:6);

            c0 = omega(1) * omega(1);
            c1 = omega(2) * omega(2);
            c2 = omega(3) * omega(3);
            c3 = c0 + c1 + c2;
            o = (0);
            h = (0.5);
            i = (1);

            if (c3 < Constants.epsilon)
                ux = (0.5) * upsilon(1);
                uy = (0.5) * upsilon(2);
                uz = (0.5) * upsilon(3);

                %  clang-format off
                J = [o, o, o, h, o, o;
                    o, o, o, o, h, o;
                    o, o, o, o, o, h;
                    o, o, o, o, o, o;
                    i, o, o, o, uz, -uy;
                    o, i, o, -uz, o, ux;
                    o, o, i, uy, -ux, o];
                %  clang-format on
                return;
            end

            c4 = sqrt(c3);
            c5 = (1.0) / c4;
            c6 = (0.5) * c4;
            c7 = sin(c6);
            c8 = c5 * c7;
            c9 = pow(c3, -3.0 / 2.0);
            c10 = c7 * c9;
            c11 = (1.0) / c3;
            c12 = cos(c6);
            c13 = (0.5) * c11 * c12;
            c14 = c7 * c9 * omega(1);
            c15 = (0.5) * c11 * c12 * omega(1);
            c16 = -c14 * omega(2) + c15 * omega(2);
            c17 = -c14 * omega(3) + c15 * omega(3);
            c18 = omega(2) * omega(3);
            c19 = -c10 * c18 + c13 * c18;
            c20 = c5 * omega(1);
            c21 = (0.5) * c7;
            c22 = c5 * omega(2);
            c23 = c5 * omega(3);
            c24 = -c1;
            c25 = -c2;
            c26 = c24 + c25;
            c27 = sin(c4);
            c28 = -c27 + c4;
            c29 = c28 * c9;
            c30 = cos(c4);
            c31 = -c30 + (1);
            c32 = c11 * c31;
            c33 = c32 * omega(3);
            c34 = c29 * omega(1);
            c35 = c34 * omega(2);
            c36 = c32 * omega(2);
            c37 = c34 * omega(3);
            c38 = pow(c3, -5.0 / 2.0);
            c39 = (3) * c28 * c38 * omega(1);
            c40 = c26 * c9;
            c41 = -c20 * c30 + c20;
            c42 = c27 * c9 * omega(1);
            c43 = c42 * omega(2);
            c44 = pow(c3, -2);
            c45 = (2) * c31 * c44 * omega(1);
            c46 = c45 * omega(2);
            c47 = c29 * omega(3);
            c48 = c43 - c46 + c47;
            c49 = (3) * c0 * c28 * c38;
            c50 = c9 * omega(1) * omega(3);
            c51 = c41 * c50 - c49 * omega(3);
            c52 = c9 * omega(1) * omega(2);
            c53 = c41 * c52 - c49 * omega(2);
            c54 = c42 * omega(3);
            c55 = c45 * omega(3);
            c56 = c29 * omega(2);
            c57 = -c54 + c55 + c56;
            c58 = (-2) * c56;
            c59 = (3) * c28 * c38 * omega(2);
            c60 = -c22 * c30 + c22;
            c61 = -c18 * c39;
            c62 = c32 + c61;
            c63 = c27 * c9;
            c64 = c1 * c63;
            c65 = (2) * c31 * c44;
            c66 = c1 * c65;
            c67 = c50 * c60;
            c68 = -c1 * c39 + c52 * c60;
            c69 = c18 * c63;
            c70 = c18 * c65;
            c71 = c34 - c69 + c70;
            c72 = (-2) * c47;
            c73 = (3) * c28 * c38 * omega(3);
            c74 = -c23 * c30 + c23;
            c75 = -c32 + c61;
            c76 = c2 * c63;
            c77 = c2 * c65;
            c78 = c52 * c74;
            c79 = c34 + c69 - c70;
            c80 = -c2 * c39 + c50 * c74;
            c81 = -c0;
            c82 = c25 + c81;
            c83 = c32 * omega(1);
            c84 = c18 * c29;
            c85 = (-2) * c34;
            c86 = c82 * c9;
            c87 = c0 * c63;
            c88 = c0 * c65;
            c89 = c9 * omega(2) * omega(3);
            c90 = c41 * c89;
            c91 = c54 - c55 + c56;
            c92 = -c1 * c73 + c60 * c89;
            c93 = -c43 + c46 + c47;
            c94 = -c2 * c59 + c74 * c89;
            c95 = c24 + c81;
            c96 = c9 * c95;
            J(1, 1) = o;
            J(1, 2) = o;
            J(1, 3) = o;
            J(1, 4) = -c0 * c10 + c0 * c13 + c8;
            J(1, 5) = c16;
            J(1, 6) = c17;
            J(2, 1) = o;
            J(2, 2) = o;
            J(2, 3) = o;
            J(2, 4) = c16;
            J(2, 5) = -c1 * c10 + c1 * c13 + c8;
            J(2, 6) = c19;
            J(3, 1) = o;
            J(3, 2) = o;
            J(3, 3) = o;
            J(3, 4) = c17;
            J(3, 5) = c19;
            J(3, 6) = -c10 * c2 + c13 * c2 + c8;
            J(4, 1) = o;
            J(4, 2) = o;
            J(4, 3) = o;
            J(4, 4) = -c20 * c21;
            J(4, 5) = -c21 * c22;
            J(4, 6) = -c21 * c23;
            J(5, 1) = c26 * c29 + (1);
            J(5, 2) = -c33 + c35;
            J(5, 3) = c36 + c37;
            J(5, 4) = upsilon(1) * (-c26 * c39 + c40 * c41) + upsilon(2) * (c53 + c57) + ...
            upsilon(3) * (c48 + c51);
            J(5, 5) = upsilon(1) * (-c26 * c59 + c40 * c60 + c58) + ...
                upsilon(2) * (c68 + c71) + upsilon(3) * (c62 + c64 - c66 + c67);
            J(5, 6) = upsilon(1) * (-c26 * c73 + c40 * c74 + c72) + ...
                upsilon(2) * (c75 - c76 + c77 + c78) + upsilon(3) * (c79 + c80);
            J(6, 1) = c33 + c35;
            J(6, 2) = c29 * c82 + (1);
            J(6, 3) = -c83 + c84;
            J(6, 4) = upsilon(1) * (c53 + c91) + ...
                upsilon(2) * (-c39 * c82 + c41 * c86 + c85) + ...
                upsilon(3) * (c75 - c87 + c88 + c90);
            J(6, 5) = upsilon(1) * (c68 + c79) + upsilon(2) * (-c59 * c82 + c60 * c86) + ...
                upsilon(3) * (c92 + c93);
            J(6, 6) = upsilon(1) * (c62 + c76 - c77 + c78) + ...
                upsilon(2) * (c72 - c73 * c82 + c74 * c86) + ...
                upsilon(3) * (c57 + c94);
            J(7, 1) = -c36 + c37;
            J(7, 2) = c83 + c84;
            J(7, 3) = c29 * c95 + (1);
            J(7, 4) = upsilon(1) * (c51 + c93) + upsilon(2) * (c62 + c87 - c88 + c90) + ...
                upsilon(3) * (-c39 * c95 + c41 * c96 + c85);
            J(7, 5) = upsilon(1) * (-c64 + c66 + c67 + c75) + upsilon(2) * (c48 + c92) + ...
                upsilon(3) * (c58 - c59 * c95 + c60 * c96);
            J(7, 6) = upsilon(1) * (c71 + c80) + upsilon(2) * (c91 + c94) + ...
                upsilon(3) * (-c73 * c95 + c74 * c96);
        end

        function J = Dx_exp_x_at_0()
            %  Returns derivative of exp(x) wrt. x_i at x=0.
            % 
            o = (0);
            h = (0.5);
            i = (1);

            % // clang-format off
            J = [o, o, o, h, o, o;
                o, o, o, o, h, o;
                o, o, o, o, o, h;
                o, o, o, o, o, o;
                i, o, o, o, o, o;
                o, i, o, o, o, o;
                o, o, i, o, o, o];
            % // clang-format on
        end
        function retval = Dxi_exp_x_matrix_at_0(i)
            %  Returns derivative of exp(x).matrix() wrt. ``x_i at x=0``.
            % 
            retval = Se3.generator(i);
        end
        function retval = exp(a)
            %  Group exponential
            %
            %  This functions takes in an element of tangent space (= twist ``a``) and
            %  returns the corresponding element of the group SE(3).
            %
            %  The first three components of ``a`` represent the translational part
            %  ``upsilon`` in the tangent space of SE(3), while the last three components
            %  of ``a`` represents the rotation vector ``omega``.
            %  To be more specific, this function computes ``expmat(hat(a))`` with
            %  ``expmat(.)`` being the matrix exponential and ``hat(.)`` the hat-operator
            %  of SE(3), see below.
            %
            omega = a(4:6);
            [so3, theta] = So3.expAndTheta(omega);
            Omega = So3.hat(omega);
            Omega_sq = Omega * Omega;
            if (theta < Constants.epsilon)
                V = so3.matrix();
                %  Note: That is an accurate expansion!
            else
                theta_sq = theta * theta;
                V = (eye(3) + ((1) - cos(theta)) / (theta_sq) * Omega + (theta - sin(theta)) / (theta_sq * theta) * Omega_sq);
            end
            retval = Se3( so3, V * a(1:3));
        end

        function retval = fitToSE3(T)
            %  Returns closest SE3 given arbitrary 4x4 matrix.
            % 
            retval = Se3( So3.fitToSO3(T(1:3,1:3)), T(1:3,4));
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
            retval = Se3.hat(e);
        end
        function retval = hat(a)
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
        function retval = lieBracket( a, b)
            %  Lie bracket
            %
            %  It computes the Lie bracket of SE(3). To be more specific, it computes
            %
            %    ``[omega_1, omega_2]_se3 := vee([hat(omega_1), hat(omega_2)])``
            %
            %  with ``[A,B] := AB-BA`` being the matrix commutator, ``hat(.) the
            %  hat-operator and ``vee(.)`` the vee-operator of SE(3).
            %
            upsilon1 = a(1:3);
            upsilon2 = b(1:3);
            omega1 = a(4:6);
            omega2 = b(4:6);

            retval = [cross(omega1, upsilon2) + cross(upsilon1, omega2);
                cross(omega1, omega2)];
        end
        function retval = rotX(x)
            %  Construct x-axis rotation.
            % 
            retval = Se3(So3.rotX(x), [0; 0; 0]);
        end

        function retval = rotY(y)
            %  Construct y-axis rotation.
            % 
            retval = Se3(So3.rotY(y), [0; 0; 0]);
        end

        function retval = rotZ(z)
            %  Construct z-axis rotation.
            % 
            retval = Se3(So3.rotZ(z), [0; 0; 0]);
        end

        function upsilon_omega = vee(Omega)
            %  vee-operator
            % 
            %  It takes 4x4-matrix representation ``Omega`` and maps it to the
            %  corresponding 6-vector representation of Lie algebra.
            % 
            %  This is the inverse of the hat()-operator, see above.
            % 
            %  Precondition: ``Omega`` must have the following structure:
            % 
            %                 |  0 -f  e  a |
            %                 |  f  0 -d  b |
            %                 | -e  d  0  c |
            %                 |  0  0  0  0 | .
            % 
            upsilon_omega = [Omega(1:3,4); So3.vee(Omega(1:3,1:3))];
        end
    end

    methods
        function self = Se3(varargin)
            % internally represented by a unit complex number z
            if (nargin == 1 && isequal(class(varargin{1}),'Se3'))
                self.translation = varargin{1}.getTranslation();
                self.so3 = varargin{1}.getSo3();
            elseif (nargin == 1 && numel(varargin{1}) == 6)
                % input is an tangent vector from a Lie algebra
                % we must immediately convert it to a Lie group using
                % the exponential map
                %self.translation = [varargin{1}(1); varargin{1}(2); varargin{1}(3)];
                %self.so3 = So3.exp(varargin{1}(4:6));
                self = Se3.exp(varargin{1});
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

        function res = Adj(self)
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
            %  Returns derivative of  this * exp(x)  wrt x at x=0.
            % 
            q = self.so3.unit_quaternion();
            [qw, qx, qy, qz] = q.parts();
            c0 = (0.5) * qw;
            c1 = (0.5) * qz;
            c2 = -c1;
            c3 = (0.5) * qy;
            c4 = (0.5) * qx;
            c5 = -c4;
            c6 = -c3;
            c7 = qw * qw;
            c8 = qx * qx;
            c9 = qy * qy;
            c10 = -c9;
            c11 = qz * qz;
            c12 = -c11;
            c13 = (2) * qw;
            c14 = c13 * qz;
            c15 = (2) * qx;
            c16 = c15 * qy;
            c17 = c13 * qy;
            c18 = c15 * qz;
            c19 = c7 - c8;
            c20 = c13 * qx;
            c21 = (2) * qy * qz;
            J = [0, 0, 0, c0, c2, c3;
                0, 0, 0, c1, c0, c5;
                0, 0, 0, c6, c4, c0;
                0, 0, 0, c5, c6, c2;
                c10 + c12 + c7 + c8, -c14 + c16, c17 + c18, 0, 0, 0;
                c14 + c16, c12 + c19 + c9, -c20 + c21, 0, 0, 0;
                -c17 + c18, c20 + c21, c10 + c11 + c19, 0, 0, 0];
        end

        function retval = inverse(self)
            %  Returns group inverse.
            invR = self.so3.inverse();
            retval = Se3(invR, invR*(self.translation * -1));
        end

        function upsilon_omega = log(self)
            %  Logarithmic map
            % 
            %  Computes the logarithm, the inverse of the group exponential which maps
            %  element of the group (rigid body transformations) to elements of the
            %  tangent space (twist).
            % 
            %  To be specific, this function computes ``vee(logmat(.))`` with
            %  ``logmat(.)`` being the matrix logarithm and ``vee(.)`` the vee-operator
            %  of SE(3).
            % 
            %  For the derivation of the logarithm of SE(3), see
            %  J. Gallier, D. Xu, "Computing exponentials of skew symmetric matrices
            %  and logarithms of orthogonal matrices", IJRA 2002.
            %  https://pdfs.semanticscholar.org/cfe3/e4b39de63c8cabd89bf3feff7f5449fc981d.pdf
            %  (Sec. 6., pp. 8)
            upsilon_omega = zeros(6,1);
            omega_and_theta = self.so3.logAndTheta();
            theta = omega_and_theta.theta;
            upsilon_omega(4:6) = omega_and_theta.tangent;
            Omega = So3.hat(upsilon_omega(4:6));

            if (abs(theta) < Constants.epsilon)
                V_inv = eye(3) - (0.5) * Omega + (1. / 12.) * (Omega * Omega);
                upsilon_omega(1:3) = V_inv * self.translation;
            else
                half_theta = (0.5) * theta;
                V_inv = (eye(3) - (0.5) * Omega + ...
                    ((1) - theta * cos(half_theta) / ((2) * sin(half_theta))) / (theta * theta) * (Omega * Omega));
                upsilon_omega(1:3) = V_inv * self.translation;
            end
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
            retval = sprintf("Se3: (%3g, %3g, %3g, %s)", self.translation(1), ...
                self.translation(2), self.translation(3), self.so3.print());
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
            if (isequal(class(other),'Se3') == true)
                retval = Se3(self.so3 * other.so3, self.getTranslation() + self.so3 * other.getTranslation());
            elseif (isscalar(other) == true)
                %retval = Se3( self.so3 * other, self.getTranslation() * other);
                retval = Se3.exp( self.log() * other);
            elseif numel(other) == 3
                %  Group action on 3-points.
                %
                %  This function rotates and translates a three dimensional point ``p`` by the
                %  SE(3) element ``bar_T_foo = (bar_R_foo, t_bar)`` (= rigid body
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
                fprintf(1,"Could not multiply Se3 object.");
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

        function retval = unit_quaternion(self)
            % Accessor of unit quaternion.
            %
            retval = self.so3.unit_quaternion();
        end
    end
end
