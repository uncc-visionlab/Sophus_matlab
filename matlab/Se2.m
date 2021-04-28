classdef Se2 < Constants
    %  SE2 base type - implements SE2 class but is storage agnostic.
    %
    %  SE(2) is the group of rotations  and translation in 2d. It is the semi-direct
    %  product of SO(2) and the 2d Euclidean vector space.  The class is represented
    %  using a composition of SO2Group  for rotation and a 2-vector for translation.
    %
    %  SE(2) is neither compact, nor a commutative group.
    %
    %  See SO2Group for more details of the rotation representation in 2d.
    %
    properties (Constant)
        % Degrees of freedom of manifold, number of dimensions in tangent space
        % (two for translation, three for rotation).
        DoF = 3;
        % Number of internal parameters used (tuple for complex, two for
        % translation).
        num_parameters = 4;
        % Group transformations are 3x3 matrices.
        N = 3;
    end

    properties
        so2
        translation
    end

    methods (Static)
        % ////////////////////////////////////////////////////////////////////////////
        % // public static functions
        % ////////////////////////////////////////////////////////////////////////////

        function res = d_lieBracketab_by_d_a(b)
            %  Derivative of Lie bracket with respect to first element.
            %
            %  This function returns ``D_a [a, b]`` with ``D_a`` being the
            %  differential operator with respect to ``a``, ``[a, b]`` being the lie
            %  bracket of the Lie algebra se3.
            %  See ``lieBracket()`` below.
            %
            zero = 0;
            upsilon2 = b(1:2);
            theta2 = b(3);
            res = [zero, theta2, -upsilon2(2); -theta2, zero, upsilon2(1); zero, zero, zero];
        end

        function retval =  exp(a)
            %  Group exponential
            %
            %  This functions takes in an element of tangent space (= twist ``a``) and
            %  returns the corresponding element of the group SE(2).
            %
            %  The first two components of ``a`` represent the translational part
            %  ``upsilon`` in the tangent space of SE(2), while the last three components
            %  of ``a`` represents the rotation vector ``omega``.
            %  To be more specific, this function computes ``expmat(hat(a))`` with
            %  ``expmat(.)`` being the matrix exponential and ``hat(.)`` the hat-operator
            %  of SE(2), see below.
            %
            theta = a(3);
            so2 = So2.exp(theta);
            if (abs(theta) < Constants.epsilon)
                theta_sq = theta * theta;
                sin_theta_by_theta = 1. - (1. / 6.) * theta_sq;
                one_minus_cos_theta_by_theta = 0.5 * theta - (1. / 24.) * theta * theta_sq;
            else
                sin_theta_by_theta = imag(so2.unit_complex()) / theta;
                one_minus_cos_theta_by_theta = (1.) - real(so2.unit_complex()) / theta;
            end
            trans = [sin_theta_by_theta * a(1) - one_minus_cos_theta_by_theta * a(2), ...
                one_minus_cos_theta_by_theta * a(1) + sin_theta_by_theta * a(2)];
            retval = Se2( so2, trans);
        end

        function retval = fitToSE2(T)
            % /// Returns closest SE2 given arbitrary 3x3 matrix.
            % ///
            retval = Se2( So2.fitToSO2(T(1:2,1:2)), T(1:2,3));
        end

        function retval = generator(i)
            %  Returns the ith infinitesimal generators of SE(2).
            %
            %  The infinitesimal generators of SE(2) are:
            %
            %          |  0  0  1 |
            %    G_0 = |  0  0  0 |
            %          |  0  0  0 |
            %
            %          |  0  0  0 |
            %    G_1 = |  0  0  1 |
            %          |  0  0  0 |
            %
            %          |  0 -1  0 |
            %    G_2 = |  1  0  0 |
            %          |  0  0  0 |
            %  Precondition: ``i`` must be in 0, 1 or 2.
            %
            if (i < 0 || i > 2)
                fprintf(1,"i should be an integer in range [0,2].\n");
            end
            e = zeros(1,N);
            e(i+1) = 1;
            retval = Se2.hat(e);
        end

        function retval =  hat(a)
            %  hat-operator
            %
            %  It takes in the 3-vector representation (= twist) and returns the
            %  corresponding matrix representation of Lie algebra element.
            %
            %  Formally, the ``hat()`` operator of SE(3) is defined as
            %
            %    ``hat(.): R^3 -> R^{3x33},  hat(a) = sum_i a_i * G_i``  (for i=0,1,2)
            %
            %  with ``G_i`` being the ith infinitesimal generator of SE(2).
            %
            retval = [So2.hat(a(3)); [a(1); a(2)]; 0 0 0];
        end
        function retval = lieBracket(a,b)
            %  Lie bracket
            %
            %  It computes the Lie bracket of SE(2). To be more specific, it computes
            %
            %    ``[omega_1, omega_2]_se2 := vee([hat(omega_1), hat(omega_2)])``
            %
            %  with ``[A,B] := AB-BA`` being the matrix commutator, ``hat(.) the
            %  hat-operator and ``vee(.)`` the vee-operator of SE(2).
            %
            upsilon1 = a(1:2);
            upsilon2 = b(1:2);
            theta1 = a(3);
            theta2 = b(3);

            retval = [-theta1 * upsilon2(2) + theta2 * upsilon1(2), ...
                theta1 * upsilon2(1) - theta2 * upsilon1(1), (0)];
        end
        function upsilon_theta = log(other)
            %  Logarithmic map
            %
            %  Computes the logarithm, the inverse of the group exponential which maps
            %  element of the group (rigid body transformations) to elements of the
            %  tangent space (twist).
            %
            %  To be specific, this function computes ``vee(logmat(.))`` with
            %  ``logmat(.)`` being the matrix logarithm and ``vee(.)`` the vee-operator
            %  of SE(2).
            %
            upsilon_theta = zeros(1,N);
            so2 = other.getSo2();
            theta = so2.log_();
            upsilon_theta(3) = theta;
            halftheta = (0.5) * theta;

            z = so2.unit_complex();
            real_minus_one = real(z) - 1.;
            if (abs(real_minus_one) < Constants.epsilon)
                halftheta_by_tan_of_halftheta = 1. - (1. / 12) * theta * theta;
            else
                halftheta_by_tan_of_halftheta = -(halftheta * imag(z)) / (real_minus_one);
            end
            V_inv = [halftheta_by_tan_of_halftheta, halftheta; -halftheta, halftheta_by_tan_of_halftheta];
            upsilon_theta(1:2) = V_inv * other.getTranslation();
        end
        function retval = vee(omega)
            %  vee-operator
            %
            %  It takes the 3x3-matrix representation ``Omega`` and maps it to the
            %  corresponding 3-vector representation of Lie algebra.
            %
            %  This is the inverse of the hat-operator, see above.
            %
            %  Precondition: ``Omega`` must have the following structure:
            %
            %                 |  0 -d  a |
            %                 |  d  0  b |
            %                 |  0  0  0 | .
            %
            if (norm(omega(3,:),1) > Constants.epsilon)
                fprintf(1,"Omega: \n% %f", omega);
            end
            retval = [omega(1:2,3); Se2.vee(omega(1:2,1:2))];
        end
    end

    methods
        function self = Se2(varargin)
            % internally represented by a unit complex number z
            if (nargin == 1 && isequal(class(varargin{1}),'Se2'))
                self.translation = varargin{1}.getTranslation();
                self.so2 = varargin{1}.getSo2();
            elseif (nargin == 1 && numel(varargin{1}) == 3)
                % input is an tangent vector
                self.translation = varargin{1}(1:2);
                self.so2 = So2.exp(varargin{1}(3));
            elseif (nargin == 2 && numel(varargin{1}) == 4 && size(varargin{1},1) == 2 && size(varargin{1},2) == 2 && ...
                    numel(varargin{2}) == 2 && size(varargin{2},1) == 2 && size(varargin{2},2) == 1)
                % input is a 2x2 rotation matrix and translation vector
                self.so2 = So2(varargin{1}(1:2,1:2));
                self.translation = [varargin{2}(1); varargin{2}(2)];
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
            res = Se2.eye(N);
            res(1:2,1:2) = self.so2.matrix();
            res(1, 3) = self.translation(2);
            res(2, 3) = -self.translation(1);
        end

        function retval = inverse(self)
            %  Returns group inverse.
            invR = self.so2.inverse();
            retval = Se2(invR, invR*(self.translation * -1));
        end

        function retval = log_(self)
            %  Logarithmic map
            %
            %  Returns tangent space representation (= twist) of the instance.
            %
            retval = Se2.log(self.z);
        end

        function normalize(self)
            % /**
            % * \brief Normalize SO2 element
            % *
            % * It re-normalizes the SO2 element.
            % */
            self.so2.normalize();
        end

        function retval = matrix(self)
            %  Returns 3x3 matrix representation of the instance.
            %
            %  It has the following form:
            %
            %    | R t |
            %    | o 1 |
            %
            %  where ``R`` is a 2x2 rotation matrix, ``t`` a translation 2-vector and
            %  ``o`` a 2-column vector of zeros.
            %
            retval = [self.matrix2x3(); 0 0 1];
        end

        function retval = matrix2x3(self)
            %  Returns the significant first two rows of the matrix above.
            %
            retval = [self.so2.matrix(), self.translation];
        end

        function retval = print(self)
            retval = sprintf("Se2: (%f,%f,%f)",self.translation(1),self.translation(2), ...
                self.so2.log_());
        end

        function retval = mtimes(self, other)
            % overload function for the * operator in MATLAB
            if (isequal(class(other),'Se2') == true)
                retval = Se2(self * other);
            elseif numel(other) == 2
                %  Group action on 2-points.
                %
                %  This function rotates and translates a two dimensional point ``p`` by the
                %  SE(2) element ``bar_T_foo = (bar_R_foo, t_bar)`` (= rigid body
                %  transformation):
                %
                %    ``p_bar = bar_R_foo * p_foo + t_bar``.
                %
                retval = self.so2 * other + self.translation;
            elseif numel(other) == 3
                % Group action on homogeneous 2-points. See above for more details.
                % 
                retval = [self.so2 * other(1:2) + other(3) * self.translation; other(3)];
            elseif false
                % /// Group action on lines.
                % ///
                % /// This function rotates and translates a parametrized line
                % /// ``l(t) = o + t * d`` by the SE(2) element:
                % ///
                % /// Origin ``o`` is rotated and translated using SE(2) action
                % /// Direction ``d`` is rotated using SO(2) action
                % ///
                % TODO:
                % SOPHUS_FUNC Line operator*(Line const& l) const {
                % return Line((*this) * l.origin(), so2() * l.direction());
                % }
            else
                fprintf(1,"Could not multiply Se2 object.");
            end
        end

        function retval = rotationMatrix(self)
            %  Returns rotation matrix.
            %
            retval = self.so2.matrix();
        end

        function setComplex(self, z)
            %  Takes in complex number, and normalizes it.
            %
            %  Precondition: The complex number must not be close to zero.
            %
            self.so2.setComplex(z);
        end

        function setRotationMatrix(self, R)
            %  Sets ``so2`` using ``rotation_matrix``.
            %
            %  Precondition: ``R`` must be orthogonal and ``det(R)=1``.
            %
            self.so2.setComplex(complex(0.5 * (R(1, 1) + R(2, 2)) + 1i * 0.5 * (R(2, 1) - R(1, 2))));
        end

        function retval = getSo2(self)
            %  Accessor of SO2 group.
            %
            retval = self.so2;
        end

        function retval = getTranslation(self)
            %  Accessor of translation vector
            %
            retval = self.translation;
        end

        function retval = unit_complex(self)
            %  Accessor of unit complex number.
            %
            retval = self.so2.unit_complex();
        end
    end
end
