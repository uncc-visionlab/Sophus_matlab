classdef So3 < RotationMatrix
    %  SO3 base type - implements SO3 class but is storage agnostic.
    %
    %  SO(3) is the group of rotations in 3d. As a matrix group, it is the set of
    %  matrices which are orthogonal such that ``R * R' = I`` (with ``R'`` being the
    %  transpose of ``R``) and have a positive determinant. In particular, the
    %  determinant is 1. Internally, the group is represented as a unit quaternion.
    %  Unit quaternion can be seen as members of the special unitary group SU(2).
    %  SU(2) is a double cover of SO(3). Hence, for every rotation matrix ``R``,
    %  there exists two unit quaternion: ``(r, v)`` and ``(r, -v)``, with ``r`` the
    %  real part and ``v`` being the imaginary 3-vector part of the quaternion.
    %
    %  SO(3) is a compact, but non-commutative group. First it is compact since the
    %  set of rotation matrices is a closed and bounded set. Second it is
    %  non-commutative since the equation ``R_1 * R_2 = R_2 * R_1`` does not hold in
    %  general. For example rotating an object by some degrees about its ``x``-axis
    %  and then by some degrees about its y axis, does not lead to the same
    %  orienation when rotation first about ``y`` and then about ``x``.
    %
    %  Class invariant: The 2-norm of ``unit_quaternion`` must be close to 1.
    %  Technically speaking, it must hold that:
    %
    %    ``|unit_quaternion().squaredNorm() - 1| <= Constants<Scalar>::epsilon()``.

    properties (Constant)
        % Degrees of freedom of group, number of dimensions in tangent space.
        DoF = 3;
        % Number of internal parameters used (quaternion is a 4-tuple).
        num_parameters = 4;
        % Group transformations are 3x3 matrices.
        N = 3;
    end

    properties
        unit_quaternion_
    end

    methods (Static)
        %////////////////////////////////////////////////////////////////////////////
        %// public static functions
        %////////////////////////////////////////////////////////////////////////////


        function J = Dx_exp_x(omega)
            %  Returns derivative of exp(x) wrt. x.
            %
            c0 = omega(1) * omega(1);
            c1 = omega(2) * omega(2);
            c2 = omega(3) * omega(3);
            c3 = c0 + c1 + c2;

            if (c3 < Constants.epsilon)
                J = So3.Dx_exp_x_at_0();
                return;
            end

            c4 = sqrt(c3);
            c5 = 1.0 / c4;
            c6 = 0.5 * c4;
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
            c20 = (0.5) * c5 * c7;
            J(1, 1) = -c0 * c10 + c0 * c13 + c8;
            J(1, 2) = c16;
            J(1, 3) = c17;
            J(2, 1) = c16;
            J(2, 2) = -c1 * c10 + c1 * c13 + c8;
            J(2, 3) = c19;
            J(3, 1) = c17;
            J(3, 2) = c19;
            J(3, 3) = -c10 * c2 + c13 * c2 + c8;
            J(4, 1) = -c20 * omega(1);
            J(4, 2) = -c20 * omega(2);
            J(4, 3) = -c20 * omega(3);
        end
        function J = Dx_exp_x_at_0()
            %  Returns derivative of exp(x) wrt. x_i at x=0.
            %
            %  clang-format off
            J = [(0.5), (0), (0);
                (0), (0.5), (0);
                (0), (0), (0.5);
                (0), (0), (0)];
            %  clang-format on
        end

        function retval = Dxi_exp_x_matrix_at_0(i)
            %  Returns derivative of exp(x).matrix() wrt. ``x_i at x=0``.
            %
            retval = generator(i);
        end

        function retval = exp(omega)
            %  Group exponential
            %
            %  This functions takes in an element of tangent space (= rotation vector
            %  ``omega``) and returns the corresponding element of the group SO(3).
            %
            %  To be more specific, this function computes ``expmat(hat(omega))``
            %  with ``expmat(.)`` being the matrix exponential and ``hat(.)`` being the
            %  hat()-operator of SO(3).
            %
            [retval, ~] = So3.expAndTheta(omega);
        end

        function [so3_out, theta] = expAndTheta(omega)
            %  As above, but also returns ``theta = |omega|`` as out-parameter.
            %
            %  Precondition: ``theta`` must not be ``nullptr``.
            %
            theta_sq = sum(omega.^2);
            if (theta_sq < Constants.epsilon * Constants.epsilon)
                theta = (0);
                theta_po4 = theta_sq * theta_sq;
                imag_factor = (0.5) - (1.0 / 48.0) * theta_sq + ...
                    (1.0 / 3840.0) * theta_po4;
                real_factor = (1) - (1.0 / 8.0) * theta_sq + ...
                    (1.0 / 384.0) * theta_po4;
            else
                theta = sqrt(theta_sq);
                half_theta = (0.5) * (theta);
                sin_half_theta = sin(half_theta);
                imag_factor = sin_half_theta / (theta);
                real_factor = cos(half_theta);
            end
            so3_out = So3(quaternion(real_factor, imag_factor * omega(1), ...
                imag_factor * omega(2), imag_factor * omega(3)));
            if ((sum(so3_out.unit_quaternion().compact().^2) - (1)) >= Constants.epsilon)
                fprintf(1,"SO3::exp failed! omega: (%f,%f,%f), real: %f, img: %f\n", ...
                    omega', real_factor, imag_factor);
            end
        end
        function retval = fitToSO3(R)
            %  Returns closest SO3 given arbitrary 3x3 matrix.
            %
            retval = So3.makeRotationMatrix(R);
        end
        function retval = generator(i)
            %  Returns the ith infinitesimal generators of SO(3).
            %
            %  The infinitesimal generators of SO(3) are:
            %
            %  ```
            %          |  0  0  0 |
            %    G_0 = |  0  0 -1 |
            %          |  0  1  0 |
            %
            %          |  0  0  1 |
            %    G_1 = |  0  0  0 |
            %          | -1  0  0 |
            %
            %          |  0 -1  0 |
            %    G_2 = |  1  0  0 |
            %          |  0  0  0 |
            %  ```
            %
            %  Precondition: ``i`` must be 0, 1 or 2.
            %
            if (i < 0 | i > 2)
                fprintf(1,"i should be an integer in range [0,2].\n");
            end
            e = zeros(1,N);
            e(i+1) = 1;
            retval = hat(e);
        end
        function retval =  hat(omega)
            %  hat-operator
            %
            %  It takes in the 3-vector representation ``omega`` (= rotation vector) and
            %  returns the corresponding matrix representation of Lie algebra element.
            %
            %  Formally, the hat()-operator of SO(3) is defined as
            %
            %    ``hat(.): R^3 -> R^{3x3},  hat(omega) = sum_i omega_i * G_i``
            %    (for i=0,1,2)
            %
            %  with ``G_i`` being the ith infinitesimal generator of SO(3).
            %
            %  The corresponding inverse is the vee()-operator, see below.
            %
            retval = [0, -omega(3), omega(2); omega(3), 0, -omega(1); -omega(2), omega(1) 0];
        end
        function retval = lieBracket(omega1, omega2)
            %  Lie bracket
            %
            %  It computes the Lie bracket of SO(3). To be more specific, it computes
            %
            %    ``[omega_1, omega_2]_so3 := vee([hat(omega_1), hat(omega_2)])``
            %
            %  with ``[A,B] := AB-BA`` being the matrix commutator, ``hat(.)`` the
            %  hat()-operator and ``vee(.)`` the vee()-operator of SO3.
            %
            %  For the Lie algebra so3, the Lie bracket is simply the cross product:
            %
            %  ``[omega_1, omega_2]_so3 = omega_1 x omega_2.``
            %
            retval = cross(omega1, omega2);
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
        function retval = sampleUniform()
            %  Draw uniform sample from SO(3) manifold.
            %  Based on: http://planning.cs.uiuc.edu/node198.html
            %
            u1 = rand(1);
            u2 = rand(1)*2*Constants.M_PI;
            u3 = rand(1)*2*Constants.M_PI;
            a = sqrt(1 - u1);
            b = sqrt(u1);
            retval = So3(quaternion(b * cos(u3), a * sin(u2), a * cos(u2), b * sin(u3)));
        end
        function retval = vee(omega)
            %  vee-operator
            %
            %  It takes the 3x3-matrix representation ``Omega`` and maps it to the
            %  corresponding vector representation of Lie algebra.
            %
            %  This is the inverse of the hat()-operator, see above.
            %
            %  Precondition: ``Omega`` must have the following structure:
            %
            %                 |  0 -c  b |
            %                 |  c  0 -a |
            %                 | -b  a  0 |
            %
            retval = [omega(3, 2); omega(1,3); omega(2,1)];
        end
    end

    methods
        function self = So3(varargin)
            % internally represented by a unit complex number z
            if (nargin == 1 && isequal(class(varargin{1}),'So3'))
                self = So3(varargin{1}.unit_quaternion());
            elseif (nargin == 1 && isequal(class(varargin{1}),'quaternion'))
                % input is a quaternion
                self.unit_quaternion_ = varargin{1};
                %elseif (nargin == 1 && numel(varargin{1}) == 1 && imag(varargin{1}) > 0)
                %    % input is a complex number
                %    self.z = complex(varargin{1});
            elseif (nargin == 1 && numel(varargin{1}) == 3)
                self = So3.exp([varargin{1}(1); varargin{1}(2); varargin{1}(3)]);
            elseif (nargin == 1 && numel(varargin{1}) == 9 && size(varargin{1},1) == 3 && size(varargin{1},2) == 3)
                % input is a 3x3 rotation matrix
                self.unit_quaternion_ = quaternion(rotm2quat(varargin{1}));
            elseif false
                % TODO:
                % SOPHUS_ENSURE(isOrthogonal(R), "R is not orthogonal:\n %",
                % R * R.transpose());
                % SOPHUS_ENSURE(R.determinant() > Scalar(0), "det(R) is not positive: %",
                % R.determinant());
                %    self.real = complex(0.5*(varargin{1}(1,1) + varargin{1}(2,2)) + ...
                %        1i*(0.5*(varargin{1}(2,1) + varargin{1}(1,2))));
                %    if (abs(det(varargin{1}) - 1) <= Constants.epsilon)
                %        fprintf(1,"det(R) should be (close to) 1.\n")
                %    end
                %elseif (nargin == 2)
                %    % input is a pair of arguments arg1 + i*arg2
                %    self.z = complex(varargin{1} + 1i*varargin{2});
            elseif nargin == 0
                % empty initializer creates the identity rotation in So3
                self.unit_quaternion_ = quaternion(1, 0, 0, 0);
            else                
                fprintf(1,"Error could not construct So3 object.\n");
            end
        end

        function retval = Adj(self)
            %  Adjoint transformation
            %
            %  This function return the adjoint transformation ``Ad`` of the group
            %  element ``A`` such that for all ``x`` it holds that
            %  ``hat(Ad_A * x) = A * hat(x) A^{-1}``. See hat-operator below.
            %
            %  For SO(3), it simply returns the rotation matrix corresponding to ``A``.
            %
            retval = self.matrix();
        end

        function retval = angleX(self)
            %  Extract rotation angle about canonical X-axis
            %
            R = self.matrix();
            Rx = R(2:3,2:3);
            retval = So2(RotationMatrix.makeRotationMatrix(Rx)).log();
        end

        function retval = angleY(self)
            %  Extract rotation angle about canonical Y-axis
            %
            R = self.matrix();
            Ry = [R(1,1), R(3,1); R(1,3), R(3,3)];
            retval = So2(RotationMatrix.makeRotationMatrix(Ry)).log();
        end

        function retval = angleZ(self)
            %  Extract rotation angle about canonical X-axis
            %
            R = self.matrix();
            Rz = R(1:2,1:2);
            retval = So2(RotationMatrix.makeRotationMatrix(Rz)).log();
        end

        function J = Dx_this_mul_exp_x_at_0(self)
            %  Returns derivative of  this * SO3::exp(x)  wrt. x at x=0.
            %
            [qw, qx, qy, qz] = self.unit_quaternion().parts();
            c0 = (0.5) * qw;
            c1 = (0.5) * qz;
            c2 = -c1;
            c3 = (0.5) * qy;
            c4 = (0.5) * qx;
            c5 = -c4;
            c6 = -c3;
            J = [c0, c2, c3; c1, c0 c5; c6, c4, c0; c5 c6 c2];
        end

        function retval = params(self)
            %  Returns internal parameters of SO(3).
            %
            %  It returns (q.real, q.imag[0], q.imag[1], q.imag[2]), with q being the
            %  unit quaternion.
            %
            retval = self.unit_quaternion().compact()';
        end

        function retval = inverse(self)
            %  Returns group inverse.
            %
            retval = So3(self.unit_quaternion_.conj());
        end

        function retval = log(self)
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
            retval = self.logAndTheta().tangent;
        end

        function J = logAndTheta(self)
            %  As above, but also returns ``theta = |omega|``.
            %
            [qw, qx, qy, qz] = self.unit_quaternion().parts();
            squared_n = sum([qx qy qz].^2);

            %  Atan-based log thanks to
            %
            %  C. Hertzberg et al.:
            %  "Integrating Generic Sensor Fusion Algorithms with Sound State
            %  Representation through Encapsulation of Manifolds"
            %  Information Fusion, 2011

            if (squared_n < Constants.epsilon * Constants.epsilon)
                %  If quaternion is normalized and n=0, then w should be 1;
                %  w=0 should never happen here!
                if (abs(qw) < Constants.epsilon)
                    fprintf(1,"Quaternion (%f,%fi,%fj,%fk) should be normalized!", ...
                        self.unit_quaternion().compact());
                end
                squared_w = qw * qw;
                two_atan_nbyw_by_n = ...
                    (2) / qw - (2.0/3.0) * (squared_n) / (qw * squared_w);
                J.theta = (2) * squared_n / qw;
            else
                n = sqrt(squared_n);
                if (abs(qw) < Constants.epsilon)
                    if (qw > (0))
                        two_atan_nbyw_by_n = Constants.pi() / n;
                    else
                        two_atan_nbyw_by_n = -Constants.pi() / n;
                    end
                else
                    two_atan_nbyw_by_n = (2) * atan(n / qw) / n;
                end
                J.theta = two_atan_nbyw_by_n * n;
            end
            J.tangent = two_atan_nbyw_by_n * [qx; qy; qz];
        end

        function normalize(self)
            %  It re-normalizes ``unit_quaternion`` to unit length.
            %
            %  Note: Because of the class invariant, there is typically no need to call
            %  this function directly.
            %
            length = norm(self.unit_quaternion());
            if (length < Constants.epsilon)
                fprintf(1,"Quaternion (%f,%fi,%fj,%fk) should not be close to zero!\n",self.unit_quaternion_.compact());
            end
            self.unit_quaternion_ = self.unit_quaternion_/length;
        end

        function retval = toString(self)
            retval = sprintf("So3: (%f,%fi,%fj,%fk)", self.unit_quaternion_.compact())
        end

        function retval = matrix(self)
            %  Returns 3x3 matrix representation of the instance.
            %
            %  For SO(3), the matrix representation is an orthogonal matrix ``R`` with
            %  ``det(R)=1``, thus the so-called "rotation matrix".
            %
            retval = self.unit_quaternion_.rotmat('point');
        end

        function retval = mtimes(self, other)
            % overload function for the * operator in MATLAB
            if (isequal(class(other),'So3') == true)
                % Group multiplication, which is rotation concatenation.
                % We can assume that the squared-norm is close to 1 since we deal with a
                % unit complex number. Due to numerical precision issues, there might
                % be a small drift after pose concatenation. Hence, we need to renormalizes
                % the complex number here.
                % Since squared-norm is close to 1, we do not need to calculate the costly
                % square-root, but can use an approximation around 1 (see
                % http://stackoverflow.com/a/12934750 for details).
                %[aw, ax, ay, az] = self.unit_quaternion().parts();
                %[bw, bx, by, bz] = other.unit_quaternion().parts();
                %quat = quaternion(aw * bw - ax * bx - ay * by - az * bz, ...
                %    aw * bx + ax * bw + ay * bz - az * by, ...
                %    aw * by + ay * bw + az * bx - ax * bz, ...
                %    aw * bz + az * bw + ax * by - ay * bx);
                %retval = So3( quat);
                retval = So3( self.unit_quaternion() * other.unit_quaternion());
            elseif (isscalar(other) == true)
                retval = So3.exp(self.log() * other);
            elseif nargin == 1 && numel(other) == 9 && size(other,1) == 3 && size(other,2) == 3
                retval = self * So3(other);
            elseif numel(other) == 3
                %  Group action on 3-points.
                %
                %  This function rotates a 3 dimensional point ``p`` by the SO3 element
                %   ``bar_R_foo`` (= rotation matrix): ``p_bar = bar_R_foo * p_foo``.
                %
                %  Since SO3 is internally represented by a unit quaternion ``q``, it is
                %  implemented as ``p_bar = q * p_foo * q^{*}``
                %  with ``q^{*}`` being the quaternion conjugate of ``q``.
                %
                %  Geometrically, ``p``  is rotated by angle ``|omega|`` around the
                %  axis ``omega/|omega|`` with ``omega := vee(log(bar_R_foo))``.
                %
                %  For ``vee``-operator, see below.
                %
                [qw, qx, qy, qz] = self.unit_quaternion().parts();
                qvec = [qx; qy; qz];
                uv = cross(qvec, other);
                uv = uv + uv;
                retval = other + qw * uv + cross(qvec, uv);
                % TODO:
                %  Group action on homogeneous 3-points. See above for more details.
                % template <typename HPointDerived,
                %         typename = typename std::enable_if<
                %             IsFixedSizeVector<HPointDerived, 4>::value>::type>
                % SOPHUS_FUNC HomogeneousPointProduct<HPointDerived> operator*(
                %   Eigen::MatrixBase<HPointDerived> const& p) const {
                % const auto rp = *this * p.template head<3>();
                % return HomogeneousPointProduct<HPointDerived>(rp(0), rp(1), rp(2), p(3));
                % }
                %
                %  Group action on lines.
                %
                %  This function rotates a parametrized line ``l(t) = o + t * d`` by the SO3
                %  element:
                %
                %  Both direction ``d`` and origin ``o`` are rotated as a 3 dimensional point
                %
                % SOPHUS_FUNC Line operator*(Line const& l) const {
                % return Line((*this) * l.origin(), (*this) * l.direction());
                % }
            else
                fprintf(1,"Could not multiply So3 object.\n");
            end
        end

        function setQuaternion(self, quaternion_)
            %  Takes in complex number / tuple and normalizes it.
            %
            %  Precondition: The complex number must not be close to zero.
            %
            self.unit_quaternion_ = quaternion_;
            self.normalize();
        end

        function retval = unit_quaternion(self)
            % Accessor of unit quaternion.
            %
            retval = self.unit_quaternion_;
        end
    end
end
