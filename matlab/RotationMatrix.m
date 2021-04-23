classdef RotationMatrix < Constants
    %/// Rotation matrix helper functions.

    methods (Static)
        function retval = isOrthogonal(R)
            % /// Takes in arbitrary square matrix and returns true if it is
            % /// orthogonal.
            [N, M] = size(sR);
            % TODO:
            % static_assert(N == M, "must be a square matrix");
            % static_assert(N >= 2, "must have compile time dimension >= 2");
            retval = norm(R * R' - eye(N)) < Constants.epsilon;
        end
        function retval = isScaledOrthogonalAndPositive(sR)
            % /// Takes in arbitrary square matrix and returns true if it is
            % /// "scaled-orthogonal" with positive determinant.
            % ///
            [N, M] = size(sR);
            det = det(sR);
            if (det <=0)
                retval = false;
            end
            scale_sqr = det^(2. / N);
            % TODO:
            % static_assert(N == M, "must be a square matrix");
            % static_assert(N >= 2, "must have compile time dimension >= 2");
            retval = norm((sR * sR' - scale_sqr * eye(N)),Inf) < ...
                sqrt(Constants.epsilon);
        end

        function retval = makeRotationMatrix(R)
            % /// Takes in arbitrary square matrix (2x2 or larger) and returns closest
            % /// orthogonal matrix with positive determinant.
            [N, M] = size(sR);
            % TODO:
            %  static_assert(N == M, "must be a square matrix");
            %  static_assert(N >= 2, "must have compile time dimension >= 2");
            [U,S,V] = svd(R);
            %  Determine determinant of orthogonal matrix U*V'.
            d = det(U * V');
            %  Starting from the identity matrix D, set the last entry to d (+1 or
            %  -1),  so that det(U*D*V') = 1.
            Diag = eye(N);
            Diag(N, N) = d;
            retval = U * Diag * V';
        end
    end
end