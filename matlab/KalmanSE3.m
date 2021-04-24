classdef KalmanSE3 < handle
    % Extended Kalman filter for 3D pose using Lie Group SE(3)

    properties
        ekf_SE3
        state_SE3
        state_SE3_dot
    end
    methods (Static)
        function trajectory = demo()
            deg_per_sec = 95*pi/180;
            %initialState = [0, 0, 0, 0, 0, 0, ...
            %    0, 0, 1, (deg_per_sec/3)^1/3, (deg_per_sec/3)^(1/3), (deg_per_sec/3)^(1/3)]';
            initialState = [0, 0, 0, 0, 0, 0, ...
                0, 1, 5, deg_per_sec, 0, 0]';
            ekf_se3 = KalmanSE3(initialState);
            %gm = importGeometry('BracketTwoHoles.stl');
            %pdegplot(gm);
            box = collisionBox(3,1,2);            
            x_k0 = ekf_se3.getState();            
            box.Pose = x_k0(1).matrix();
            hold off;
            show(box);
            xlim([-10 10])
            ylim([-10 10])
            zlim([-10 10])
            title('Box');
            hold on;
            drawnow;
            start_time = 0;
            end_time = 2;
            frames_per_second = 10;
            dt = 1/frames_per_second;
            trajectory = x_k0(1);
            for t=start_time:dt:(end_time-dt)
                x_k0 = ekf_se3.getState();
                x_k1 = ekf_se3.f(x_k0, 0, dt);
                ekf_se3.setState(x_k1);
                transform = x_k1(1).matrix();
                %det(transform(1:3,1:3))
                box.Pose = transform;                
                show(box);
                title('Box');
                pause(dt);
                trajectory = [trajectory, x_k1(1)];
            end
            show(box);
            title('Box');
            pause(0.1);
            
        end
    end
    methods
        function self = KalmanSE3(varargin)
            if nargin == 0
                self.state_SE3 = Se3(zeros(6,1));
                self.state_SE3_dot = Se3(zeros(6,1));
            elseif nargin == 1 && size(varargin{1},1) == 12
                self.state_SE3 = Se3(varargin{1}(1:6));
                self.state_SE3_dot = Se3(varargin{1}(7:12));
            elseif nargin == 2 && size(varargin{1},1) == 6 && size(varargin{2},1) == 6
                self.state_SE3 = [varargin{1}; varargin{2}];
            else
                fprintf(1,"Could not initialize the KalmanSE3 filter.\n");
                return;
            end
            self.state_SE3.matrix()
            self.state_SE3_dot.matrix()
            %self.f(self.getState(), 0, 1);
            %trackingEKF()
        end

        function x = getState(self)
            x = [self.state_SE3; self.state_SE3_dot];
        end

        function setState(self, x)
            self.state_SE3 = x(1);
            self.state_SE3_dot = x(2);
        end

        function x_predicted = f(self, x, u, dt)
            %! Predict given current pose and velocity
            % Constant velocity model: update position based on last known velocity
            x_predicted(1) = self.state_SE3 * (self.state_SE3_dot * dt);
            x_predicted(2) = self.state_SE3_dot
            %x_predicted(1:6) = Se3.log(Se3.exp(x(1)) * Se3.exp(x(2) * dt));
            % Update velocity directly based on measurements
            %Se3.hat(x(1).log())
            %Se3.hat(x_predicted.log())
        end
        % // new interface
        % void predict(const S& x, const double dt, S& out)
        % {
        % std::cout << "SystemModel::predict" << std::endl;
        % out.x = S::SE3::log(S::SE3::exp(dt * x.v) * S::SE3::exp(x.x));
        % out.v = x.v;
        % }

        function res = plus(self, other)
            res = self.x + other.x;
        end
    end
end