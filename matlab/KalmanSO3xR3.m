classdef KalmanSO3xR3 < handle
    % Extended Kalman filter for 3D pose using Lie Group SE(3)

    properties
        ekf_SO3xR3
        state_SO3xR3
        state_SO3xR3_dot
    end
    methods (Static)
        function trajectory = kalman_lie()
            initialState = [0, 0, 0, 0, 0, 0, ...
               .05, 0.1, 0, 0, 0.1, 0]';           
            ekf_so3xR3 = KalmanSO3xR3(initialState);
            x = ekf_so3xR3.getState();
            %x_mat = x(1).matrix()';
            %xdot_mat = x(2).matrix()';
            % // Generate a trajectory
            traj(1:20) = So3xR3();
            v_init = initialState(7:12);
            traj(1) = So3xR3(initialState(1:6));
            for i=2:size(traj,2)
                traj(i) = So3xR3.exp(traj(i-1).log()) * So3xR3.exp(v_init);
            end
            x_mat_str = sprintf('%0.6g ',x(1).log()');
            xdot_mat_str = sprintf('%0.6g ',x(2).pseudo_log()');
            fprintf(1,"Initial State: %s, %s\n", x_mat_str, xdot_mat_str);
            fd = fopen('output.txt','w+');
            fprintf(fd,"#x_pred;x_mes;x_ekf\n");
            fclose(fd);
            box = collisionBox(3,1,2);            
            x_k0 = ekf_so3xR3.getState();            
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
                x_k0 = ekf_so3xR3.getState();
                x_k1 = ekf_so3xR3.f(x_k0, 0, dt);
                ekf_so3xR3.setState(x_k1);
                transform = x_k1(1).matrix()
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
        
        function trajectory = demo()
            deg_per_sec = 45*pi/180;
            initialState = [3, 3, 3, 0, 0, 0, ...
               3, 0, 0, 0, deg_per_sec, deg_per_sec]';
            %pose = Se3([1, 1, 1, 0, 0, 0]);
            %delta_pose = Se3([0, 0, 0, deg_per_sec, 0, 0]);
            %rotate_around_pose = Se3([pose.getTranslation() - delta_pose.getSo3() * pose.getTranslation(); deg_per_sec; 0; 0]);
            %new_pose = rotate_around_pose * rotate_around_pose * pose;
            %new_pose.matrix()
            ekf_so3xR3 = KalmanSO3xR3(initialState);
            %gm = importGeometry('BracketTwoHoles.stl');
            %pdegplot(gm);
            box = collisionBox(3,1,2);            
            x_k0 = ekf_so3xR3.getState();            
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
            end_time = 1;
            frames_per_second = 10;
            dt = 1/frames_per_second;
            trajectory = x_k0(1);
            for t=start_time:dt:(end_time-dt)
                x_k0 = ekf_so3xR3.getState();
                x_k1 = ekf_so3xR3.f(x_k0, 0, dt);
                ekf_so3xR3.setState(x_k1);
                transform = x_k1(1).matrix()
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
        function self = KalmanSO3xR3(varargin)
            if nargin == 0
                self.state_SO3xR3 = So3xR3(zeros(6,1));
                self.state_SO3xR3_dot = So3xR3(zeros(6,1));
            elseif nargin == 1 && size(varargin{1},1) == 12
                self.state_SO3xR3 = So3xR3(varargin{1}(1:6));
                self.state_SO3xR3_dot = So3xR3(varargin{1}(7:12));
            elseif nargin == 2 && size(varargin{1},1) == 6 && size(varargin{2},1) == 6
                self.state_SO3xR3 = [varargin{1}; varargin{2}];
            else
                fprintf(1,"Could not initialize the KalmanSE3 filter.\n");
                return;
            end
            self.state_SO3xR3.matrix()
            self.state_SO3xR3_dot.matrix()
            %self.f(self.getState(), 0, 1);
            %trackingEKF()
        end

        function x = getState(self)
            x = [self.state_SO3xR3; self.state_SO3xR3_dot];
        end

        function setState(self, x)
            self.state_SO3xR3 = x(1);
            self.state_SO3xR3_dot = x(2);
        end

        function J = computeJacobian(self)
            J_se3 = self.state_SO3xR3.Dx_this_mul_exp_x_at_0();
            J_se3_dot = zeros(7,6);
            J = [J_se3; J_se3_dot];
        end
        
        function x_predicted = f(self, x, u, dt)
            %! Predict given current pose and velocity
            % Constant velocity model: update position based on last known velocity
            pose = self.state_SO3xR3;
            % compute delta pose over time frame dt
            delta_pose = self.state_SO3xR3_dot * dt;
            % decouple rotation and translation components by removing the
            % motion of the frame origin due to the rotational component
            % that effects the translation. The results is that rotational forces do not impact positional
            % updates.
            newTrans = pose.getTranslation() - delta_pose.getSo3() * pose.getTranslation() + delta_pose.getTranslation();
            delta_pose.setTranslation(newTrans);
            % left composition of transformations
            x_predicted(1) = delta_pose * self.state_SO3xR3;
            % constant velocity assumption
            x_predicted(2) = self.state_SO3xR3_dot;
        end
        % // new interface
        % void predict(const S& x, const double dt, S& out)
        % {
        % std::cout << "SystemModel::predict" << std::endl;
        % out.x = S::SE3::log(S::SE3::exp(dt * x.v) * S::SE3::exp(x.x));
        % out.v = x.v;
        % }
    end
end