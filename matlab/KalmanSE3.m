classdef KalmanSE3 < handle
    % Extended Kalman filter for 3D pose using Lie Group SE(3)
    
    properties
        ekf_SE3
        state_SE3
        state_SE3_dot
    end
    methods (Static)
        function trajectory = kalman_lie()
            initialState = [0, 0, 0, 0, 0, 0, ...
                .05, 0.1, 0, 0, 0.1, 0]';
            ekf_se3 = KalmanSE3(initialState);
            x = ekf_se3.getState();
            %x_mat = x(1).matrix()';
            %xdot_mat = x(2).matrix()';
            % // Generate a trajectory
            num_steps = 20;
            traj(1:num_steps) = Se3();
            v_init = initialState(7:12);
            traj(1) = Se3(initialState(1:6));
            velSe3 = Se3(v_init);
            traj2 = zeros(6,num_steps);
            vel2 = zeros(6,num_steps);
            traj2(:,1) = initialState(1:6);
            vel2(:,1) = initialState(7:12);
            for i=2:size(traj,2)
                %traj(i) = Se3.exp(traj(i-1).log()) * Se3.exp(v_init);
                traj(i) = traj(i-1) * velSe3
                newpose = Se3.exp(traj2(:,i-1)) * Se3.exp(vel2(:,1));
                traj2(:,i) = newpose.log();
                x_mat_str = sprintf('%0.6g ',traj(i).log()');
                fprintf(1,"%s\n",x_mat_str);
            end
            %return
            fd = fopen('output.txt','w+');
            fprintf(fd,"#x_pred;x_mes;x_ekf\n");
            x_mat_str = sprintf('%0.6g ',x(1).log()');
            xdot_mat_str = sprintf('%0.6g ',[x(2).getTranslation(); x(2).getSo3().log()]');
            fprintf(1,"Initial State: %s, %s\n", x_mat_str, xdot_mat_str);
            x_mat_str = sprintf('%0.6g;',x(1).log()');
            %x_mat_str(end)=';';
            fprintf(fd,"%s%s%s",x_mat_str,x_mat_str,x_mat_str);
            x_mes = Se3([0,0,0,0,0,0]');
            for i=2:size(traj,2)
                x_ref = traj(i);
                
                % Use reference velocity for testing
                % velocity_measurement.addPosition(x_mes);
                % velocity_measurement.addPosition(x_ref);
                
                % Simulate system (constant velocity model)
                % x = sys.f(x, u);
                % std::cout << "New State: " << x.x.transpose() << ", " << x.v.transpose() << std::endl;
                     

                % Predict state for current time-step using the filters
                x_k = ekf_se3.getState();
                x_kplus1 = ekf_se3.predict(x_k, 1.);
                x_kplus1_p =  x_k(2) * x_k(1);
                x_kplus1_pp = (Se3.exp(1. * x(2).log()) * Se3.exp(x(1).log()));
                x_kplus1_ppp = (Se3.exp(1 * vel2(:,1)) * Se3.exp(traj2(:,i-1)));
                x_kplus1_p_vec = x_kplus1_p.log();
                x_kplus1_pp_vec = x_kplus1_pp.log();
                x_kplus1_ppp_vec = x_kplus1_ppp.log();
                ekf_se3.setState(x_kplus1);
                %x_kplus1 = ekf_se3.getState();
                x_k_pos_str = sprintf('%0.6g ',x_k(1).log()');
                x_kplus1_pos_str = sprintf('%0.6g ',x_kplus1(1).log()');
                %x_kplus1_pos_str = sprintf('%0.6g ',x_kplus1_ppp(1).log()');
                x_kplus1_vel_str = sprintf('%0.6g ',x_kplus1(2).log()');
                
                fprintf(1,"Pred State: %s, %s\n", x_kplus1_pos_str, x_kplus1_vel_str);
                
                % Sensor update every 3 iteration, predict the rest of the time
                % if (i % 2 == 0)
                % {
                %     std::cout << "UPDATING VELOCITY" << std::endl;
                %     v_mes = velocity_measurement.v;
                %     ekf.update(velocity_measurement, v_mes);
                % }
                
                if (mod(i,3) == 0)
                    %std::cout << "UPDATING POSITION" << std::endl;
                    %x_mes = noise.addNoise(x_ref);
                    %ekf.update(pos_measurement, x_mes);
                end
                %auto cov = pos_measurement.getCovariance();
                %std::cout << "covariance:\n"
                %          << ekf.getCovariance().matrix() << std::endl;
                x_mes_str = sprintf('%0.6g ',x_mes(1).log()');
                %std::cout << "x_pred: " << x_ref.transpose() << std::endl;
                %std::cout << "x_mes: " << x_mes.transpose() << std::endl;
                %std::cout << "x_ekf: " << x_ekf.x.transpose() << std::endl;
                %csv.write({x_ref, x_mes, x_ekf.x});
                fprintf(fd,"%s%s%s",x_k_pos_str,x_k_pos_str,x_kplus1_pos_str);
            end
            fclose(fd);
            
            %State x_future;
            x_future = ekf_se3.predict(ekf_se3.getState(), 5);
            x_future_pos_str = sprintf('%0.6g ',x_future(1).log()');
            fprintf(1,"Future (dt=5): %s\n", x_future_pos_str);
        end
        
        function trajectory = demo()
            rad_per_sec = 45*pi/180;
            initialState = [3, 3, 3, 0, 0, 0, ...
                3, 0, 0, 0, rad_per_sec, rad_per_sec]';
            %pose = Se3([1, 1, 1, 0, 0, 0]);
            %delta_pose = Se3([0, 0, 0, deg_per_sec, 0, 0]);
            %rotate_around_pose = Se3([pose.getTranslation() - delta_pose.getSo3() * pose.getTranslation(); deg_per_sec; 0; 0]);
            %new_pose = rotate_around_pose * rotate_around_pose * pose;
            %new_pose.matrix()
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
            end_time = 1;
            frames_per_second = 10;
            dt = 1/frames_per_second;
            trajectory = x_k0(1);
            for t=start_time:dt:(end_time-dt)
                x_k0 = ekf_se3.getState();
                x_k1 = ekf_se3.f(x_k0, 0, dt);
                ekf_se3.setState(x_k1);
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
        
        function J = computeJacobian(self)
            J_se3 = self.state_SE3.Dx_this_mul_exp_x_at_0();
            J_se3_dot = zeros(7,6);
            J = [J_se3; J_se3_dot];
        end
        
        function x_predicted = f(self, x, u, dt)
            %! Predict given current pose and velocity
            % Constant velocity model: update position based on last known velocity
            pose = self.state_SE3;
            % compute delta pose over time frame dt
            delta_pose = self.state_SE3_dot * dt;
            % decouple rotation and translation components by removing the
            % motion of the frame origin due to the rotational component
            % that effects the translation. The results is that rotational forces do not impact positional
            % updates.
            newTrans = pose.getTranslation() - delta_pose.getSo3() * pose.getTranslation() + delta_pose.getTranslation();
            delta_pose.setTranslation(newTrans);
            % left composition of transformations
            x_predicted(1) = delta_pose * self.state_SE3;
            % constant velocity assumption
            x_predicted(2) = self.state_SE3_dot;
        end
        
        % // new interface
        function out = predict(self, x, dt)
            fprintf(1,"SystemModel::predict\n");
            out(1) =  x(1) * x(2);
            %out(1) = (Se3.exp(dt * x(2).log()) * Se3.exp(x(1).log()));
            out(2) = x(2);
        end
    end
end