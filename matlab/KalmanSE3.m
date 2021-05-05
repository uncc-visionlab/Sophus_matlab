
% Resources:
% Kalman for SE(3) - https://hal.archives-ouvertes.fr/hal-01122406/document
% On-manifold optimization - https://arxiv.org/pdf/1512.02363.pdf
% SE(3) vs. SO(3) x R3 - https://arxiv.org/pdf/1805.02543.pdf
% SE(3) On-manifold with a retration (to make it SO(3) x R3) - http://rpg.ifi.uzh.ch/docs/TRO16_forster.pdf
% SE(3) Lie Group theory and Jacobians of the SE(3) exponential and log
% maps - http://asrl.utias.utoronto.ca/~tdb/bib/barfoot_ser17.pdf
% SE(3) Tutorial and (some) useful derivatives - https://ingmec.ual.es/~jlblanco/papers/jlblanco2010geometry3D_techrep.pdf
% Kalman filtering - https://docs.ufpr.br/~danielsantos/ProbabilisticRobotics.pdf
% SE(3) and other parameterizations of 3D pose change - http://campar.in.tum.de/pub/tbirdal20183dv/tbirdal20183dv.pdf
classdef KalmanSE3 < handle
    % Extended Kalman filter for 3D pose using Lie Group SE(3)
    
    properties
        state_SE3
        state_SE3_dot
        num_diff
        dt
    end
    methods (Static)
        function trajectory = kalman_lie()
            clc;
            initialState = [0, 0, 0, 0, 0, 0, ...
                0.05, 0.1, 0, 0, 0.1, 0]';
            ekf_se3 = KalmanSE3(initialState);
            x = ekf_se3.getState();
            %x_mat = x(1).matrix()';
            %xdot_mat = x(2).matrix()';
            % Generate a trajectory
            num_steps = 7;
            %traj(1:num_steps) = Se3();
            v_init = initialState(7:12);
            traj(1) = Se3(initialState(1:6));
%            velSe3 = Se3(v_init);
            traj2 = zeros(6,num_steps);
            vel2 = zeros(6,num_steps);
            traj2(:,1) = initialState(1:6);
            vel2(:,1) = initialState(7:12);
            dt = 1;
            ekf_se3.setTimeStep(dt);

%             for i=2:size(traj,2)
%                 %traj(i) = Se3.exp(traj(i-1).log()) * Se3.exp(v_init);
%                 traj(i) = (traj(i-1) * dt) * velSe3;
%                 newpose = Se3.exp(traj2(:,i-1)) * Se3.exp(vel2(:,1));
%                 traj2(:,i) = newpose.log();
%                 x_mat_str = sprintf('%0.6g ',traj(i).log()');
%                 fprintf(1,"traj[%d]: %s\n",i,x_mat_str);
%             end

            fd = fopen('output.txt','w+');
            fprintf(fd,"#x_pred;x_mes;x_ekf\n");
            x_mat_str = sprintf('%0.6g ',x(1:6)');
            xdot_mat_str = sprintf('%0.6g ',x(7:12)');
            fprintf(1,"Initial State: %s, %s\n", x_mat_str, xdot_mat_str);
            fprintf(fd,"%s%s%s",x_mat_str,x_mat_str,x_mat_str);
            x_mes = Se3([0,0,0,0,0,0]');
            for i=2:num_steps
                x_ref = Se3.exp(traj2(:,i-1)) * Se3.exp(initialState(7:12) * dt);
                
                % Use reference velocity for testing
                % velocity_measurement.addPosition(x_mes);
                % velocity_measurement.addPosition(x_ref);
                
                % Simulate system (constant velocity model)
                % x = sys.f(x, u);
                % std::cout << "New State: " << x.x.transpose() << ", " << x.v.transpose() << std::endl;
                     

                % Predict state for current time-step using the filters
                x_k = ekf_se3.getState();
                x_kplus1 = ekf_se3.predict(x_k);
                ekf_se3.getJacobian(x_k);
                %x_kplus1_p =  (x_k(2) * dt) * x_k(1);
                %x_kplus1_pp = (Se3.exp(dt * x(2).log()) * Se3.exp(x(1).log()));
                %x_kplus1_ppp = (Se3.exp(dt * vel2(:,1)) * Se3.exp(traj2(:,i-1)));
                %x_kplus1_p_vec = x_kplus1_p.log();
                %x_kplus1_pp_vec = x_kplus1_pp.log();
                %x_kplus1_ppp_vec = x_kplus1_ppp.log();
                ekf_se3.setState(x_kplus1);
                %x_kplus1 = ekf_se3.getState();
                x_k_pos_str = sprintf('%0.6g ',x_k(1:6)');
                x_kplus1_pos_str = sprintf('%0.6g ',x_kplus1(1:6)');
                %x_kplus1_pos_str = sprintf('%0.6g ',x_kplus1_ppp(1:6)');
                x_kplus1_vel_str = sprintf('%0.6g ',x_kplus1(1:6)');
                
                fprintf(1,"Pred State: %s, %s\n", x_kplus1_pos_str, x_kplus1_vel_str);
                
                % Sensor update every 3 iteration, predict the rest of the time
                % if (i % 2 == 0)
                % {
                %     std::cout << "UPDATING VELOCITY" << std::endl;
                %     v_mes = velocity_measurement.v;
                %     ekf.update(velocity_measurement, v_mes);
                % }
                
                if (mod(i,3) == 0)
                    %fprintf(1, "UPDATING POSITION\n");
                    %x_mes = noise.addNoise(x_ref);
                    %ekf.update(pos_measurement, x_mes);
                end
                %auto cov = pos_measurement.getCovariance();
                %std::cout << "covariance:\n"
                %          << ekf.getCovariance().matrix() << std::endl;
                x_kplus1_corrected = ekf_se3.getState();

                x_pred_str = sprintf('%0.6g ',x_ref(1).log()');
                x_meas_str = sprintf('%0.6g ',x_mes(1).log()');
                x_kplus1_corrected_str = sprintf('%0.6g ',x_kplus1_corrected(1:6)');
                
                %fprintf(1,"x_pred: %s\n", x_kplus1_pos_str);
                fprintf(1,"x_pred: %s\n", x_pred_str);
                fprintf(1,"x_mes: %s\n", x_meas_str);
                fprintf(1,"x_ekf: %s\n", x_kplus1_corrected_str);
                %csv.write({x_ref, x_mes, x_ekf.x});
                fprintf(fd,"%s%s%s\n",x_k_pos_str, x_meas_str, x_kplus1_corrected_str);
            end
            fclose(fd);
            
            %State x_future;
            ekf_se3.setTimeStep(5);
            x_future = ekf_se3.predict(ekf_se3.getState());
            x_future_pos_str = sprintf('%0.6g ',x_future(1:6)');
            fprintf(1,"Future (dt=5): %s\n", x_future_pos_str);
        end
        
        function trajectory = demo()
            clc;
            rad_per_sec = 45*pi/180;
            initialState = [3, 3, 3, 0, 0, 0, ...
                3, 0, 0, rad_per_sec, 0, 0]';
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
            box.Pose = Se3(x_k0(1:6)).matrix();
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
            frames_per_second = 3;
            dt = 1/frames_per_second;
            ekf_se3.setTimeStep(dt);
            trajectory = x_k0(1:6);
            for t=start_time:dt:(end_time-dt)
                x_k0 = ekf_se3.getState();
                %x_k1 = ekf_se3.f(x_k0);
                x_k1 = ekf_se3.predict(x_k0);

                %J2 = ekf_se3.num_diff.df(x_k0)
                ekf_se3.getJacobian(x_k0);
                ekf_se3.setState(x_k1);
                transform = Se3(x_k1(1:6)).matrix()
                %det(transform(1:3,1:3))
                box.Pose = transform;
                show(box);
                title('Box');
                pause(dt);
                trajectory = [trajectory, x_k1(1:6)];
            end
            show(box);
            title('Box');
            pause(0.1);
        end
    end
    methods
        function self = KalmanSE3(varargin)
            if nargin == 0
                self.state_SE3 = zeros(6,1);
                self.state_SE3_dot = zeros(6,1);
            elseif nargin == 1 && size(varargin{1},1) == 12
                self.state_SE3 = varargin{1}(1:6);
                self.state_SE3_dot = varargin{1}(7:12);
            elseif nargin == 2 && size(varargin{1},1) == 6 && size(varargin{2},1) == 6
                self.state_SE3 = varargin{1};
                self.state_SE3_dot = varargin{2};
            else
                fprintf(1,"Could not initialize the KalmanSE3 filter.\n");
                return;
            end
            %self.num_diff = NumericalDiff(@self.f, 12, 12, 'Forward');
            self.num_diff = NumericalDiff(@self.f, 12, 12, 'Central');
            self.dt = 1.0;
            %self.state_SE3.matrix()
            %self.state_SE3_dot.matrix()
            %self.f(self.getState(), 0, 1);
            %trackingEKF()
        end
        
        function x = getState(self)
            x = [self.state_SE3; self.state_SE3_dot];
        end
        
        function setState(self, x)
            self.state_SE3 = x(1:6);
            self.state_SE3_dot = x(7:12);
        end
        
        function setTimeStep(self, dt)
            self.dt = dt;
        end
        
        function J = getJacobian(self, x)
            fprintf(1,"SystemModel::getJacobian\n");
            pose_SE3 = Se3(x(1:6));
            deltapose_SE3 = Se3(x(7:12) * self.dt);
            %newpose_SE3 = deltapose_SE3 * pose_SE3;
            newpose_SE3 = self.transform(x);
            %J_se3 = Se3.Dx_exp_x(newpose_SE3.log());
            %J_se32 = pose_SE3.Dx_this_mul_exp_x_at_0();
            %J_so3 = So3.right_JacobianExpmap(newpose_SE3.getSo3().log());
            %J_se3_dot = Se3.Dx_exp_x(deltapose_SE3.log());
            %J_se3_dot2 = deltapose_SE3.Dx_this_mul_exp_x_at_0();
            %J_r_pose = Se3.right_JacobianExpmap(newpose_SE3.log())
            J_r_pose = Se3.right_JacobianExpmap(x(1:6));
            J_r_deltapose = Se3.right_JacobianExpmap(x(7:12));
            J1 = [J_r_pose, J_r_deltapose; zeros(6,6), eye(6,6)]
            %J_rInv = Se3.right_JacobianLogmap(newpose_SE3.log())
            J2 = self.num_diff.df(x)
            J_r_pose(1:6,1:6)-J2(1:6,1:6)
            %J_se3 = [J_se3, eye(6); zeros(6,6), eye(6)]
            %J_so3xr3 = [[eye(3), J_so3; zeros(3,3), J_so3], eye(6); zeros(6,6), eye(6)];
            %fprintf(1,"SystemModel numerical jacobian:\n %s\n", mat2str(J2));
            %fprintf(1,"SystemModel jacobian:\n %s\n", mat2str(J));
            J = J1;
        end
        
        function x_kplus1 = f(self, x_k)
            %if (nargin < 3)
            %    dt = 1.0;
            %end
            %! Predict given current pose and velocity
            % Constant velocity model: update position based on last known velocity            
            %pose_SE3 = Se3(x_k(1:6));
            % compute delta pose over time frame dt
            %deltapose_SE3 = Se3(x_k(7:12) * dt);
            
            % decouple rotation and translation components by removing the
            % motion of the frame origin due to the rotational component
            % that effects the translation. The results is that rotational forces do not impact positional
            % updates.
            % We can apply the SE(3) retraction from eq (21) of https://arxiv.org/pdf/1512.02363.pdf
            %newTrans = pose_SE3.getTranslation() - deltapose_SE3.getSo3() * pose_SE3.getTranslation() + deltapose_SE3.getTranslation();
            %deltapose_SE3.setTranslation(newTrans);
            % left composition of transformations
            %x_predicted(1) = delta_pose * self.state_SE3;
            %newpose_SE3 = deltapose_SE3 * self.state_SE3;
            %newpose_SE3 = deltapose_SE3 * pose_SE3;
            %newpose_SE3 = pose_SE3 * deltapose_SE3;
            newpose_SE3 = self.transform(x_k);
            x_kplus1 = zeros(12,1);
            x_kplus1(1:6,1) = newpose_SE3.log();
            % constant velocity assumption
            %x_predicted(2) = self.state_SE3_dot;
            x_kplus1(7:12,1) = x_k(7:12);
        end
        
        % // new interface
        function x_kplus1 = predict(self, x_k)
            fprintf(1,"SystemModel::predict\n");
            x_kplus1 = zeros(12,1);
            %pose_SE3 = Se3(x(1:6));
            %deltapose_SE3 = Se3(x(7:12) * self.dt);
            % compute delta pose over time frame dt
            %newpose_SE3 = deltapose_SE3 * pose_SE3;
            newpose_SE3 = self.transform(x_k);
            x_kplus1(1:6,1) = newpose_SE3.log();
            x_kplus1(7:12,1) = x_k(7:12);
        end
        
        function newpose_SE3 = transform(self, x)
            pose_SE3 = Se3(x(1:6));
            deltapose_SE3 = Se3(x(7:12) * self.dt);
            %newpose_SE3 = deltapose_SE3 * pose_SE3;
            newpose_SE3 = pose_SE3 * deltapose_SE3;
            %newpose_SE3 = Se3(x(1:6) + x(7:12) * self.dt);
        end
    end
end