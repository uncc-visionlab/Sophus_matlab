%
% Kalman Filter SE(3) using GTSAM
%
clc;
addpath("/home/arwillis/gtsam-dev/toolbox");
addpath("..");


%% Create the controls and measurement properties for our example
import gtsam.*

modelQ = noiseModel.Diagonal.Sigmas(0.1*ones(12,1));
modelR = noiseModel.Diagonal.Sigmas(0.1*ones(12,1));

%% Create an KalmanFilter object
KF = KalmanFilter(12);

%testPose = Pose3([[-1/sqrt(2) 1/sqrt(2) 0; sqrt(2), sqrt(2) 0; 0 0 1], [0 0 0]'; 0 0 0 1]);
%% Create the Kalman Filter initialization point
% initial state value
rad_per_sec = 45*pi/180;
x_initial = [0, 0, 0, 0, 1, 1, ...
    0, 0, rad_per_sec, 1, 0, 0]';
% we have confidence in our initial state
P_initial = 0.01*eye(12,12);

%% Create an KF object
state = KF.init(x_initial, P_initial);
%% Run iteration 1
%import gtsam.*
global deltaT;
deltaT = 1.0;
B = zeros(12,12);
u = ones(12,1);
%num_diff = NumericalDiff(@f, 12, 12, 'Central');
%num_diff = NumericalDiff(@f, 12, 12, 'Five-Point');
%num_diff = NumericalDiff(@f, 12, 12, 'Seven-Point');
num_diff = NumericalDiff(@f, 12, 12, 'Nine-Point');
x_k = x_initial;

Pose3.Expmap(x_k(1:6))
Pose3.Expmap(x_k(7:12)*deltaT)

for step=1:3
    %x_k = state.mean()
    T_deltak = Pose3.Expmap(deltaT*x_k(7:12))
    T_kplus1 = T_deltak.matrix()*(Pose3.Expmap(x_k(1:6)).matrix());
    x_k(1:6) = Pose3.Logmap(Pose3(T_kplus1));
    x_kNew = T_deltak.AdjointMap()*x_k(7:12);
    T_kplus1_2 = Pose3.Expmap(x_kNew).matrix()
    %x_kNew2 = Pose3.AdjointMap_(x_k(7:12))
    %transform = Pose3.Expmap(x_k(1:6)).matrix()
    transform = T_kplus1 ;
%     %Pose3.AdjointMap_(x_k(7:12))
%     F = [Pose3.ExpmapDerivative(x_k(7:12)),Pose3.ExpmapDerivative(x_k(1:6));
%         zeros(6,6),  eye(6,6)];
%     F = [Pose3.LogmapDerivative(Pose3.Expmap(-x_k(7:12)))*Pose3.LogmapDerivative(Pose3.Expmap(x_k(1:6))),Pose3.LogmapDerivative(Pose3.Expmap(x_k(1:6)));
%         zeros(6,6),  eye(6,6)];
%     df = num_diff.df(x_k);
%     dfRotPose = Rot3.ClosestTo(df(1:3,1:3));
%     dfRotVel = Rot3.ClosestTo(df(1:3,7:9));
%     %df(1:3,1:3) = dfRotPose.matrix();
%     %df(1:3,7:9) = dfRotVel.matrix();
%     df;
%     %state = KF.predict(state, df, B, u, modelQ);
%     state = f(x_k);
%     x_kplus1 = f(x_k);
%     x_k = x_kplus1;
%     %state = KF.update(state, H, z2, modelR);
end

function out = SO3hat(x)
  out =  [0.0, -x(3), +x(2); +(x3), 0.0, -x(1); -x(2), +x(1), 0.0];
end

function out = SE3hat(x)
  out =  [SO3hat(x(1:3)), [x(4);x(5);x(6)]; zeros(1,4)];
end

% allows multiplying a homogeneous point, p,  by a member of the Lie
% algebra, psi, as OpenCircle(p)*psi
% openCircle(p)*psi = SE3hat(psi)*p = T*p
function out = openCircle(tangent)
    out = [-SO3hat(tangent(4:6)), tangent(1:3)*eye(3,3); zeros(1,6)]
end


function newpose_SE3 = transformSE3(x_k)
import gtsam.*
global deltaT;
newpose_SE3 = Pose3.Expmap(x_k(1:6)).compose(Pose3.Expmap(x_k(7:12)*deltaT));
%newpose_SE3 = Pose3.Expmap(x_k(7:12)*deltaT).compose(Pose3.Expmap(x_k(1:6)));
end

function x_kplus1 = f(x_k)
import gtsam.*
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
newpose_SE3 = transformSE3(x_k);
x_kplus1 = zeros(12,1);
x_kplus1(1:6,1) = Pose3.Logmap(newpose_SE3);
% constant velocity assumption
%x_predicted(2) = self.state_SE3_dot;
x_kplus1(7:12,1) = x_k(7:12);
end

