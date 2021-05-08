%
%
% Installing GTSAM (Ubuntu 18.04)
%
% 1. Compile GTSAM
%
% GTSAM development from Github
% Compile the development github source set (-DGTSAM_BUILD_PYTHON=ON) for
% Python bindings to be compiles (did not work 5/5/2021)
% git clone https://github.com/borglab/gtsam.git
% cd gtsam
% You will need to have a github account with a valid ssh key or add one as
% decribed here: https://docs.github.com/en/github/authenticating-to-github/adding-a-new-ssh-key-to-your-github-account
% If this is the case the command below can succeed
% ./update_wrap.sh
% Otherwise you can "rm -rf wrap" and git clone
% git clone https://github.com/borglab/wrap.git
% you will need the latest python3 version of pyparsing
% pip3 install pyparsing
% mkdir build
% cd build
% cmake -DGTSAM_INSTALL_MATLAB_TOOLBOX=ON -DGTSAM_BUILD_PYTHON=OFF -DMatlab_ROOT_DIR=/usr/local/bin/matlab/R2020b -DCMAKE_INSTALL_PREFIX:PATH=$HOME/gtsam-dev -DGTSAM_TOOLBOX_INSTALL_PATH:PATH=$HOME/gtsam-dev/toolbox -DGTSAM_PYTHON_VERSION=3 ..
% make
% make check
% make install
%
% GTSAM version rc4.1
% cmake -DGTSAM_INSTALL_MATLAB_TOOLBOX=ON -DMATLAB_ROOT=/usr/local/bin/matlab/R2020b -DMEX_COMMAND=/usr/local/bin/matlab/R2020b/bin/mex -DCMAKE_INSTALL_PREFIX:PATH=$HOME/gtsam-4.1 -DGTSAM_TOOLBOX_INSTALL_PATH:PATH=$HOME/gtsam-4.1/toolbox -DGTSAM_BUILD_PYTHON=ON -DGTSAM_PYTHON_VERSION=3 ..
%
% 2. Move or delete MATLAB's standard C++ libraries so MATLAB uses the
% systems updated C++ versions of these libraries.
%
% Make old standard C++ libs unavailable to MATLAB
% mv libstdc++.so.6 libstdc++.so.6.bak
% mv libstdc++.so.6.0.25 libstdc++.so.6.0.25.bak
%
% 3. Allow Ubuntu to find your custom location for the GTSAM shared libraries
%
% Before running MATLAB when using GTSAM
% export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/home/arwillis/gtsam-dev/lib
%
% 4. Tell MATLAB where the GTSAM code and C++ interface library is located
% In MATLAB, to use GTSAM start your code with
% addpath("/home/arwillis/gtsam-dev/toolbox")
% import gtsam.*
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
    0, 0, rad_per_sec, 0, 0, 0]';
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
num_diff = NumericalDiff(@f, 12, 12, 'Five-Point');

for step=1:3
    x_k = state.mean()
    transform = Pose3.Expmap(x_k(1:6)).matrix()
    x_kplus1 = f(x_k);
    %Pose3.AdjointMap_(x_k(7:12))
    F = [Pose3.ExpmapDerivative(x_k(7:12)), eye(6,6);
        zeros(6,6),  eye(6,6)]
    df = num_diff.df(state.mean());
    dfRotPose = Rot3.ClosestTo(df(1:3,1:3));
    dfRotVel = Rot3.ClosestTo(df(1:3,7:9));
    df(1:3,1:3) = dfRotPose.matrix();
    df(1:3,7:9) = dfRotVel.matrix();
    df
    state = KF.predict(state, df, B, u, modelQ);
    %state = KF.update(state, H, z2, modelR);
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

