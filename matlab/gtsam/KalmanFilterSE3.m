

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
addpath("/home/arwillis/gtsam-dev/toolbox")


%% Create the controls and measurement properties for our example
import gtsam.*

%pose = Pose3.Expmap(initialState(1:6))
%J_pose = Pose3.ExpmapDerivative(initialState(1:6))
%velocity = Pose3.Expmap(initialState(7:12))
%J_velocity = Pose3.ExpmapDerivative(initialState(7:12))

Q = 0.01*eye(2,2);
H = eye(12,12);
z2 = [2.0, 0.0]';
z3 = [3.0, 0.0]';
modelR = noiseModel.Diagonal.Sigmas(0.1*ones(12,1));
R = 0.01*eye(2,2);

%% Create the set of expected output TestValues
expected0 = [0.0, 0.0]';
P00 = 0.01*eye(2,2);

expected1 = [1.0, 0.0]';
P01 = P00 + Q;
I11 = inv(P01) + inv(R);

expected2 = [2.0, 0.0]';
P12 = inv(I11) + Q;
I22 = inv(P12) + inv(R);

expected3 = [3.0, 0.0]';
P23 = inv(I22) + Q;
I33 = inv(P23) + inv(R);

%% Create an KalmanFilter object
KF = KalmanFilter(12);

%% Create the Kalman Filter initialization point
% initial state value
deg_per_sec = 45*pi/180;
x_initial = [0, 0, 0, 0, 0, 0, ...
    1, 0, 0, deg_per_sec, 0, 0]';

% we have confidence in our initial state
P_initial = 0.01*eye(12,12); 

%% Create an KF object
state = KF.init(x_initial, P_initial);
%EQUALITY('expected0,state.mean', expected0,state.mean);
%EQUALITY('expected0,state.mean', P00,state.covariance);

%% Run iteration 1
%import gtsam.*
deltaT = 1.0;
poseExp = state.mean();
F = [Pose3.ExpmapDerivative(poseExp(1:6)), eye(6,6) * deltaT;
    zeros(6,6),  eye(6,6)]
B = zeros(12,12);
u = ones(12,1);
modelQ = noiseModel.Diagonal.Sigmas(0.1*ones(12,1));

state = KF.predict(state, F, B, u, modelQ);
%EQUALITY('expected1,state.mean', expected1,state.mean);
%EQUALITY('P01,state.covariance', P01,state.covariance);
%state = KF.update(state, H, z1, modelR);
%EQUALITY('expected1,state.mean', expected1,state.mean);
%EQUALITY('I11,state.information', I11,state.information);

%% Run iteration 2
%import gtsam.*
poseExp = state.mean();
F = [Pose3.ExpmapDerivative(poseExp(1:6)), eye(6,6) * deltaT;
    zeros(6,6),  eye(6,6)]
state = KF.predict(state,F, B, u, modelQ);
%EQUALITY('expected2,state.mean', expected2,state.mean);
%state = KF.update(state,H,z2,modelR);
%EQUALITY('expected2,state.mean', expected2,state.mean);

%% Run iteration 3
%import gtsam.*
poseExp = state.mean();
F = [Pose3.ExpmapDerivative(poseExp(1:6)), eye(6,6) * deltaT;
    zeros(6,6),  eye(6,6)]
state = KF.predict(state,F, B, u, modelQ);
%EQUALITY('expected3,state.mean', expected3,state.mean);
%state = KF.update(state,H,z3,modelR);
%EQUALITY('expected3,state.mean', expected3,state.mean);
