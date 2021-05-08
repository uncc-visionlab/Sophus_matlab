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