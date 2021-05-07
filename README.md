# CUDARecv
An open-source parallelized direct position estimation-based GPS receiver

Stanford Navigation and Autonomous Vehicles Laboratory

## Literature References
Matthew Peretic and Grace X. Gao, Design of a Parallelized Direct Position Estimation-Based GNSS Receiver, Navigation: Journal of the Institute of Navigation, vol. 68, no. 1, pp. 21-39, Dec 2020. doi: 10.1002/navi.402. [paper](https://stanford.box.com/s/fe5obgok7kkbj9hvl1hxkbtkcoezur52)


## Development Configuration
### Host Configuration: Personal Computer with Ethernet Connection to Device
- Ubuntu 16
- CUDA 9.0.252 (obtained from Jetpack 3.3.1 -- [developer.nvidia.com](https://developer.nvidia.com/embedded/jetpack-3_3))
- Eigen 3.3.5 !!! with modification to work with CUDA !!!
    - Run [get_Eigen_Jetson.sh](cudarecv/scripts/get_Eigen_Jetson.sh) to download Eigen and apply the patch to work with CUDA.
    - This script will make one edit to Eigen/Core -- view the script to see what is changed.
- UHD >= 3.10.2 (obtained from [Ettus Research](https://github.com/EttusResearch/uhd))
    - Requires Boost -- can be installed with "sudo apt-get install libboost-all-dev"
        - Requires Mako
- Nsight Eclipse IDE

### Device Configuration: NVIDIA Jetson TX2
- Ubuntu 16
- CUDA 9.0.252 (obtained from Jetpack 3.3.1 -- [developer.nvidia.com](https://developer.nvidia.com/embedded/jetpack-3_3))
- Eigen 3.3.5 !!! with modification to work with CUDA !!!
    - Run [get_Eigen_Jetson.sh](cudarecv/scripts/get_Eigen_Jetson.sh) to download Eigen and apply the patch to work with CUDA.
    - This script will make one edit to Eigen/Core -- view the script to see what is changed.
- UHD >= 3.10.2 (obtained from [Ettus Research](https://github.com/EttusResearch/uhd))
    - Requires Boost -- can be installed with "sudo apt-get install libboost-all-dev"
        - Requires Mako

### Cross-Compilation
The NVIDIA Jetson TX2 uses the AARCH64 architecture -- different than a conventional personal computer
- If compiling on the host, Nsight must cross-compile for AARCH64
- Consider configuring Nsight to compile on the device instead



## Running CUDARecv
### Setting up the project in Nsight
- Run [nsight_config.sh](cudarecv/scripts/nsight_config.sh)
- In Nsight,
    - New > Makefile Project with Existing Code
        - Existing Code Location: the outermost cudarecv folder
        - Toolchain: Cross GCC
    - Once the project is created, right-click the project in the Project Explorer > Properties
        - Build > Builder Settings > Build location > ${workspace_loc:/cudarecv}/build
        - Build > Builder Settings > Makefile generation > "Generate Makefiles Automatically" SHOULD NOT BE CHECKED
        - Build > Target Systems > Select remote connection... > Choose and configure your ethernet-connected Jetson TX2



### Configuring a DPE Run
Three files are required to run the DPE receiver: 
- Dataset
    - Configured for Ettus Research USRP sampling at 2.5MHz (digitally downconverted from the GPS L1 band to 0Hz) with I and Q samples of 16 bits each. 20ms of samples are read-in and processed together.
    - The sampling frequency and window size are specified in [dpeflow](cudarecv/dsp/src/dpeflow.cpp) and [dpinit](cudarecv/modules/src/dpinit.cpp).
    - The data size of samples is hardcoded into [sampleblock](cudarecv/modules/src/sampleblock.cu).
- Initialization file
    - Recommend obtaining by running scalar tracking in PyGNSS for >10s.
    - If you can't get a lock on the data in PyGNSS scalar tracking (eg, particularly bad signal environment), you can conceivably make your own following the template handoff csv provided.
        - Provide the non-clock corrected receiver time estimate: rxTime 
        - Provide the initial state estimate: X_ECEF
        - Provide the location in the dataset where you want to begin reading the dataset: bytes_read
        - Provide the list of satellites you want to track: prn_list
        - Provide the channel parameters by differencing the true receiver state and the satellites' states (see discussion at end of README): rc, ri, fc, fi, 
        - Provide the TOW value, the code periods elapsed since the start of the run when the TOW arrived, and the code periods elapsed when beginning reading the dataset: TOW, cp_timestamp, cp 
    - File settings can be found in [dpeflow](cudarecv/dpe/src/dpeflow.cpp).
- RINEX file
    - Designed for Version 2.10 NAVIGATION DATA
    - Recommend downloading nistDDD0.YYn filetype from [CDDIS](ftp://cddis.nasa.gov/gnss/data/daily/)
        - For year 20YY and day of the year DDD, go to folder 20YY > DDD > YYn > ctrl+F "nist"
        - Unzip the download
    - File settings can be found in [dpeflow](cudarecv/dpe/src/dpeflow.cpp).
- (Optional) DPE grid
    - A position grid can be read-in (instead of generated within the code) through a csv file of the ENU-delta-t offsets (see rngrid3.csv for example).
    - File settings can be found in [dpeflow](cudarecv/dsp/src/dpeflow.cpp).


### Processing a Dataset
- Run the project. A prompt 'gnss$' will appear in the console.
- Type 'newflow dpe' and hit enter
- Type 'loadflow 0' and hit enter
- Type 'startflow 0' and hit enter




### Use Our Sample Datasets
Move the [demofiles](demofiles) folder to your Ubuntu Desktop and run CUDARecv following the "Processing a Dataset" section above.
- [static_opensky](demofiles/static_opensky_20180704_190000_usrp6_2500kHz.dat) is a 45-second GPS dataset simulated for July 4, 2018 at 19:00 UTC Time. The satellite PRNs simulated are 2, 3, 6, 12, 17, 19, 24, and 28. The simulated receiver is stationary and located in Urbana-Champaign, IL, USA.
- [handoff_params](demofiles/handoff_params_usrp6.csv) provides the channel initialization data from scalar tracking in PyGNSS.
- [nist1860](demofiles/nist1860.18n) is the CDDIS navigation data for the simulated dataset, providing the GPS ephemeris.
- [rngrid3](demofiles/rngrid3.csv) is a custom position-time-domain manifold grid for the BCM module. 
If there are errors with opening these files, provide their absolute paths in [dpeflow.cpp](cudarecv/dsp/src/dpeflow.cpp).



### Generating Your Own Handoff Parameters
- Specify the details of your dataset in PyGNSS's [setting.py](pygnss/setting.py).
- Run PyGNSS's [1_Data_reduct_scalar.py](pygnss/1_Data_reduct_scalar.py).
- Use the "handoff_params.csv" file generated.



## CUDARecv Development and You
- The DPE flow is intended to run on *minutes* of data at most -- things like TOW and ephemeris do not update, and no tracking is performed to update the channel parameters.
- The [DataLogger](cudarecv/modules/src/datalogger.cu) and [SampleBlock](cudarecv/modules/src/sampleblock.cu) modules have timeouts of 1.5s which will crash the flow when exceeded. If your hardware does not complete one iteration in this time, extend the timeouts in those two files.
- Error catching has not been fully implemented. 
- Configuration for the DPE flow is done in [main.cu](cudarecv/cudarecv/src/main.cu) and [dpeflow.cpp](cudarecv/dsp/src/dpeflow.cpp).
- The "Flow-Module" paradigm is detailed in [Matthew Peretic's MS Thesis](http://hdl.handle.net/2142/105723).
- When creating a new Flow:
    - Use the functions (or create new ones) defined in [flow.cu](cudarecv/dsp/src/flow.cu).
    - Add the file name to the [CMakeLists in dsp](cudarecv/dsp/src/CMakeLists.txt).
    - All Inputs, Outputs, and Parameters to constituent modules should be connected in your new flow file.
    - Add your flow to [flowmgr.h](cudarecv/dsp/inc/flowmgr.h):
        - Add a #include to your new flow's header file.
        - Add your flow using the RegFlow function.
    - "push back" your flow in "_regis" in [flowmgr.cpp](cudarecv/dsp/src/flowmgr.cpp).
    - Create a class declaration in [dsp.h](cudarecv/dsp/inc/dsp.h).
    - Include your flow's header in [main.cu](cudarecv/cudarecv/src/main.cu).
- When creating a new Module:
    - Use the functions (or create new ones) defined in [module.cpp](cudarecv/modules/src/module.cpp).
    - Add the file name to the [CMakeLists in modules](cudarecv/modules/src/CMakeLists.txt).
- The syntax for console commands can be found in [cmdFlow.h](cudarecv/console/inc/cmdFlow.h).
- Module parameters can be changed through the console after 'loadflow' and before 'startflow' by the command 'setparam *flow* *module* *param* *value*'.
    - *flow*: flow number, ex 0
    - *module*: module name, ex SampleBlock
    - *param*: parameter name, ex SamplingFrequency
    - *value*: parameter value, ex 2500000
- The [guhd](cudarecv/guhd) folder contains some starter code that may be useful for USRP streaming.


## Commentary on Channel Parameters
- rc, ri, fc, fi can be back-calculated by differencing the receiver's true state and the satellites' true states. 
- These parameters are the receiver's track of the transmissions sent by each satellite
- The better these estimates, the better the replica created by the DPE receiver, and the more accurate the correlation scores
- Consider for your work whether you want to study solutions to the manifold (poor initial state, good initial channel parameters) or study how to generate the replica (poor initial channel parameters). 
    - The authors hold that the latter question is more challenging, but is the bigger roadblock to broader use of DPE. 
    - Acquiring good estimates of channel parameters is a tracking problem, and may be addressed with tracking loops like conventional GNSS receivers use. However, this mitigates some of the benefits of using DPE.
    - For a primer and the authors' views on this problem, a survey is presented in [this paper by Peretic and Gao](https://web.stanford.edu/~gracegao/publications//journal//2020%20Navigation_DPE%20receiver%20I_Matt%20Peretic.pdf)



# PyGNSS

An open-source sequential GPS receiver with scalar tracking (conventional two-steps) and DPE (one-step) localization algorithms.

Stanford Navigation and Autonomous Vehicles Laboratory

## Literature References
E. Wycoff, Y. Ng, and G. X. Gao, "Python GNSS Receiver: An Object-Oriented Software Platform Suitable for Multiple Receivers," *GPS World Magazine*, February 2015. [Online](https://web.stanford.edu/~gracegao/publications//magazine//052-057_GPS_Feb15_v2.pdf)


## Development Configuration

### Personal Computer
- Python 2.7
- Numpy 1.15.1

## Running PyGNSS
- Run one of the numerically-indexed .py files in [pygnss](pygnss)
    - [0_Data_reduction.py](pygnss/0_Data_reduction.py) performs scalar tracking initialization then hands off and runs DPE.
    - [1_Data_reduct_scalar.py](pygnss/1_Data_reduct_scalar.py) performs scalar tracking only and can print DPE handoff parameters to a file.
    - [2_Generate_ephemerides.py](pygnss/2_Generate_ephemerides.py) acquires and stores ephemerides from the dataset being processed.
    - [3_Data_reduct_dp.py](pygnss/3_Data_reduct_dp.py) performs DPE only, loading in external channel parameters for initialization.
- Results are saved in folders 'post-simulator' and 'pre-simulator', which are created on-run.

## Configuring PyGNSS
- Set parameters such as sampling frequency, etc. accordingly in [setting.py](pygnss/setting.py)
- Change the flow of the run itself in the .py files.


# Development Team
Grace Gao Research Group
- Lead PI: Dr. Grace X. Gao
- Ramya Bhamidipati
- Shubhendra Chauhan
- Arthur Hsi-Ping Chu
- Enyu Luo
- Yuting Ng
- Matthew Peretic
- ...and contributions through many helpful discussions with other members of the research group!


# Copyright 
Copyright (c) 2020 Stanford Navigation and Autonomous Systems Laboratory

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
