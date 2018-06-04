%> @mainpage The OpenShoe Matlab Implemenation
%>
%> \section sec1 Introduction
%> This is the documentation for the OpenShoe Matlab Implemenation. The
%> OpenShoe Matlab Implementation is a Matlab script (library) for processing
%> inertial measurement unit (IMU) data using a Kalma filter based zero-velocity 
%> aided inertial navigation system algorithm. The main (skeleton) file that 
%> should be called to run the algorithm is \a main.m. All settings for the 
%> algorithm are done in the file \a settings.m. The processing of the data is
%> done by the functions zero_velocity_detector and ZUPTaidedINS, located 
%> in the files \a zero_velocity_detector.m and \a ZUPTaidedINS.m, respectively. 
%> The result from running the algorithm is plotted by calling the script 
%> \a view_data.m. 
%>
%> \section sec2 The zero-velocity aided inertial navigation algorithm
%> The zero-velocity aided inertial navigation system algorithm is
%> implemented using a complimentary feedback filtering structure. That is, 
%> the inertial navigation system works as the backbone of the system and 
%> a Kalman filter is used to estimate (track) the perturbations (errors) 
%> in the inertial navigation system. When a zero-velocity observation 
%> is done, the Kalman filter estimates the current perturbations (errors)
%> in the navigation state estimate of the inertial navigation system; the
%> estimated perturbations are feedback into the inertial navigation system 
%> to correct its internal states. The user can choose between four 
%> different state space models to be used in the Kalman filter. The 
%> default model is a nine state model, having the position, velocity, and
%> attitude perturbations as states.  The user can then choose between 
%> adding sensor biases errors and/or scale factors errors as additional 
%> states to the default model.        
%>
%> To determine when a zero-velocity update should be applied in the Kalman
%> filter, a zero-velocity detection algorithm is used. The zero-velocity 
%> detection algorithm monitors the signal measured by the inertial 
%> measurement unit, and based upon the prior information about the signal 
%> at different motion dynamics it chooses between the hypotheses that the 
%> navigation system is stationary or is moving. The user can choose 
%> between four different zero-velocity detection algorithms, the SHOE 
%> (\a GLRT.m) detector, the acceleration moving variance (\a MV.m) 
%> detector, the acceleration magnitude (\a MAG.m) detector, and the 
%> angular rate energy (\a ARE.m) detector. Details about these detectors 
%> can be found in the papers   
%> 
%> \li <A href="http://dx.doi.org/10.1109/TBME.2010.2060723">Zero-Velocity Detection -- An Algorithm Evaluation</A> 
%> \li <A href="http://dx.doi.org/10.1109/IPIN.2010.5646936">Evaluation of Zero-Velocity Detectors for Foot-Mounted Inertial Navigation Systems</A>
%>
%> \section sec3 Algorithm settings and configuartions
%> All settings for the algorithm and the inertial measurement data file
%> that should be used are controlled from the function settings.m. The 
%> default settings are such that the algorithm produces a good output 
%> when processing, the with the program, provided inertial measurement 
%> unit data and with the default state space model. For other data sets or when using 
%> the higher order state space models the settings may have to be tuned in
%> order for the algorithm to produce a good result. Information on how to 
%> tune the system can be found in the paper       
%>
%> \li <A href="http://dx.doi.org/10.1109/IPIN.2010.5646939">Performance characterisation of foot-mounted ZUPT-aided INSs and other related systems</A>
%>
%> \section sec4 The inertial measurement unit data
%> The inertial measurement unit data that comes with the OpenShoe Matlab
%> Implementation code has been recorded using a MicroStrain 3DX-GX2 inertial
%> measurement unit, with a dynamic range of +-18g and 1200 deg/s, and a 
%> sample rate of 250 Hz. The inertial measurement unit was mounted in the 
%> sole of the right side shoe of the user, and the user was walking in a 
%> closed loop trajectory where he returned to his starting position within
%> +-1 cm. The gate speed for the different data sets are 5 km/h and 7 km/h.   