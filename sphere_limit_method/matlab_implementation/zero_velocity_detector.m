%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%> @file zero_velocity_detector.m
%>
%> @brief Functions for implementing different zero-velocity detection 
%> algorithms. 
%>
%> @details Functions for implementing different zero-velocity detection 
%> algorithms, as well as a wrapper function for easy to use purpose. The
%> settings used by the wrapper function "zero_velocity_detector()" is  
%> specified in the file \a setting.m. Details about the detectors can be 
%> found in papers
%> 
%> \li <A href="http://dx.doi.org/10.1109/TBME.2010.2060723">Zero-Velocity Detection -- An Algorithm Evaluation</A> 
%> \li <A href="http://dx.doi.org/10.1109/IPIN.2010.5646936">Evaluation of Zero-Velocity Detectors for Foot-Mounted Inertial Navigation Systems</A>
%>   
%>
%> @authors Isaac Skog, John-Olof Nilsson
%> @copyright Copyright (c) 2013 OpenShoe, ISC License (open source)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  funtion [zupt T] = zero_velocity_detector(u) 
%
%>
%> @brief Wrapper function for running the zero-velocity detection 
%> algorithms. 
%>
%> @details A wrapper function that runs the zero-velocity detection 
%> algorithm that was specified in the file \a setting.m. 
%>
%> @param[out]  zupt       Vector with the detector decsions. [ true = zero velocity, false = moving]    
%> @param[out]  T          The test statistics of the detector 
%> @param[in]   u          The IMU data vector.     
%>
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [zupt T]=zero_velocity_detector(u)


% Global struct holding the simulation settings
global simdata;

% Allocate memmory
zupt=zeros(1,length(u));

% Run the desired detector type. Each detector return a vector with their 
% calculated test statistics T. 
switch simdata.detector_type
    
    case 'GLRT'
        T=GLRT(u);    
    
    case 'MV'
        T=MV(u);
        
    case 'MAG'
        T=MAG(u);
        
    case 'ARE'
        T=ARE(u);
        
    otherwise
        disp('The choosen detector type not recognized. The GLRT detector is used')
        T=GLRT(u);
end

% Check if the test statistics T are below the detector threshold. If so, 
% chose the hypothesis that the system has zero velocity 
W=simdata.Window_size;
for k=1:length(T)
    if T(k)<simdata.gamma
       zupt(k:k+W-1)=ones(1,W); 
    end    
end

% Fix the edges of the detector statistics
T=[max(T)*ones(1,floor(W/2)) T max(T)*ones(1,floor(W/2))];
end





%% SUBFUNCTIONS 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  funtion T = GLRT(u)
%
%> @brief Function that runs the generalized likelihood test (SHOE detector). 
%>
%> @param[out]  T          The test statistics of the detector 
%> @param[in]   u          The IMU data vector.     
%>
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function T=GLRT(u)

global simdata;

g=simdata.g;
sigma2_a=simdata.sigma_a^2;
sigma2_g=simdata.sigma_g^2;
W=simdata.Window_size;


N=length(u);
T=zeros(1,N-W+1);

for k=1:N-W+1
   
    ya_m=mean(u(1:3,k:k+W-1),2);
    
    for l=k:k+W-1
        tmp=u(1:3,l)-g*ya_m/norm(ya_m);
        T(k)=T(k)+u(4:6,l)'*u(4:6,l)/sigma2_g+tmp'*tmp/sigma2_a;    
    end    
end

T=T./W;

end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  funtion T = MV(u)
%
%> @brief Function that runs the acceleration moving variance detector. 
%>
%> @param[out]  T          The test statistics of the detector 
%> @param[in]   u          The IMU data vector.     
%>
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function T=MV(u)

global simdata;

sigma2_a=simdata.sigma_a^2;
W=simdata.Window_size;

N=length(u);
T=zeros(1,N-W+1);



for k=1:N-W+1
    
    ya_m=mean(u(1:3,k:k+W-1),2);
    
    for l=k:k+W-1
        tmp=u(1:3,l)-ya_m;
        T(k)=T(k)+tmp'*tmp;    
    end    
end

T=T./(sigma2_a*W);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  funtion T = MAG(u)
%
%> @brief Function that runs the acceleration magnitude detector. 
%>
%> @param[out]  T          The test statistics of the detector 
%> @param[in]   u          The IMU data vector.     
%>
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function T=MAG(u)

global simdata;
g=simdata.g;
sigma2_a=simdata.sigma_a^2;
W=simdata.Window_size;

N=length(u);
T=zeros(1,N-W+1);


for k=1:N-W+1
    for l=k:k+W-1
        T(k)=T(k)+(norm(u(1:3,l))-g)^2;    
    end    
end

T=T./(sigma2_a*W);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  funtion T = ARE(u)
%
%> @brief Function that runs the angular rate energy detector. 
%>
%> @param[out]  T          The test statistics of the detector 
%> @param[in]   u          The IMU data vector.     
%>
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function T=ARE(u)


global simdata;

sigma2_g=simdata.sigma_g^2;
W=simdata.Window_size;

N=length(u);
T=zeros(1,N-W+1);


for k=1:N-W+1
    for l=k:k+W-1
        T(k)=T(k)+norm(u(4:6,l))^2;    
    end    
end

T=T./(sigma2_g*W);
end