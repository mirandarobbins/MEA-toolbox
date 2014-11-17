%% movie example 
% This example calls  "AD_MED64_LFP_movie_Miranda.m", 
% ...which then calls "sweep2movie_LFPInterp3D.m" internally
data.raw(:,51,:)=0; % this hides the stim channel when it is = Ch 51

AD_MED64_LFP_movie_Miranda(data)