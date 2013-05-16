close all;
clear;
clc;

%====================================%
% setting up the scaling factor
%====================================%

%-constants
g=9.81

%- particle props
dp=0.0008
rhoP=200
Vp=dp^3*pi/6
np=19531

%- scaling
dParcel=dp*2
VParcel=dParcel^3*pi/6
nParcel=np*Vp/VParcel
dragForceScale=dParcel/dp
