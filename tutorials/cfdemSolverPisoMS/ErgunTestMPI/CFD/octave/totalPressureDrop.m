close all;
clear;
clc;

%====================================%
% simulation data 1
%====================================%
rhoG = 5			% density in kg/m3
%path = '../probes/0/p';
path = '../postProcessing/probes/0/p';
columns=22;
headerlines=4;
data = loaddata(path,columns,headerlines);
data=transpose(data);
[x,y]=size(data)
dp_sim = (data(:,2)-data(:,y))*rhoG;
t_sim = data(:,1);
%fprintf('final pressureDrop of sim = %f Pa\n',dp_sim(length(dp_sim)) )

%====================================%
% analytical calculation
%====================================%

%===================
% Ergun Equation
%===================
fprintf('\ncalc Ergun eqn:\n')

% calc density of clump
% Note: for drag calc the clump diameter is used-> reduced density of clump
dp = 0.002			% clump diameter
ds_=0.001*0.3530494*2 		% diameter of sphere
ns_=10;              		% nr of spheres per clump
rhos_=1000;          		% density of spheres
VclumpToVbody=ns_*(ds_/dp)^3 	% ratio of the tot volume of the spheres to the volume of the body with d_
rhoP = VclumpToVbody*rhos_ 	% density of clump

nP_=2500;			% nr of clumps
phip = 1			% sphericity
Ustart = 0.
Uend = 2
timeStepSize = 0.001;            % time interval of pressure data
Tstart = 0;
Tend = t_sim(length(t_sim));
deltaU=(Uend-Ustart)/((Tend-Tstart)/timeStepSize);
U = Ustart+deltaU:deltaU:Uend;  % velocity over time
L = 0.0212;			% length of bed filled with particles
D = 0.0276;			% diameter of bed

% calc epsilon from bedheight
%Vpartiles=dp^3*pi/6*nP_;	    % tot clump vol
Vpartiles=ds_^3*pi/6*nP_*ns_;	% tot particle vol
Vfilled=D^2*pi/4*L;		%volume filled with particles
epsilon = 1-Vpartiles/Vfilled   % void fraction

%epsilon = 0.656968

Ua = U / epsilon;		% physical velocity
nuG = 1.5*10^-5			% kinemat Visk in m2/s
muG = nuG*rhoG			% dynam visc in Pa s

dpErgun= L * (
                150*((1-epsilon)^2/epsilon^3)*((muG.*U)/(phip*dp)^2) 
              +1.75*((1-epsilon)/epsilon^3)*((rhoG.*U.^2)/(phip*dp))
        );

fprintf('NOTE: this pressure is divided by density (according to CFD solver)\n')
fprintf('so the result does not depend on pressure\n')

%fprintf('final pressure drop (Ergun eqn)= %f Pa\n',dpErgun)

%==================================
% min fluidization velocity in m/s
%==================================
%rhoP = 2000                      % particle density in kg/m3
g = 9.81                        % gravity m/s2

Umf = dp^2*(rhoP-rhoG)*g/(150*muG)*(epsilon^3*phip^2)/(1-epsilon);
ReMF = Umf*dp*rhoG/muG;
if(ReMF<20)
    fprintf('applying eqn1 for Umf.\n')
elseif(ReMF>20 && ReMF<1000)
    fprintf('applying eqn1 for Umf.\n')
elseif (ReMF>=1000)
    fprintf('applying eqn2 for Umf.\n')
    Umf = sqrt(dp*(rhoP-rhoG)*g/(1.75*rhoG)*epsilon^3*phip);
    ReMF = Umf*dp*rhoG/muG;
end
Umf
ReMF

dpUmf= L * (
                150*((1-epsilon)^2/epsilon^3)*((muG.*Umf)/(phip*dp)^2) 
              +1.75*((1-epsilon)/epsilon^3)*((rhoG.*Umf.^2)/(phip*dp))
        );

%====================================%
% plot data
%====================================%

fig=figure(1)
plot(U,dpErgun,U,dp_sim,[Umf,Uend],dpUmf*ones(1,2))
title("Ergun pressure drop vs. simulation")
a=strcat("analytical (Ergun), Umf=",num2str(Umf),", dpUmf=",num2str(dpUmf));
legend(a,"simulation")
xlabel("velocity in [m/s]")
ylabel("pressure drop [Pa]")
axis([0,Uend,0,dpErgun(length(dpErgun))])

%print('cfdemSolverPiso_settlingTest.eps','-deps2')
print -color "cfdemSolverPisoMS_ErgunTestMPI.eps"
SimName="ErgunTestMPI_sphereOfSpheres"
print(fig,strcat("figure_",SimName,".png"));

