clear all
clc
close all

graphics_toolkit gnuplot

time = dlmread('./CFD/postProcessing/flowRatePatch/0/surfaceRegion.dat','\t',4,0)(:,1);
phi = dlmread('./CFD/postProcessing/flowRatePatch/0/surfaceRegion.dat','\t',4,0)(:,2);
alphaAir = dlmread('./CFD/postProcessing/patchAverage/0/surfaceRegion.dat','\t',4,0)(:,2);
alphaWater = dlmread('./CFD/postProcessing/patchAverage/0/surfaceRegion.dat','\t',4,0)(:,3);
alphaOil = dlmread('./CFD/postProcessing/patchAverage/0/surfaceRegion.dat','\t',4,0)(:,4);

volflowAir = phi.*alphaAir;
volflowWater = phi.*alphaWater;
volflowOil = phi.*alphaOil;

figure
hold on
plot(time,volflowAir*1000,'linewidth',5)
plot(time,volflowWater*1000,'r','linewidth',5)
plot(time,volflowOil*1000,'k','linewidth',5)
axis([0 20 0 0.45])
xlabel('Time, s','FontSize',14,'FontWeight','Bold')
ylabel('Volumetric flow rate, dm^3/s','FontSize',14,'FontWeight','Bold')
h=legend('Air  ', 'Water  ', 'Oil  ');
set(h,'FontWeight','Bold')
set(gca,'FontSize',14)
set(gca,'FontWeight','bold')

print -dpng volFlow.png
pause







