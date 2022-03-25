clear all;
clc;

%% Script to plot the angular velocities and the particle positions

A = importdata('../../DEM/post/angular_velocity_no_coupling.txt',' ',1);
B = importdata('../../DEM/post/position_no_coupling.txt',' ',1);
C = importdata('../../DEM/post/angular_velocity.txt',' ',1);
D = importdata('../../DEM/post/position.txt',' ',1);

pos1 = B.data();
omega1 = A.data();
pos2 = D.data();
omega2 = C.data();

time = omega1(:,1);
omegax1 = omega1(:,2);
omegay1 = omega1(:,3);
omegaz1 = omega1(:,4);

time2 = omega2(:,1);
omegax2 = omega2(:,2);
omegay2 = omega2(:,3);
omegaz2 = omega2(:,4);

figure
plot(time,omegaz1,'-.-',time2,omegaz2,'Linewidth',1.5)
xlabel('Time (s)')
ylabel('Angular velocity (1/s)')
legend('One-way coupling','Two-way coupling')
axis([0 0.5 0 11])
set(gca,'FontSize',12)
print('angular_velocity_compare.eps')
