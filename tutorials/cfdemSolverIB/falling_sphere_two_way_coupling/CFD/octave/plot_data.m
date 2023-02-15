clear all;
clc;

% plot the angular velocity and the particle position

A = importdata('../../DEM/post/angular_velocity.txt',' ',1);
B = importdata('../../DEM/post/position.txt',' ',1);

omega = A.data();
pos = B.data();

time1 = omega(:,1);
omegax = omega(:,2);
omegay = omega(:,3);
omegaz = omega(:,4);

time2 = pos(:,1);
posx = pos(:,2);
posy = pos(:,3);
posz = pos(:,4);

figure(1)
plot(time1,omegaz,'-.-','Linewidth',1.5)
xlabel('Time (s)')
ylabel('Angular velocity (1/s)')
legend('Two-way torque coupling')
axis([0 0.5 0 11])
set(gca,'FontSize',12)
print('angular_velocity.eps')

figure(2)
plot(time2,posy,'-.-','Linewidth',1.5)
xlabel('Time (s)')
ylabel('Y-position (m)')
legend('Two-way torque coupling')
axis([0 0.5 0 0.1])
set(gca,'FontSize',12)
print('position.eps')

