clc, clf, clear, close all

load x_traj_total_control;
load x_traj_total_uncontrol;
t = 0:0.01:20;

control_quiver = x_total_control(:,1:100:end);
uncontrol_quiver = x_total_uncontrol(:,1:100:end);

figure(1)
plot(x_total_control(1,:),x_total_control(2,:))
hold on
plot(x_total_uncontrol(1,:),x_total_uncontrol(2,:))
plot(x_total_control(1,1),x_total_control(2,1),'*')
arrowh(x_total_control(1,:),x_total_control(2,:),'b',[],[2,5,8]);
arrowh(x_total_uncontrol(1,:),x_total_uncontrol(2,:),'r',[],5:5:95);

legend('controlled','uncontrolled','initial point','location','northwest','fontsize',13)
xlabel('$x_1$','interpreter','latex','fontsize',18)
ylabel('$x_2$','interpreter','latex','fontsize',18)
grid on

figure(2)
plot(t,x_total_control)
xlabel('$t$','interpreter','latex','fontsize',18)
ylabel('$x$','interpreter','latex','fontsize',18)
legend('$x_1$','$x_2$','$x_3$','$x_4$','interpreter','latex','fontsize',12)
grid on

figure(3)
plot(t,x_total_uncontrol)
xlabel('$t$','interpreter','latex','fontsize',18)
ylabel('$x$','interpreter','latex','fontsize',18)
legend('$x_1$','$x_2$','$x_3$','$x_4$','interpreter','latex','location','northwest','fontsize',12)
grid on