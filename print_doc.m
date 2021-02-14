close all

D = importdata('Wing_Polar_Graph_0CD.csv');

% plot WEISSINGER vs XFLR5 --> Cd vs Cl
figure(10)
plot(D.data(:,1),D.data(:,2),'k','LineWidth',5);
hold on
plot(Cd_vec,Cl_vec,'or','LineWidth',5);
grid on
grid minor
xlabel('$C_{D}$','Interpreter','latex');
ylabel('$C_{L}$','Interpreter','latex');
legend('XFLR5','WEISSINGER','Interpreter','latex');