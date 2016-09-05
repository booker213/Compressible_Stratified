clc
close all


T = dlmread('quad_16_1000.txt');
figure 
set(gca,'fontsize', 18)
semilogy( T(:,1), T(:,2))
title('Plot of the drift from initial energy')
xlabel('Wave periods')
ylabel('Error in energy')