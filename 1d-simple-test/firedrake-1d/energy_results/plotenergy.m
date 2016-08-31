clc
close all


T = dlmread('quad_3_256.txt');
figure 
set(gca,'fontsize', 14)
semilogy( T(:,1), T(:,2))
title('Plot of the drift from initial energy')
xlabel('Wave periods')
ylabel('Error in energy')