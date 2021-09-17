% Simulation results: RMS error vs. # slices
load incoherent % nSlices x permutations
rmsErrorInc = sqz(mean(rmsError,2)); % incoherent method, average RMS 
% error over 50 permutations

load conventional % nSlices x permutations
rmsErrorConv = rmsError; % conventional method, RMS error 

figure, plot(1:5, rmsErrorConv,'LineWidth',4)
xticks(1:5), set(gca, 'fontsize',24), axis([1 5 0 inf])
s = ['RMS Error (' char(176) 'C)'];
xlabel('Number of Slices'), ylabel(s)
hold on, plot(1:5,rmsErrorInc,'LineWidth',4)
legend('Conventional','Incoherent','Location','NorthWest')
legend('boxoff')

