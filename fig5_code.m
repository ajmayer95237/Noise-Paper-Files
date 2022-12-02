ms = 0:1:30;
C = 1;
mpt = 20;
epsilon = .01;
figure(1)
for i = 1:4
    plot(ms,((ms-mpt).^2+epsilon).^(1/2)+(ms-mpt),'LineWidth',1)
    hold on
    epsilon = epsilon*10;
end
legend('\epsilon = 0.01','\epsilon = 0.1', '\epsilon = 1', '\epsilon = 10')
xlabel('Free mRNA (m)')
ylabel('Transition rate (min^{-1}) into droplet phase ')