


clear all
% close all 

load('Results/XTRG.mat')



figure(1)
hold on, box on, grid on
[beta,I] = sort(beta);
ener = ener(I);
free_energy = free_energy(I);
plot((1./beta), ener, '^')
temp = 1./beta;

for i = 1:length(temp)-1
    cv(i) = (ener(i+1) - ener(i))/(temp(i+1)-temp(i));
    xx(i) = (temp(i+1)+temp(i))/2;
end

figure(2)
hold on, box on, grid on
plot((xx),log10(cv),'o')

% figure(3)
% hold on, box on, grid on
% plot(log(temp), -temp.*(free_energy)/5, '-o')

% %% EXACT MODEL XY
L = double(N);
J = 1;

for k = 1:L
    eps(k) = J*cos(k*pi/(L+1));
end

temp = linspace(50, 1e-3, 10000);
beta = 1./temp;

for j = 1:length(temp)
    zz = 1;
    for k = 1:L
        zz = zz*(1 + exp(-beta(j)*eps(k)));
    end
    Z(j) = zz;
end

f = -temp/L.*log(Z);
logz = -log(Z)/L;
beta_f = beta.*f;
for i = 1:length(temp)-1
    ener_beta(i) = (logz(i+1)-logz(i))/(beta(i+1)-beta(i));
    xxb(i) = (beta(i+1)+beta(i))/2;
end
tempb = 1./xxb;
for i = 1:length(tempb)-1
    cv(i) = (ener_beta(i+1) - ener_beta(i))/(tempb(i+1)-tempb(i));
    xx(i) = (tempb(i+1)+tempb(i))/2;  
end
% can compute the energy 



figure(1)
plot(1./xxb, ener_beta,'-')

figure(2)
plot((xx), log10(cv),'-')
legend('XTRG: $\chi = 50$','analytical','Interpreter','Latex')
title( strcat('1D XY model, N=',num2str(N)),'Interpreter','Latex')
xlabel('log(T)','Interpreter','Latex')
ylabel('$\log(c_v)$','Interpreter','Latex')
set(gca, 'fontsize',20)
set(gca,'TickLabelInterpreter','latex')

% figure(3)
% plot(log(temp), f)
