clear; close all; clc;
folder2 = '/Users/rongfan/Desktop/Growth_Model/Latex';

%%
FileName = fullfile(folder2,'gamma.png');
figure('Position',[0 0 400 200]); 
f1 = @(x) 1;
fplot(f1,[1,2],'linewidth',1);
hold on;
f2 = @(x) x+1;
fplot(f2,[1,2],'linewidth',1);
hold on;
f3 = @(x) 2*x;
fplot(f3,[1,2],'linewidth',1);
hold on; 
xline(1,'--')
hold on; 
xline(2,'--')
xlim([0.75 2.5])
ylim([0 5])
xticks([1,2]); xticklabels({'N-1','N'})
yticks([1,2]); yticklabels({'1','e^{B(N-1)}'})
legend({'\gamma_K','\gamma_L(i,0)','\gamma_H(i,0)'})
title('Task productivity')
saveas(gcf,FileName); 


%%
FileName = fullfile(folder2,'kappa_N.png');
figure('Position',[0 0 300 250]); 
f = @(x) -(x-2)^2+4;
fplot(f,[0,1],'linewidth',1)
xlabel('Scientist \epsilon'); ylabel('Innovation \kappa')
xticks([0,1]); yticks(3); yticklabels('\kappa(1)')
title('Innovation function')
saveas(gcf,FileName); 

%%
a = 1; b = 1.5; theta = 0.5;
FileName = fullfile(folder2,'growth1.png');
figure('Position',[0 0 300 250]); 
f = @(x) 0;
h1 = fplot(f,[0,b+1],'linewidth',1,'Color','b');
hold on; 
f = @(x) a;
h2 = fplot(f,[0,(1-theta)*a],'linewidth',1,'Color','r');
hold on;
f = @(x) (1/theta)*(a-x);
h3 = fplot(f,[(1-theta)*a,a],'linewidth',1,'Color','r');
hold on; 
f = @(x) 0;
h4 = fplot(f,[a,b+1],'linewidth',1,'Color','r');
xlabel('\rho'); ylabel('g')
xticks([0.2,a]); xticklabels({'\SigmaB\mu^N','b/\mu^h'});yticks(0);
ylim([0,2]);
legend([h1 h2],'bg^N','Bg^h')
title('Growth Rate with High \mu^N')
saveas(gcf,FileName); 


FileName = fullfile(folder2,'growth2.png');
figure('Position',[0 0 300 250]); 
f = @(x) b-a;
h1 = fplot(f,[0,theta*a+(1-theta)*b],'linewidth',1,'Color','b');
hold on; 
f = @(x) (1/theta)*(b-x);
h2 = fplot(f,[theta*a+(1-theta)*b,b],'linewidth',1,'Color','b');
hold on;
f = @(x) 0;
h3 = fplot(f,[b,b+1],'linewidth',1,'Color','b');
hold on; 
f = @(x) a
h4 = fplot(f,[0,(1-theta)*b],'linewidth',1,'Color','r');
hold on;
f = @(x) (1/theta)*(a-x)+(b-a)*(1-theta)/theta;
h5 = fplot(f,[(1-theta)*b,theta*a+(1-theta)*b],'linewidth',1,'Color','r');
hold on; 
f = @(x) 0;
h6 = fplot(f,[theta*a+(1-theta)*b,b+1],'linewidth',1,'Color','r');
xlabel('\rho'); ylabel('g')
xticks([a,b]); xticklabels({'b/\mu^h','\SigmaB\mu^N'});yticks(0);
ylim([0,2]);
legend([h1 h4],'bg^N','Bg^h')
title('Growth Rate with low \mu^N')
saveas(gcf,FileName); 

%%
FileName = fullfile(folder2,'kappa_NI.png');
figure('Position',[0 0 300 250]); 
f = @(x) -(x-2)^4+4;
fplot(f,[0,1],'linewidth',1,'Color','b')
hold on;
f = @(x) -(x+2)^2-3;
fplot(f,[0,1],'linewidth',1,'Color','r')
xlabel('Scientist \epsilon'); ylabel('Innovation \kappa')
xticks([0,1]); yticks(3); yticklabels('\kappa(f1)')
legend('\kappa^N','\kappa^I')
title('Innovation function')
saveas(gcf,FileName); 

%%
FileName = fullfile(folder2,'I1.png');
figure('Position',[0 0 300 250]); 
f = @(x) 1;
h1 = fplot(f,[0,4],'r','LineWidth',1);
hold on;
f = @(x) -((x-4)/3)^2+1;
h2 = fplot(f,[4,7],'r','LineWidth',1);
hold on;
h3 = fplot(f,[1,4],'b','LineWidth',1);
hold on;
f = @(x) 1;
h4 = fplot(f,[4,8],'b','LineWidth',1);
hold on; 
xline(4,'--')
str = {'S_I>0','S_N>0'};text(3.5,0.5,str);
str = {'No','innovation'};text(0.25,0.75,str);
str = {'No','automation'};text(6,0.75,str);
str = {'Abundant capital'};text(0.25,0.25,str);
str = {'Scarce capital'};text(5,0.25,str);
xlabel('R');
xticks([1,4,7])
xticklabels({'\rho_{min}','\rho^*','\rho_{max}'})
ylabel('$$\tilde{I}$$', 'Interpreter', 'LaTeX')
yticks(1)
saveas(gcf,FileName); 


%%
FileName = fullfile(folder2,'I2.png');
figure('Position',[0 0 300 250]); 
f = @(x) 1;
h1 = fplot(f,[0,4],'r','LineWidth',1);
hold on;
f = @(x) -((x-4)/3)^2+1;
h2 = fplot(f,[4,7],'r','LineWidth',1);
hold on;
h3 = fplot(f,[1,4],'b','LineWidth',1);
hold on;
f = @(x) 1;
h4 = fplot(f,[4,8],'b','LineWidth',1);
hold on; 
f = @(x) 1;
h5 = fplot(f,[0,4],'k','LineWidth',2);
hold on;
f = @(x) -((x-4)/1.5)^2+1;
h6 = fplot(f,[4,5.5],'k','LineWidth',2);
hold on;
f = @(x) 0;
h7 = fplot(f,[5.5,8],'k','LineWidth',2);
xline(4,'--')
str = {'S_I>0','S_N>0'};text(3.5,0.5,str);
str = {'No','innovation'};text(0.25,0.75,str);
str = {'No','automation'};text(6,0.75,str);
str = {'Abundant capital'};text(0.25,0.25,str);
str = {'Scarce capital'};text(5,0.25,str);
xlabel('R');
xticks([1,4,7])
xticklabels({'\rho_{min}','\rho^*','\rho_{max}'})
ylabel('$$\tilde{I}$$', 'Interpreter', 'LaTeX')
yticks(1)
legend(h5, 'I^*')
saveas(gcf,FileName); 

%%
FileName = fullfile(folder2,'I3.png');
figure('Position',[0 0 300 250]); 
x = [5 5.1 5.4];
f = @(x) 1;
h1 = fplot(f,[0,4],'r','LineWidth',1);
hold on;
f = @(x) -((x-4)/3)^2+1;
h2 = fplot(f,[4,7],'r','LineWidth',1);
hold on;
h3 = fplot(f,[1,4],'b','LineWidth',1);
hold on;
f = @(x) 1;
h4 = fplot(f,[4,8],'b','LineWidth',1);
hold on; 
xline(4,'--')
f = @(x) 1;
h5 = fplot(f,[0,4],'k','LineWidth',2);
hold on;
f = @(x) -((x-4)/1.5)^2+1;
y(1) = f(x(1));
h6 = fplot(f,[4,5.5],'k','LineWidth',2);
hold on;
f = @(x) 0;
h7 = fplot(f,[5.5,8],'k','LineWidth',2);
hold on;
f = @(x) 1;
h8 = fplot(f,[0,4],'m','LineWidth',2);
hold on;
f = @(x) -((x-4)/2)^2+1;
y(2) = f(x(2)); y(3) = f(x(3))
h9 = fplot(f,[4,6],'m','LineWidth',2);
hold on;
f = @(x) 0;
h10 = fplot(f,[6,8],'m','LineWidth',2);
hold on; 
plot([0 x(1)],[y(1) y(1)],'k','LineWidth',1) 
text(x(1)+0.2,y(1),'Initial automation')
xlabel('R');
xticks([1,4,7])
xticklabels({'\rho_{min}','\rho^*','\rho_{max}'})
ylabel('$$\tilde{I}$$', 'Interpreter', 'LaTeX')
yticks(1)
legend([h5,h8], {'I^*','I^*'''})
saveas(gcf,FileName); 

FileName = fullfile(folder2,'I4.png');
figure('Position',[0 0 300 250]); 
x = [5 5.1 5.4];
f = @(x) 1;
h1 = fplot(f,[0,4],'r','LineWidth',1);
hold on;
f = @(x) -((x-4)/3)^2+1;
h2 = fplot(f,[4,7],'r','LineWidth',1);
hold on;
h3 = fplot(f,[1,4],'b','LineWidth',1);
hold on;
f = @(x) 1;
h4 = fplot(f,[4,8],'b','LineWidth',1);
hold on; 
xline(4,'--')
f = @(x) 1;
h5 = fplot(f,[0,4],'k','LineWidth',2);
hold on;
f = @(x) -((x-4)/1.5)^2+1;
y(1) = f(x(1));
h6 = fplot(f,[4,5.5],'k','LineWidth',2);
hold on;
f = @(x) 0;
h7 = fplot(f,[5.5,8],'k','LineWidth',2);
hold on;
f = @(x) 1;
h8 = fplot(f,[0,4],'m','LineWidth',2);
hold on;
f = @(x) -((x-4)/2)^2+1;
y(2) = f(x(2)); y(3) = f(x(3))
h9 = fplot(f,[4,6],'m','LineWidth',2);
hold on;
f = @(x) 0;
h10 = fplot(f,[6,8],'m','LineWidth',2);
hold on; 
plot([0 x(1)],[y(1) y(1)],'k','LineWidth',1) 
hold on;
plot([0 x(2)],[y(2) y(2)],'m','LineWidth',1) 
text(x(2)+0.2,y(2),'Without human capital')
hold on;
plot([0 x(3)],[y(3) y(3)],'m','LineWidth',1) 
text(x(3)+0.2,y(3),'With human capital')
xlabel('R');
xticks([1,4,7])
xticklabels({'\rho_{min}','\rho^*','\rho_{max}'})
ylabel('$$\tilde{I}$$', 'Interpreter', 'LaTeX')
yticks(1)
legend([h5,h8], {'I^*','I^*'''})
saveas(gcf,FileName); 