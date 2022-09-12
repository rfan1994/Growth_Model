cd '/Users/rongfan/Desktop/Growth_Model/Matlab'
clear; close all; clc;
folder1 = '/Users/rongfan/Desktop/Growth_Model/Matlab/Static';
Color(1,:) = [0 0.4470 0.7410];
Color(2,:) = [0.8500 0.3250 0.0980];
Color(3,:) = [0.9290 0.6940 0.1250];
Color(4,:) = [0.4940 0.1840 0.5560];
Color(5,:) = [0.4660 0.6740 0.1880];
Color(6,:) = [0.3010 0.7450 0.9330];
Color(7,:) = [0.6350 0.0780 0.1840];


%%
close all;
cd '/Users/rongfan/Desktop/Growth_Model/Matlab/Static';
dinfo = dir('BGP*');
names = [];
data = [];
for K = 1:length(dinfo)
    filename = dinfo(K).name;     %just the name
    table = readtable(filename);
    names = [names,table.Properties.VariableNames];
    data = [data,table2array(table)];
end
data = data(data(:,3)==1,:);

rho = 0.0152; theta = 0.8;

figure;
for n = [1,3]
    rr = data(data(:,1)==n,22);
    I_tilde = data(data(:,1)==n,16);
    g = data(data(:,1)==n,32);
    h(n) = plot(rr,I_tilde,'LineWidth',1);
    hold on; 
    [~,i] = min(abs(rho+theta*g-rr));
    line([rr(i),rr(i)],[0,I_tilde(i)],'Color','k','LineWidth',1,'LineStyle','--')
    hold on; 
    line([rr(i),0],[I_tilde(i),I_tilde(i)],'Color','k','LineWidth',1,'LineStyle','--')
    hold on; 
end
legend([h(1),h(3)], {'Low ES, High CA','High ES, Low CA'},'Location','Northeast')
xlabel('Interest rate'); ylabel('Automation')
xlim([0.03,0.09]); ylim([0,1])
xticks([])

%%
for n = [1,3]
    figure;
    rr = data(data(:,1)==n,22);
    I_tilde = data(data(:,1)==n,16);
    g = data(data(:,1)==n,32);
    h(1) = plot(rr,I_tilde,'LineWidth',1); 
    hold on; 
    [~,i] = min(abs(rho+theta*g-rr));
    line([rr(i),rr(i)],[0,I_tilde(i)],'Color','k','LineWidth',1,'LineStyle','--')
    hold on; 
    line([rr(i),0],[I_tilde(i),I_tilde(i)],'Color','k','LineWidth',1,'LineStyle','--')
    xlim([0.03,0.055]); ylim([0,1])
    hold on; 
    rr = data(data(:,1)==n+3,22);
    I_tilde = data(data(:,1)==n+3,16);
    g = data(data(:,1)==n+3,32);
    h(2) = plot(rr,I_tilde,'LineWidth',1);
    hold on; 
    [~,i] = min(abs(rho+theta*g-rr));
    line([rr(i),rr(i)],[0,I_tilde(i)],'Color','k','LineWidth',1,'LineStyle','--')
    hold on; 
    line([rr(i),0],[I_tilde(i),I_tilde(i)],'Color','k','LineWidth',1,'LineStyle','--')
    legend([h(1),h(2)],  {'Before wave','After wave'},'Location','Northeast')
    xlabel('Interest rate'); ylabel('Automation')
    xlim([0.03,0.09]); ylim([0,1])
    xticks([])
end 
