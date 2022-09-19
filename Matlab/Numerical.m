cd '/Users/rongfan/Desktop/Growth_Model/Matlab'
clear; close all; clc;
folder1 = '/Users/rongfan/Desktop/Growth_Model/Matlab/20220918_2133/Two_Group';
folder2 = '/Users/rongfan/Desktop/Growth_Model/Latex';
Color(1,:) = [0 0.4470 0.7410];
Color(2,:) = [0.8500 0.3250 0.0980];
Color(3,:) = [0.9290 0.6940 0.1250];
Color(4,:) = [0.4940 0.1840 0.5560];
Color(5,:) = [0.4660 0.6740 0.1880];
Color(6,:) = [0.3010 0.7450 0.9330];
Color(7,:) = [0.6350 0.0780 0.1840];


%%
close all;
cd(folder1)
dinfo = dir('IR*1.txt');
names = [];
data1 = [];
for K = 1:length(dinfo)
    filename = dinfo(K).name;     %just the name
    table = readtable(filename);
    names = [names,table.Properties.VariableNames];
    data1 = [data1,table2array(table)];
end

dinfo = dir('IR*3.txt');
data2 = [];
for K = 1:length(dinfo)
    filename = dinfo(K).name;     %just the name
    table = readtable(filename);
    data2 = [data2,table2array(table)];
end

eta = 0.1136;

variable = {'Automation','Labor share','Interest rate','Wage inequality',...
            'Automation rate','Innovaton rate','Human capital growth rate','Human capital gap'};

y1 = data1(:,[3,37,33,36,19,18,20,5]);       
y1(:,2) = data1(:,37)+data1(:,38)+eta;

y2 = data2(:,[3,37,33,36,19,18,20,5]);   
y2(:,2) = data2(:,37)+data2(:,38)+eta;


T = data1(1:31,1);

figure('position',[0,0,1000,300])
for n = 1:8
    subplot(2,4,n)
    plot(T,y1(1:31,n),'linewidth',1)
    hold on;
    plot(T,y2(1:31,n),'linewidth',1)
    title(variable(n))
    if (n==1) 
        legend('Without human capital','With human capital','Location','southeast')
    end
end
filename = append(folder2,'/Transition.png');       
exportgraphics(gcf,filename)




