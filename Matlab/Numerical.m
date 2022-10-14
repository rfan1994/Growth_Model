%% Balanced Growth 
cd '/Users/rongfan/Desktop/Growth_Model/Matlab'
clear; close all; clc;

folder1 = '/Users/rongfan/Desktop/Growth_Model/Matlab/20220921_1159/Two_Group';
folder2 = '/Users/rongfan/Desktop/Growth_Model/Latex';
Color(1,:) = [0 0.4470 0.7410];
Color(2,:) = [0.8500 0.3250 0.0980];
Color(3,:) = [0.9290 0.6940 0.1250];
Color(4,:) = [0.4940 0.1840 0.5560];
Color(5,:) = [0.4660 0.6740 0.1880];
Color(6,:) = [0.3010 0.7450 0.9330];
Color(7,:) = [0.6350 0.0780 0.1840];

cd(folder1)
dinfo = dir('BGP*.txt');
names = [];
data = [];
for K = 1:3
    filename = dinfo(K).name;     %just the name
    table = readtable(filename);
    names = [names,table.Properties.VariableNames];
    data = [data,table2array(table)];
end

tax_names = [];
tax = [];
filename = 'Policy.txt';
table = readtable(filename);
tax_names = [names,table.Properties.VariableNames];
tax = [tax,table2array(table)];

tau_I = tax(1:length(tax)/2+0.5,4);
tau_hH = tax([1,length(tax)/2+1.5:length(tax)],5);

variable = {'Automation','Labor share','Wage premium','Welfare inequality'...
            'Total welfare','Innovation rate','Human capital growth rate','Human capital gap'};
y_I = data(1:length(tax)/2+0.5,[16,28,25,12,10,34,36,18]);  
y_hH = data([1,length(tax)/2+1.5:length(tax)],[16,28,25,12,10,34,36,18]);
epsilon_H = 0.3; epsilon_L = 0.7; theta = 0.9;
y_I(:,5) = ((epsilon_H*data(1:length(tax)/2+0.5,10)+epsilon_L*data(1:length(tax)/2+0.5,11))./...
           (epsilon_H*data(1,10)+epsilon_L*data(1,11))).^(1/(1-theta))-1d0;
y_hH(:,5) = ((epsilon_H*data([1,length(tax)/2+1.5:length(tax)],10)+epsilon_L*data([1,length(tax)/2+1.5:length(tax)],11))./...
            (epsilon_H*data(1,10)+epsilon_L*data(1,11))).^(1/(1-theta))-1d0;
t = 1:11;
figure('position',[0,0,1000,300])
for n = 1:8
    subplot(2,4,n)
    plot(tau_I(t),y_I(t,n),'linewidth',1)
    hold on;
    plot(tau_hH(t),y_hH(t,n),'linewidth',1)
    if (n>4) 
        xlabel('Tax rate')
    end
    title(variable(n))
    if (n==1) 
        legend('Automation tax','Traning tax','Location','northwest')
    end
end

filename = append(folder2,'/Tax.png');       
exportgraphics(gcf,filename)


%% Transition
cd '/Users/rongfan/Desktop/Growth_Model/Matlab'
clear; close all; clc;
folder1 = '/Users/rongfan/Desktop/Growth_Model/Matlab/20221014_1327/Two_Group';
folder2 = '/Users/rongfan/Desktop/Growth_Model/Latex';
Color(1,:) = [0 0.4470 0.7410];
Color(2,:) = [0.8500 0.3250 0.0980];
Color(3,:) = [0.9290 0.6940 0.1250];
Color(4,:) = [0.4940 0.1840 0.5560];
Color(5,:) = [0.4660 0.6740 0.1880];
Color(6,:) = [0.3010 0.7450 0.9330];
Color(7,:) = [0.6350 0.0780 0.1840];

cd(folder1)
dinfo = dir('IR*1_1.txt');
names = [];
data1_1 = [];
for K = 1:length(dinfo)
    filename = dinfo(K).name;     %just the name
    table = readtable(filename);
    names = [names,table.Properties.VariableNames];
    data1_1 = [data1_1,table2array(table)];
end

dinfo = dir('IR*1_3.txt');
data1_3 = [];
for K = 1:length(dinfo)
    filename = dinfo(K).name;     %just the name
    table = readtable(filename);
    data1_3 = [data1_3,table2array(table)];
end

dinfo = dir('IR*2_1.txt');
data2_1 = [];
for K = 1:length(dinfo)
    filename = dinfo(K).name;     %just the name
    table = readtable(filename);
    data2_1 = [data2_1,table2array(table)];
end

dinfo = dir('IR*2_3.txt');
data2_3 = [];
for K = 1:length(dinfo)
    filename = dinfo(K).name;     %just the name
    table = readtable(filename);
    data2_3 = [data2_3,table2array(table)];
end

dinfo = dir('IR*3_1.txt');
data3_1 = [];
for K = 1:length(dinfo)
    filename = dinfo(K).name;     %just the name
    table = readtable(filename);
    data3_1 = [data3_1,table2array(table)];
end

dinfo = dir('IR*3_3.txt');
data3_3 = [];
for K = 1:length(dinfo)
    filename = dinfo(K).name;     %just the name
    table = readtable(filename);
    data3_3 = [data3_3,table2array(table)];
end

eta = 0.1136;
variable = {'Automation','Labor share','Interest rate','Wage premium',...
            'Automation rate','Innovaton rate','Human capital growth rate','Human capital gap'};

y1_1 = data1_1(:,[3,37,33,36,19,18,20,5]);       
y1_1(:,2) = data1_1(:,37)+data1_1(:,38)+eta;

y1_3 = data1_3(:,[3,37,33,36,19,18,20,5]);   
y1_3(:,2) = data1_3(:,37)+data1_3(:,38)+eta;

y2_1 = data2_1(:,[3,37,33,36,19,18,20,5]);       
y2_1(:,2) = data2_1(:,37)+data2_1(:,38)+eta;

y2_3 = data2_3(:,[3,37,33,36,19,18,20,5]);   
y2_3(:,2) = data2_3(:,37)+data2_3(:,38)+eta;
 
y3_1 = data3_1(:,[3,37,33,36,19,18,20,5]);       
y3_1(:,2) = data3_1(:,37)+data3_1(:,38)+eta;
 
y3_3 = data3_3(:,[3,37,33,36,19,18,20,5]);   
y3_3(:,2) = data3_3(:,37)+data3_3(:,38)+eta;

t = 1:81;
T = data1_1(t,1);

figure('position',[0,0,1000,300])
for n = 1:8
    subplot(2,4,n)
    plot(T,sgolayfilt(y1_1(t,n),3,11),'linewidth',1)
    hold on;
    plot(T,sgolayfilt(y1_3(t,n),3,11),'linewidth',1)
    title(variable(n))
    xlim([T(1),T(end)])
    if (n==5||n==6)
        ylim([4,8])
    elseif (n==7)
        ylim([0.25,0.5])
    end 
    if (n==1) 
        legend('Without human capital','With human capital','Location','southeast')
    end
end
filename = append(folder2,'/Transition.png');       
exportgraphics(gcf,filename)


figure('position',[0,0,1000,300])
for n = 1:8
    subplot(2,4,n)
    plot(T,y1_1(t,n),'linewidth',1)
    hold on;
    plot(T,y2_1(t,n),'linewidth',1)
    title(variable(n))
    xlim([T(1),T(end)])
    if (n==5||n==6)
        ylim([4,8])
    elseif (n==7)
        ylim([0.25,0.5])
    end 
    if (n==1) 
        legend('Benchmark','Automation tax','Location','southeast')
    end
end
filename = append(folder2,'/Transition2.png');       
exportgraphics(gcf,filename)


figure('position',[0,0,1000,300])
for n = 1:8
    subplot(2,4,n)
    plot(T,sgolayfilt(y1_3(t,n),3,11),'linewidth',1)
    hold on;
    plot(T,sgolayfilt(y2_3(t,n),3,11),'linewidth',1)
    hold on;
    plot(T,sgolayfilt(y3_3(t,n),3,11),'linewidth',1)
    title(variable(n))
    xlim([T(1),T(end)])
    if (n==5||n==6)
        ylim([4,8])
    elseif (n==7)
        ylim([0.25,0.5])
    end 
    if (n==1) 
        legend('Benchmark','Automation tax','Training tax','Location','southeast')
    end
end
filename = append(folder2,'/Transition3.png');       
exportgraphics(gcf,filename)




