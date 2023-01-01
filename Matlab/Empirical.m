clear; close all; clc;
folder1 =  '/Users/rongfan/Desktop/Growth_Model/STATA/EU_KLEMS';
folder2 = '/Users/rongfan/Desktop/Growth_Model/Latex';

figure('position',[0,0,800,250])
for n = 1:2
    switch n
        case 1
            filename = append(folder1,'/dH_shares1_LP.txt');
            table = readcell(filename);
        case 2
            filename = append(folder1,'/VAConLC_LP.txt');
            table = readcell(filename);
    end
    coef = table(3:4,3:8);
    coef = string(coef);

    T = 1:6;
    for i = 1:6
        coef(1,i) = regexp(coef(1,i),'\d+(\.)?(\d+)?','match');
        coef(2,i) = regexp(coef(2,i),'\d+(\.)?(\d+)?','match');
    end
    coef = str2double(coef(:,:));

    subplot(1,2,n)
    patch([T fliplr(T)],[coef(1,:)-1.036*coef(2,:) fliplr(coef(1,:)+1.036*coef(2,:))],'k','FaceAlpha',0.1,'EdgeColor','w');
    patch([T fliplr(T)],[coef(1,:)-1.96*coef(2,:) fliplr(coef(1,:)+1.96*coef(2,:))],'k','FaceAlpha',0.05,'EdgeColor','w');
    hold on;
    plot(T,coef(1,:),'linewidth',1);
    hold on;
    yline(0,'linewidth',1)
    xticks(1:6) 
    xlabel('Year')
    switch n
        case 1
            title('Share of Skilled Workers')
        case 2
            title('Labor Composition Growth')
    end
end
filename = append(folder2,'/LP.png');
exportgraphics(gcf,filename)
close all;


%%
clear; close all; clc;
folder1 =  '/Users/rongfan/Desktop/Growth_Model/STATA/ONET';
folder2 = '/Users/rongfan/Desktop/Growth_Model/Latex';

figure('position',[0,0,800,300])
for n = 1:2
    switch n
        case 1
            filename = append(folder1,'/LV_trend.csv');            
            table = readcell(filename);
        case 2
            filename = append(folder1,'/IM_trend.csv');
            table = readcell(filename);
    end
    table = readcell(filename);
    T = cell2mat(table(1,2:end));
    T = str2double(regexp(T,'\d+','match'));
    y = cell2mat(table(2:8,2:end));
    names = table(2:8,1);

    subplot(1,2,n)
    for i = 1:7
        dy = y(i,:)-y(i,1);
        plot(T,dy,'linewidth',1);
        hold on;
    end
    xlim([519 751])
    xticks([540 600 660 720])
    xticklabels({'2005','2010','2015','2020'})
    xlabel('Year')
    ylim([-1 1.8])
    if (n==1)
        ylabel('Cumulative change')
    end 
    
    switch n
        case 1
            title('Level')             
        case 2
            title('Importance')
            legend(names,'location','northeast') 

    end  
end
filename = append(folder2,'/trend.png');       
exportgraphics(gcf,filename)


%%
clear; close all; clc;
folder1 =  '/Users/rongfan/Desktop/Growth_Model/STATA/ONET';
folder2 = '/Users/rongfan/Desktop/Growth_Model/Latex';
for n = 1:2
    switch n
        case 1
            filename = append(folder1,'/LV_trend1.csv');
            table = readcell(filename);
        case 2
            filename = append(folder1,'/IM_trend1.csv');
            table = readcell(filename);
    end
    table = readcell(filename);
    T = cell2mat(table(1,3:end));
    T = str2double(regexp(T,'\d+','match'));
    y = cell2mat(table(2:end,2:end));
    names = table(2:8,1);

    figure('position',[0,0,1000,500])
    for i = 1:7
        subplot(2,4,i)
        for j = 1:3
            dy = y(i+7*(j-1),2:end)-y(i+7*(j-1),2);
            plot(T,dy,'linewidth',1);
            hold on;
        end  
        title(names(i))
        xlim([519 751])
        xticks([540 600 660 720])
        xticklabels({'2005','2010','2015','2020'})
        xlabel('Year')   
        
        switch n
            case 1
            	ylim([-0.5 1.5])
            case 2
            	ylim([-0.5 1])
        end        
        if (i==1 || i==5) 
            ylabel('Cumulative change')
        end
        if (i==1)
        	legend('Low exposure','Median exposure','High exposure')
        end 
    end

    
    switch n
        case 1
            filename = append(folder2,'/LV_trend1.png');
        case 2
            filename = append(folder2,'/IM_trend1.png');
    end  
    
    exportgraphics(gcf,filename); 
end



%%
clear; close all; clc;
folder1 =  '/Users/rongfan/Desktop/Growth_Model/STATA/ONET';
folder2 = '/Users/rongfan/Desktop/Growth_Model/Latex';
for n = 1:2
    switch n
        case 1
            filename = append(folder1,'/LV_trend2.csv');
            table = readcell(filename);
        case 2
            filename = append(folder1,'/IM_trend2.csv');
            table = readcell(filename);
    end
    table = readcell(filename);
    T = cell2mat(table(1,3:end));
    T = str2double(regexp(T,'\d+','match'));
    y = cell2mat(table(2:22,2:end));
    names = table(2:8,1);

    figure('position',[0,0,1000,500])
    for i = 1:7
        subplot(2,4,i)
        for j = 1:3
            dy = y(i+7*(j-1),2:end)-y(i+7*(j-1),2);
            plot(T,dy,'linewidth',1);
            hold on;
        end  
        title(names(i))
        xlim([519 751])
        xticks([540 600 660 720])
        xticklabels({'2005','2010','2015','2020'})
        xlabel('Year')   
        
        switch n
            case 1
            	ylim([-2 3])
            case 2
            	ylim([-1 2])
        end        
        if (i==1 || i==5) 
            ylabel('Cumulative change')
        end
        if (i==1)
        	legend('Low educaton','Median education','High education')
        end 
    end

    
    switch n
        case 1
            filename = append(folder2,'/LV_trend2.png');
        case 2
            filename = append(folder2,'/IM_trend2.png');
    end  
    
    exportgraphics(gcf,filename); 
end



%%
clear; close all; clc;
folder1 =  '/Users/rongfan/Desktop/Growth_Model/STATA/ATES';
folder2 = '/Users/rongfan/Desktop/Growth_Model/Latex';

for i = 1:3
    switch i
        case 1
            filename = append(folder1,'/l_aoe');
        case 2
            filename = append(folder1,'/l_aoe0');
        case 3
            filename = append(folder1,'/l_aoe1');
    end

table = readcell(filename);
aoe = cell2mat(table(3:end-1,end-2));
r = cell2mat(table(3:end-1,end-1));
l = cell2mat(table(3:end-1,end));

figure('position',[0,0,800,250])
subplot(1,2,1);
scatter(aoe,r,'filled')
hold on;
coefficients = polyfit(aoe,r,1);
xFit = linspace(min(aoe),max(aoe), 1000);
yFit = polyval(coefficients , xFit);
plot(xFit,yFit,'linewidth',1.5,'color', [0.8500 0.3250 0.0980]);
xlim([0 100])
xlabel('Automation Exposure')
ylabel('Percentage')
title('Training-Working Ratio')

subplot(1,2,2);
scatter(aoe,l,'filled')
hold on;
coefficients = polyfit(aoe,l,1);
xFit = linspace(min(aoe),max(aoe), 1000);
yFit = polyval(coefficients , xFit);
plot(xFit,yFit,'linewidth',1.5,'color', [0.8500 0.3250 0.0980]);
xlim([0 100])
xlabel('Automation Exposure')
ylabel('Percentage')
title('Training Parcitipation Rate')
    switch i
        case 1
            filename = append(folder2,'/train.png'); 
        case 2
            filename = append(folder2,'/train0.png'); 
        case 3
            filename = append(folder2,'/train1.png'); 
    end    
exportgraphics(gcf,filename)
end

filename = append(folder1,'/l_aoe0');
table = readcell(filename);
aoe0 = cell2mat(table(3:end-1,end-2));
l0 = cell2mat(table(3:end-1,end));

filename = append(folder1,'/l_aoe1');
table = readcell(filename);
aoe1 = cell2mat(table(3:end-1,end-2));
l1 = cell2mat(table(3:end-1,end));

figure('Position',[0 0 400 250]); 
scatter(aoe0,l0,'filled');
hold on;
scatter(aoe0,l1,'filled');
hold on;
coefficients = polyfit(aoe0,l0,1);
xFit = linspace(min(aoe0),max(aoe0), 1000);
yFit = polyval(coefficients , xFit);
p1 = plot(xFit,yFit,'linewidth',1.5,'color', [0 0.4470 0.7410]);
hold on
coefficients = polyfit(aoe1,l1,1);
xFit = linspace(min(aoe1),max(aoe1), 1000);
yFit = polyval(coefficients, xFit);
p2 = plot(xFit,yFit,'linewidth',1.5,'color', [0.8500 0.3250 0.0980]);
legend([p2 p1],'Skilled workers','Unskilled workers')
xlim([10,80])
ylim([0,100])
xlabel('Automation Exposure')
ylabel('Percentage')
title('Training Parcitipation Rate')

filename = append(folder2,'/train01.png'); 
saveas(gcf,filename)


%%
clear; close all; clc;
folder1 =  '/Users/rongfan/Desktop/Growth_Model/STATA/ONET';
folder2 = '/Users/rongfan/Desktop/Growth_Model/Latex';

filename = append(folder1,'/LV_trend1.csv');
table = readcell(filename);
T = cell2mat(table(1,2:end));
T = str2double(regexp(T,'\d+','match'));

filename = append(folder1,'/LV.txt');
table = readcell(filename);
coef = table(3:96,2:8);
coef = string(coef);
names = table(2,2:8);

for i = 1:7
    for j = 1:94
        coef(j,i) = regexp(coef(j,i),'[-]?\d+(\.)?(\d+)?','match');
        coef(j,i) = regexp(coef(j,i),'[-]?\d+(\.)?(\d+)?','match');
    end
end
coef = str2double(coef(:,:));
beta = coef(1:2:end,:);
st = coef(2:2:end,:);

figure('position',[0,0,800,300])
for i = 1:7
    subplot(2,4,i)
    patch([T fliplr(T)],[beta(:,i)-1.036*st(:,i); flipud(beta(:,i)+1.036*st(:,i))]','k','FaceAlpha',0.1,'EdgeColor','w');
    hold on;
    plot(T,beta(:,i),'linewidth',1);
    title(names(i))
    xlim([519 751])
    xticks([540 600 660 720])
    xticklabels({'2005','2010','2015','2020'})
    xlabel('Year')  
end
filename = append(folder2,'/LV.png');    
exportgraphics(gcf,filename)

filename = append(folder1,'/LV_group.txt');
table = readcell(filename);
coef = table(3:284,2:8);
coef = string(coef);
names = table(2,2:8);

for i = 1:7
    for j = 1:282
        coef(j,i) = regexp(coef(j,i),'[-]?\d+(\.)?(\d+)?','match');
        coef(j,i) = regexp(coef(j,i),'[-]?\d+(\.)?(\d+)?','match');
    end
end
l = length(coef);
coef = str2double(coef(:,:));
beta1 = coef(1:2:l/3,:);
st1 = coef(2:2:l/3,:);
beta2 = coef(l/3+1:2:2*l/3,:);
st2 = coef(l/3+2:2:2*l/3,:);
beta3 = coef(2*l/3+1:2:end,:);
st3 = coef(2*l/3+2:2:end,:);

figure('position',[0,0,800,300])
for i = 1:7
    subplot(2,4,i)
    patch([T fliplr(T)],[beta1(:,i)-1.036*st1(:,i); flipud(beta1(:,i)+1.036*st1(:,i))]','k','FaceAlpha',0.1,'EdgeColor','w');
    hold on;
    plot(T,beta1(:,i),'linewidth',1);
    title(names(i))
    xlim([519 751])
    xticks([540 600 660 720])
    xticklabels({'2005','2010','2015','2020'})
    xlabel('Year')  
    sgtitle('Skill Responses to Automation: \beta_2')
end
filename = append(folder2,'/LV_group1.png');    
exportgraphics(gcf,filename)

figure('position',[0,0,800,300])
for i = 1:7
    subplot(2,4,i)
    patch([T fliplr(T)],[beta3(:,i)-1.036*st3(:,i); flipud(beta3(:,i)+1.036*st3(:,i))]','k','FaceAlpha',0.1,'EdgeColor','w');
    hold on;
    plot(T,beta3(:,i),'linewidth',1);
    title(names(i))
    xlim([519 751])
    xticks([540 600 660 720])
    xticklabels({'2005','2010','2015','2020'})
    xlabel('Year')  
    sgtitle('Differences in Skill Responses to Automation: \beta_4')
end
filename = append(folder2,'/LV_group3.png');    
exportgraphics(gcf,filename)

%%
clear; close all; clc;
folder1 =  '/Users/rongfan/Desktop/Growth_Model/STATA/NLSY97';
folder2 = '/Users/rongfan/Desktop/Growth_Model/Latex';

filename = append(folder1,'/Wage_trend.csv');
table = readcell(filename);
T0 = cell2mat(table(2:20,1));
T1 = cell2mat(table(40:54,1));
y0 = cell2mat(table(2:39,4));
y1 = cell2mat(table(40:end,4));

figure('position',[0,0,800,250])

for n = 1:2
    subplot(1,2,n)
    switch n
        case 1
            plot(T0,log(y0(1:length(T0)))-log(y0(1)),'linewidth',1);
            hold on;
            plot(T0,log(y0(length(T0)+1:end))-log(y0(length(T0)+1)),'linewidth',1);
            title('Unskilled Worker')
        case 2    
            plot(T1,log(y1(1:length(T1)))-log(y1(1)),'linewidth',1);
            hold on;
            plot(T1,log(y1(length(T1)+1:end))-log(y1(length(T1)+1)),'linewidth',1);
            title('Skilled Worker')
    end
    xlabel('Year') 
    ylabel('Cumulative Log Change in Hourly Wages')
    xlim([1997 2020])
    ylim([-0.5 2])
    if (n==1)
        legend('Low lifetime automation exposure','High lifetime automation exposure','location','best')
    end 
end

filename = append(folder2,'/wage_trend.png'); 
exportgraphics(gcf,filename); 




