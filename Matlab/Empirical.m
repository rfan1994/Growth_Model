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
    plot(T,coef(1,:),'linewidth',1);
    hold on;
    yline(0,'linewidth',1)
    patch([T fliplr(T)],[coef(1,:)-1.036*coef(2,:) fliplr(coef(1,:)+1.036*coef(2,:))],'k','FaceAlpha',0.1,'EdgeColor','w');
    patch([T fliplr(T)],[coef(1,:)-1.96*coef(2,:) fliplr(coef(1,:)+1.96*coef(2,:))],'k','FaceAlpha',0.05,'EdgeColor','w');
    xticks(1:6) 
    xlabel('Year')
    switch n
        case 1
            title('Share of skilled workers')
        case 2
            title('Labor composition growth')
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
    ylim([-1 2])
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
            	ylim([-0.5 2])
            case 2
            	ylim([-0.5 1.5])
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
        	legend('Low skill','Median skill','High skill')
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
filename = append(folder1,'/l_aoe');
table = readcell(filename);
aoe = cell2mat(table(3:end-1,4));
r = cell2mat(table(3:end-1,5));
l = cell2mat(table(3:end-1,6));

figure('position',[0,0,800,250])
ax1 = subplot(1,2,1);
scatter(aoe,r,'filled')
l1 = lsline(ax1);
set(l1,'linewidth',1.5,'Color','r')
xlim([0 100])
xlabel('Automation Exposure')
ylabel('Percentage')
title('Training-Working Ratio')

ax2 = subplot(1,2,2);
scatter(aoe,l,'filled')
l2 = lsline(ax2);
set(l2,'linewidth',1.5,'Color','r')
xlim([0 100])
xlabel('Automation Exposure')
ylabel('Percentage')
title('Training parcitipation rate')
filename = append(folder2,'/train.png');       
exportgraphics(gcf,filename)
close all;




