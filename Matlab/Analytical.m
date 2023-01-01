clear; close all; clc;
folder2 = '/Users/rongfan/Desktop/Growth_Model/Latex';

figure('Position',[0 0 300 150]); 
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
legend({'\gamma_K','\gamma_L(i,0)','\gamma_H(i,0)'},'location','northwest')
title('Task productivity')
filename = fullfile(folder2,'gamma.png');
saveas(gcf,filename); 

%%
clear; close all; clc;
folder2 = '/Users/rongfan/Desktop/Growth_Model/Latex';

global A sigma B_NH B_NL b_h h_HL k L_L L_H I;
A = 0.12; sigma = 2; B_NH = 2.26; B_NL = 1; b_h = 1;
k = 2; L_L = 0.7; L_H = 0.3;

figure('position',[0,0,800,250])
for i = 1:4
    if (i==1)
        I = 0.2; h_HL = 0.1;
        I1 = I;
    elseif (i==2)
        I = 0.4; h_HL = 0.1;
        I2 = I;
    elseif (i==3)
        I = 0.2; h_HL = 0.1;
        I1 = I;
    elseif (i==4)
        I = 0.2; h_HL = 1;
        I2 = I;
    end

    myfun = @arbitrage;
    if (i==0||i==1||i==3)
        S1 = fzero(myfun,0.6);
        S = S1;
    elseif (i==2||i==4)
        S2 = fzero(myfun,0.6);
        S = S2;
    end
    
    Eta = I^(1d0/(sigma-1d0));
    gamma_HL = exp(B_NH-B_NL+b_h*h_HL);
    Gamma_H = ((1d0-exp(-B_NH*(1d0-S)*(sigma-1d0)))/(B_NH*(sigma-1d0)))^(1d0/(sigma-1d0));   
    Gamma_L = ((exp(-B_NL*(1d0-S)*(sigma-1d0))-exp(-B_NL*(1d0-I)*(sigma-1d0)))/(B_NL*(sigma-1d0)))^(1d0/(sigma-1d0))/gamma_HL;
    y = A*((Eta*k/exp(B_NH+b_h*h_HL))^((sigma-1d0)/sigma)+(Gamma_L*L_L)^((sigma-1d0)/sigma)+(Gamma_H*L_H)^((sigma-1d0)/sigma))^(sigma/(sigma-1d0));
    R = A*Eta*(y/(A*Eta*k/exp(B_NH+b_h*h_HL)))^(1d0/sigma);
    W_H = exp(B_NH+b_h*h_HL)*A*Gamma_H*(y/(A*Gamma_H*L_H))^(1d0/sigma);
    W_L = exp(B_NH+b_h*h_HL)*A*Gamma_L*(y/(A*Gamma_L*L_L))^(1d0/sigma);

    f1 = @(i) (R/A)^(1-sigma)*i;
    f2 = @(i) f1(I)+(W_L/A)^(1-sigma)*(exp(B_NL*(sigma-1)*i)-exp(B_NL*(sigma-1)*I))/(B_NL*(sigma-1));
    f3 = @(i) f2(S)+(W_H/A/exp(b_h*h_HL))^(1-sigma)*(exp(B_NH*(sigma-1)*i)-exp(B_NH*(sigma-1)*S))/(B_NH*(sigma-1));
    if (i==1||i==3)
        COLOR = [0 0.1 0.6; 0 0.5 0.6; 0 0.9 0.6];
        subplot(1,2,(i+1)/2)
    elseif (i==2||i==4)
        COLOR = [0.6 0 0.1; 0.6 0 0.5; 0.6 0 0.9];
    end
    
    if (i==1||i==3)
        h1 = fplot(f1,[0,I],'linewidth',1,'color',COLOR(1,:));
        hold on;
    elseif (i==2||i==4)
        h2 = fplot(f1,[0,I],'linewidth',1,'color',COLOR(1,:));
        hold on;
    end    
    line([I,I],[0,f1(I)],'linewidth',1,'color',COLOR(1,:),'lineStyle','--');
    hold on;
    fplot(f2,[I,S],'linewidth',1,'color',COLOR(2,:));
    hold on;
    line([S,S],[0,f2(S)],'linewidth',1,'color',COLOR(2,:),'lineStyle','--');
    hold on;
    fplot(f3,[S,1],'linewidth',1,'color',COLOR(3,:));
    hold on; 
    if (i==2)
        annotation('textbox',[0.2,0.5,0,0],'HorizontalAlignment','center','string','Capital')
        annotation('textbox',[0.3,0.7,0,0],'HorizontalAlignment','center','string','Low-skill labor')
        annotation('textbox',[0.4,0.9,0,0],'HorizontalAlignment','center','string','High-skill labor')
        xticks([0,I1,I2,S1,S2,1])
        xticklabels({'N-1','I','I''','S','S''','N'})
        xlabel('Task index')
        yticks([0 1])
        ylabel('Cumulative output share')
        title('Increase in automation')
        legend([h1 h2],'Benchmark','Increase in automation','location','northwest')
    elseif (i==4)
        annotation('textbox',[0.6,0.4,0,0],'HorizontalAlignment','center','string','Capital')
        annotation('textbox',[0.7,0.65,0,0],'HorizontalAlignment','center','string','Low-skill labor')
        annotation('textbox',[0.85,0.9,0,0],'HorizontalAlignment','center','string','High-skill labor')
        xticks([0,I,S2,S1,1])
        xticklabels({'N-1','I','S''','S','N'})
        xlabel('Task index')
        yticks([0 1])
        title('Increase in human capital')
        legend([h1 h2],'Benchmark','Increase in human capital','location','northwest')
    end
end
filename = fullfile(folder2,'comparative.png');
exportgraphics(gcf,filename)


function x = arbitrage(S)
    global A sigma B_NH B_NL b_h h_HL k L_L L_H I;
    Eta = I^(1d0/(sigma-1d0));
    gamma_HL = exp(B_NH-B_NL+b_h*h_HL);
    Gamma_H = ((1d0-exp(-B_NH*(1d0-S)*(sigma-1d0)))/(B_NH*(sigma-1d0)))^(1d0/(sigma-1d0));   
    Gamma_L = ((exp(-B_NL*(1d0-S)*(sigma-1d0))-exp(-B_NL*(1d0-I)*(sigma-1d0)))/(B_NL*(sigma-1d0)))^(1d0/(sigma-1d0))/gamma_HL;
    y = A*((Eta*k/exp(B_NH+b_h*h_HL))^((sigma-1d0)/sigma)+(Gamma_L*L_L)^((sigma-1d0)/sigma)+(Gamma_H*L_H)^((sigma-1d0)/sigma))^(sigma/(sigma-1d0));
    w_H = A*Gamma_H*(y/(A*Gamma_H*L_H))^(1d0/sigma);
    w_L = A*Gamma_L*(y/(A*Gamma_L*L_L))^(1d0/sigma);
    x = w_H/w_L-exp((B_NH-B_NL)*S+b_h*h_HL);
end

