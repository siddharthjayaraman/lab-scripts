%% LOAD DATA
mdsgene=importdata('zscore.txt');


%% 2D MDS calculation

mdsgene.dist=pdist(mdsgene.data');
mdsgene.distsq=squareform(mdsgene.dist);
mdsgene.mds=mdscale(mdsgene.distsq,2);

%% Samples groups and colour coding
mdsgene.sgrp.samples=mdsgene.textdata(1,2:end);
[mdsgene.sgrp.mouse,rem]=strtok(mdsgene.sgrp.samples,'_');
[mdsgene.sgrp.tryp,rem]=strtok(rem,'_');
[mdsgene.sgrp.repl,rem]=strtok(rem,'_');
[mdsgene.sgrp.day,~]=strtok(rem,'_');

mdsgene.samples=strcat(mdsgene.sgrp.mouse,'_',mdsgene.sgrp.tryp,'_',mdsgene.sgrp.day);
mdsgene.legend=strcat(mdsgene.sgrp.mouse,'/',mdsgene.sgrp.tryp,' -',mdsgene.sgrp.day);


for i=1:89
    if strncmp(mdsgene.samples{1,i},'C57_control',11)
        mdsgene.colour(i,1:3)= [1 1 0.4];
    elseif strcmp(mdsgene.samples{1,i},'C57_927_D3')
        mdsgene.colour(i,1:3)= [1 0.8 0.4];
    elseif strcmp(mdsgene.samples{1,i},'C57_927_D6')
        mdsgene.colour(i,1:3)= [1 0.6 0.4];
    elseif strcmp(mdsgene.samples{1,i},'C57_927_D10')        
        mdsgene.colour(i,1:3)= [1 0.6 0];
        
    elseif strcmp(mdsgene.samples{1,i},'C57_247_D3')
        mdsgene.colour(i,1:3)= [0.6 1 0.6];
    elseif strcmp(mdsgene.samples{1,i},'C57_247_D6')
        mdsgene.colour(i,1:3)= [0.4 1 0.4];
    elseif strcmp(mdsgene.samples{1,i},'C57_247_D10')        
        mdsgene.colour(i,1:3)= [0.2 1 0.2];
    elseif strcmp(mdsgene.samples{1,i},'C57_247_D12')        
        mdsgene.colour(i,1:3)= [0 1 0];
        
    elseif strncmp(mdsgene.samples{1,i},'BALBC_control',13)
        mdsgene.colour(i,1:3)= [0 1 1];
    elseif strcmp(mdsgene.samples{1,i},'BALBC_927_D3')
        mdsgene.colour(i,1:3)= [1 0.8 1];
    elseif strcmp(mdsgene.samples{1,i},'BALBC_927_D6')
        mdsgene.colour(i,1:3)= [1 0.6 1];
    elseif strcmp(mdsgene.samples{1,i},'BALBC_927_D10')        
        mdsgene.colour(i,1:3)= [1 0.2 1];
    elseif strcmp(mdsgene.samples{1,i},'BALBC_927_D12')        
        mdsgene.colour(i,1:3)= [1 0 1];
        
    elseif strcmp(mdsgene.samples{1,i},'BALBC_247_D3')
        mdsgene.colour(i,1:3)= [0.8 0.8 1];
    elseif strcmp(mdsgene.samples{1,i},'BALBC_247_D6')
        mdsgene.colour(i,1:3)= [0.6 0.6 1];
    elseif strcmp(mdsgene.samples{1,i},'BALBC_247_D10')        
        mdsgene.colour(i,1:3)= [0.2 0.2 1];
    elseif strcmp(mdsgene.samples{1,i},'BALBC_247_D12')        
        mdsgene.colour(i,1:3)= [0 0 1];
    end
end

%% Scatter plot

%Plot legend
mdsgene.scatter.sample=mdsgene.legend;
mdsgene.scatter.sample(82:89)={'C57/Uninfected'};
mdsgene.scatter.sample(40:47)={'BALBC/Uninfected'};
mdsgene.scatter.sample=strrep(mdsgene.scatter.sample,'-D3','-D03');
mdsgene.scatter.sample=strrep(mdsgene.scatter.sample,'-D6','-D06');

%Ploting unique groups
[mdsgene.scatter.uid,~,mdsgene.scatter.uidx] = unique(mdsgene.scatter.sample);

%Print figure
figure('Position',[100,100,800,700]);
hold on;
box off;
for i=1:max(mdsgene.scatter.uidx)
    scatter(mdsgene.mds(mdsgene.scatter.uidx==i,1),mdsgene.mds(mdsgene.scatter.uidx==i,2),50,mdsgene.colour(mdsgene.scatter.uidx==i,:),'o','fill');
end
legend(mdsgene.scatter.uid,'Location','EastOutside');
title('MDS plot','fontsize',14);
set(gcf,'PaperUnits','centimeters','PaperSize',[20,10],'PaperPosition',[0 0 20 10]);
print(gcf, '-dpng', '-r300', 'MDS_Plot');



