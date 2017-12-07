%% Pie charts
map = colormap(lines);
VSGdata=importdata('reads_of_insert_job016520_1400_2000_pivot_Table.tab');
VSGdata.data(isnan(VSGdata.data)) = 0;
VSGdata.pie.colour=brewermap(9,'Set1');
VSGdata.pie.colour([1 2],:)=VSGdata.pie.colour([2 1],:);
VSGdata.pie.colour(10:11,:)=[0,0,0;1,1,1];
colormap(VSGdata.pie.colour);
VSGdata.pie.explode=[1 1 1 1 1 1 1 1 1 1 1];
VSGdata.pie.label=strrep(VSGdata.textdata(1,2:end),'_','-');
VSGdata.pie.label2={'';'';'';'';'';'';'';'';'';'';''};

%% 
figure(1)
% colormap(VSGdata.pie.colour);
% for i=1:20
%     ax=subplot(4,5,i);
%     pie(ax,VSGdata.data(:,i),VSGdata.pie.explode,VSGdata.pie.label2);
%     title({VSGdata.pie.label{1,i};''})
% end
% %saveas(gcf,'piecharts_PacBioVSG','png');
% 
% %%
[n,m] = ind2sub([5 4],(1:20)');
dd = {'Day 3', 'Day 6', 'Day 10', 'Day 12'};
for k=1:4
    axes('position',[0 (4-k)*.25 .05 .25],...
        'xcolor',[1 1 1],'ycolor',[1 1 1])
    text(.5,.5,dd{k},'units','normalized','rotation',90,...
        'horizontalalignment','center','fontsize',20)
    box off
end
for i=1:20
   axes('position',[(n(i)-1)*.19+0.05 (4-m(i))*.25 .19 .25],...
       'xcolor',[1 1 1],'ycolor',[1 1 1])
   pie(VSGdata.data(:,i),VSGdata.pie.explode,VSGdata.pie.label2);
   % make sure colormap is correct when not all variants are present!
   tf = sum(VSGdata.data(:,i)>0,2)>0;
   colormap(VSGdata.pie.colour(tf,:));
   box off
end
%% legend in separate figure
figure(2)
% find one that has all 11 variables
i = find(sum(VSGdata.data>0)==11,1,'first');
pie(VSGdata.data(:,i),VSGdata.pie.explode,VSGdata.pie.label2);
   colormap(VSGdata.pie.colour);
box off
legend(VSGdata.textdata(2:end,1),'location','eastoutside','fontsize',14);
