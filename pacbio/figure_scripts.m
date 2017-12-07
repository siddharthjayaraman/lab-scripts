%% PacBio manyscript figure script 13-Apr-2017

load('figuredata_v2.mat');

%figuredata.ccsperbase.data=blastout.stranded60.topvariant.ccsERbaseOC;
%figuredata.ccsperbase.label=blastout.stranded60.topvariant.CCSxOC;
%figuredata.nofp.data=ORFdata;
%figuredata.backgrounderror.data=blastout.stranded60.topvariant.ccsERprob3;
%figuredata.aln.aln=blastout.stranded60.topvariant.alnzoneT;
%figuredata.aln.ins=blastout.stranded60.topvariant.inszoneT;
%figuredata.aln.del=blastout.stranded60.topvariant.delzoneT;
%figuredata.aln.mis=blastout.stranded60.topvariant.miszoneT;
%figuredata.variant.all=[];
%figuredata.variant.nosingleton=[];


%% Per Base CCS Figure
%figure
%errorbar(figuredata.ccsperbase.data(:,1),figuredata.ccsperbase.data(:,2),'rx')
yyaxis left
bp = boxplot(figuredata.ccsperbase.alldata(2:17,:)','colors','k','symbol','.',...
    'outliersize',10,'widths',.8);
for k=1:numel(bp)
    set(bp(k),'linewidth',1.5);
end
set(gca,'XTick',1:16,'XTickLabel',figuredata.ccsperbase.label(2:17),...
    'fontsize',14);
xlabel('Number of full passes','fontsize',16);
ylabel('Percentage error','fontsize',16);
set(gca,'YColor','black');
yyaxis right
npos = sum(figuredata.ccsperbase.alldata>0,2);
plot(1:16,npos(2:17),'Color','red','LineWidth',3);
set(gca,'YColor','red');
ylabel('Nucleotide positions with error','fontsize',16);
axis([0 17 -55 1240])
% ydata=num2str((figuredata.ccsperbase.data(:,1)),'%.2f');
% text((1:19),figuredata.ccsperbase.data(:,1),ydata, 'horizontal','left', 'vertical','bottom','FontSize',8);
% saveas(gcf,'Base.CCSerror.rate.png','png');

%%  Number of full passes error figure
ORFrow=figuredata.ccsperbase.label(2:17,:);
figure
map = colormap(lines);
bp = bar(figuredata.nofp.data(:,[4,2]),'stacked');
bp(1).EdgeColor = 'k';
bp(1).LineWidth = 1.5;
bp(1).FaceColor = map(1,:);
bp(2).EdgeColor = 'k';
bp(2).LineWidth = 1.5;
bp(2).FaceColor = map(3,:);

legend({'No ORF', 'ORF'},'fontsize',14,'location','north')
set(gca,'ytick',0:1e4:4e4);
yyaxis right
plot(figuredata.nofp.data(:,3),'Color','red','LineWidth',3);
axis([0 17 0 50])
set(gca,'XTickLabelRotation',0,'YColor','red','XTick',1:16,'XTickLabel',ORFrow, ...
    'fontsize',14,'ytick',0:10:50);
%title({'ORF vs Reads (W.r.t. number of full passes)'; 'for Tb08.27P2.380 dominant variant'});
xlabel('Number of full passes','fontsize',16);
ylabel('Percentage ORF reads','fontsize',16);
yyaxis left
ylabel('Number of reads','fontsize',16);

%saveas(gcf,'ORF.Percentage.eps','epsc');



%% Alignment Indel fugure
scrsz = get(0,'ScreenSize');
map = colormap(lines);
figure('Position',[1 1 .75*scrsz(3) 0.75*scrsz(4)/2])
nrtot = max(figuredata.aln.aln); % total number of reads (?)
plot(figuredata.aln.aln/nrtot,'LineWidth',2,'Color','k');
hold on

xf = find(figuredata.aln.ins>0.01*nrtot); % filter less than 1% 
plot(xf,figuredata.aln.ins(xf)/nrtot,'.','markersize',18,'color',map(5,:))

xf = find(figuredata.aln.del>0.01*nrtot); % filter less than 1% 
plot(xf,figuredata.aln.del(xf)/nrtot,'.','markersize',18,'color',map(2,:))

xf = find(figuredata.aln.mis>0.01*nrtot); % filter less than 1% 
plot(xf,figuredata.aln.mis(xf)/nrtot,'.','markersize',18,'color',map(1,:))

hold off

set(gca,'fontsize',14,'ytick',0.1:0.1:1,'yticklabel',10:10:100,'box','off')
xlabel('Sequence position (bp)','fontsize',16);
ylabel('Percentage of reads','fontsize',16);

axis([0 length(figuredata.aln.aln) 0 1])
legend({'Alignment','Insertion','Deletion','Mismatch'},'fontsize',14)

%saveas(gcf,'Tb08.27P2.380.alignment.Apr2017.eps','epsc');

%axes('position',[0.5 0.5 0.3 0.2])

%%  Background error rate figure
scrsz = get(0,'ScreenSize');
map = colormap(lines);
figure('Position',[1 1 .75*scrsz(3) 0.75*scrsz(4)/2])
cut = 0.01;
xf = find(figuredata.backgrounderror.data(10,:)>cut);
plot(xf,figuredata.backgrounderror.data(10,xf),'.','markersize',18,'color',map(5,:));

hold on

xf = find(figuredata.backgrounderror.data(2,:)>cut);
plot(xf,figuredata.backgrounderror.data(2,xf),'.','markersize',18,'color',map(2,:));

xf = find(figuredata.backgrounderror.data(6,:)>cut);
plot(xf,figuredata.backgrounderror.data(6,xf),'.','markersize',18,'color',map(1,:));


% plot(figuredata.backgrounderror.data(6,:),'Color','blue');    
% plot(figuredata.backgrounderror.data(10,:),'Color','green'); 

hold off

%title('Background error rate w.r.t. Number of full passes');

set(gca,'fontsize',14,'ytick',1:1:10,'box','off')
xlabel('Sequence position (bp)','fontsize',16);
ylabel('Percentage of reads','fontsize',16);

axis([0 length(figuredata.aln.aln) 0 10])
legend({'Insertion','Deletion','Mismatch'},'fontsize',14)

%saveas(gcf,'CCS2.6.10_background error.png','png');

%%  Variant counts
figure

fd = sort(figuredata.variant.all); % sort for easier plotting of data points

plot(1+[-0.07,-0.2,0.2,0.07,0],fd(:,1),'o','markersize',10,'markerfacecolor',map(2,:),...
    'markeredgecolor','k');

hold on
plot(2+[-0.1,-0.2,0.2,0.1,0],fd(:,2),'o','markersize',10,'markerfacecolor',map(2,:),...
    'markeredgecolor','k');
plot(3+[0,-0.1,0,0.1,0],fd(:,3),'o','markersize',10,'markerfacecolor',map(2,:),...
    'markeredgecolor','k');
plot(4+[0,-0.1,0,0.1,0],fd(:,4),'o','markersize',10,'markerfacecolor',map(2,:),...
    'markeredgecolor','k');

bp = boxplot(fd,'colors','k','outliersize',1,...
    'symbol','w.');
set(bp,'linewidth',1.);

hold off
axis([0.5 4.5 0 250])
set(gca,'xticklabel',[3,6,10,12],'ytick',50:50:250,'fontsize',14,'box','off');
xlabel('Days post infection','fontsize',16);
ylabel('Number of variants','fontsize',16);

%saveas(gcf,'Number-of-variants.png','png');

%%  Variant counts without singletons
figure

fd = sort(figuredata.variant.nosingleton); % sort for easier plotting of data points

plot(1+[-0.07,-0.2,0.2,0.07,0],fd(:,1),'o','markersize',10,'markerfacecolor',map(2,:),...
    'markeredgecolor','k');

hold on
plot(2+[-0.1,-0.2,0.2,0.1,0],fd(:,2),'o','markersize',10,'markerfacecolor',map(2,:),...
    'markeredgecolor','k');
plot(3+[0,-0.1,0,0.1,0],fd(:,3),'o','markersize',10,'markerfacecolor',map(2,:),...
    'markeredgecolor','k');
plot(4+[0,-0.1,0,0.1,0],fd(:,4),'o','markersize',10,'markerfacecolor',map(2,:),...
    'markeredgecolor','k');
bp = boxplot(fd,'colors','k','outliersize',1,...
    'symbol','w.');
set(bp,'linewidth',1.);

hold off
axis([0.5 4.5 0 220])
set(gca,'xticklabel',[3,6,10,12],'ytick',50:50:200,'fontsize',14,'box','off');
xlabel('Days post infection','fontsize',16);
ylabel('Number of variants','fontsize',16);

%saveas(gcf,'Number-of-variants.nosingletons.png','png');

%% Combined variant counts
%figure
fd1 = sort(figuredata.variant.all);
fd2 = sort(figuredata.variant.nosingleton);
fd = [fd1 fd2];
grps = [1 1 1 1 2 2 2 2];
grpd = [1 2 3 4 1 2 3 4];
map = colormap(lines);
bp = boxplot(fd,{grpd; grps},'colors','k', ... %'colorgroup',[1 2 1 2 1 2 1 2], 'colors',map, ...
    'factorseparator', [1],'outliersize',1,'symbol','w.');
set(bp,'linewidth',1.);
hold on
plot(1+[-0.07,-0.2,0.2,0.07,0],fd1(:,1),'o','markersize',10,'markerfacecolor',map(1,:),...
    'markeredgecolor','k');
plot(2+[-0.07,-0.2,0.2,0.07,0],fd2(:,1),'o','markersize',10,'markerfacecolor',map(2,:),...
    'markeredgecolor','k');
plot(3+[-0.1,-0.2,0.2,0.1,0],fd1(:,2),'o','markersize',10,'markerfacecolor',map(1,:),...
    'markeredgecolor','k');
plot(4+[-0.1,-0.2,0.2,0.1,0],fd2(:,2),'o','markersize',10,'markerfacecolor',map(2,:),...
    'markeredgecolor','k');
plot(5+[0,-0.1,0,0.1,0],fd1(:,3),'o','markersize',10,'markerfacecolor',map(1,:),...
    'markeredgecolor','k');
plot(6+[0,-0.1,0,0.1,0],fd2(:,3),'o','markersize',10,'markerfacecolor',map(2,:),...
    'markeredgecolor','k');
plot(7+[0,-0.1,0,0.1,0],fd1(:,4),'o','markersize',10,'markerfacecolor',map(1,:),...
    'markeredgecolor','k');
plot(8+[0,-0.1,0,0.1,0],fd2(:,4),'o','markersize',10,'markerfacecolor',map(2,:),...
    'markeredgecolor','k');
hold off
legend({'With singletons','Without singletons'},'location','northwest')
axis([0.5 8.5 0 250])
set(gca,'xtick',1.5:2:7.5,'xticklabel',[3,6,10,12],'ytick',50:50:250,'fontsize',14,'box','off');
xlabel('Days post infection','fontsize',16);
ylabel('Number of variants','fontsize',16);