%% load data
load VSG_num_reads.mat

%% total number of reads for each vsg
figure(1)
nr = data(:,5);
h = histogram(log10(nr),'normalization','count');
set(gca,'yscale','log')
axis tight
x = .5*(h.BinEdges(1:end-2)' + h.BinEdges(2:end-1)');
y = h.Values(1:end-1)';
% regression
p = polyfit(x,log10(y),1);
x2 = [0.12; 5.25];
y2 = p(1)*x2+p(2);
hold on
plot(x2,10.^y2,'k-','linewidth',1.5);
hold off
box off
lbl = cell(5,1);
for k=1:5
    lbl{k} = sprintf('10^%d',k);
end
set(gca,'xtick',1:5,'xticklabel',lbl,'fontsize',14)
xlabel('Number of reads (n)','fontsize',16)
ylabel('Number of VSGs','fontsize',16)
% text(3.5,100,sprintf('N_{VSG}(n) \\sim %1.0f n^{%1.2f}',10^p(2), p(1)), ...
%     'fontsize',16);
text(3.7,120,sprintf('N_{VSG}(n) \\sim n^{%1.2f}', p(1)), ...
    'fontsize',16);
%%
figure(2)
map = colormap(lines);
bar(10.^(x),y)
% plot(10.^(x),y,'o','markersize',10,'markeredgecolor','k', ...
%     'markerfacecolor',map(2,:));
hold on
% regression
p = polyfit(x,log10(y),1);
x2 = [0.2; 5.4];
y2 = p(1)*x2+p(2); 
plot(10.^x2,10.^y2,'k-','linewidth',1.5);
hold off
%axis([1 5e5 5e-4 2])
box off
lbl = cell(5,1);
for k=1:5
    lbl{k} = sprintf('10^%d',k);
end
set(gca,'yscale','log','xscale','log','xtick',10.^(1:5),'xticklabel',lbl,'fontsize',14)
xlabel('Number of reads','fontsize',16)
ylabel('1-CDF','fontsize',16)
text(5e3,1,sprintf('P(x > n) \\sim n^{%1.2f}',p(1)), ...
    'fontsize',16);

