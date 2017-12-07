%% Import data and sample groups

dendrogramdata=importdata('zscore.txt');
dendrogramdata.sgrp.samples=mdsgene.textdata(1,2:end);
[dendrogramdata.sgrp.mouse,rem]=strtok(dendrogramdata.sgrp.samples,'_');
[dendrogramdata.sgrp.tryp,rem]=strtok(rem,'_');
[dendrogramdata.sgrp.repl,rem]=strtok(rem,'_');
[dendrogramdata.sgrp.day,~]=strtok(rem,'_');
dendrogramdata.sgrp.day=strrep(dendrogramdata.sgrp.day,'D3','D03');
dendrogramdata.sgrp.day=strrep(dendrogramdata.sgrp.day,'D6','D06');
dendrogramdata.sgrp.day=strrep(dendrogramdata.sgrp.day,'2','D02');
dendrogramdata.sgrp.day=strrep(dendrogramdata.sgrp.day,'D1D02','D12');
dendrogramdata.sgrp.day=strrep(dendrogramdata.sgrp.day,'16','D16');
dendrogramdata.mysamples=strcat(dendrogramdata.sgrp.mouse,'_',dendrogramdata.sgrp.tryp,'_',dendrogramdata.sgrp.day);
dendrogramdata.legend=strcat(dendrogramdata.sgrp.mouse,'/',dendrogramdata.sgrp.tryp,' -',dendrogramdata.sgrp.day);

%% Reorder samples and calculate replicate mean 
[dendrogramdata.samples.uid,~,dendrogramdata.samples.uidx]=unique(dendrogramdata.mysamples);
[dendrogramdata.samples.usid,~,dendrogramdata.samples.usidx]=unique(dendrogramdata.legend);

dendrogramdata.samples.data=[];
dendrogramdata.samples.samples={};
dendrogramdata.samples.legend={};
for i=1:19
   for j=1:89
       if dendrogramdata.samples.usidx(j,1)==i
       dendrogramdata.samples.data=[dendrogramdata.samples.data,dendrogramdata.data(:,j)];
       dendrogramdata.samples.samples=[dendrogramdata.samples.samples,dendrogramdata.mysamples(:,j)];
       dendrogramdata.samples.legend=[dendrogramdata.samples.legend,dendrogramdata.legend(:,j)];
       end
   end
end
for i=1:19
       dendrogramdata.samples.Avdata(:,i)=mean(dendrogramdata.data(:,dendrogramdata.samples.usidx==i),2);
       dendrogramdata.samples.Avsamples{1,i}=dendrogramdata.samples.usid{i};
end

%% Calculate linkage and plot
dendrogramdata.dist=linkage(dendrogramdata.samples.Avdata','average');
dendrogram(dendrogramdata.dist,'Labels',dendrogramdata.samples.Avsamples,'Orientation','right','ColorThreshold',(0.8*max(dendrogramdata.dist(:,3))));
print(gcf, '-dpng', '-r300', 'dendrogram_sampleMEAN.png');
close gcf
