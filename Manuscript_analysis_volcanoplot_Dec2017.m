%% Load data and Sample groups
volcanodata=importdata('zscore.txt');
volcanodata.sgrp.samples=volcanodata.textdata(1,2:end);
[volcanodata.sgrp.mouse,rem]=strtok(volcanodata.sgrp.samples,'_');
[volcanodata.sgrp.tryp,rem]=strtok(rem,'_');
[volcanodata.sgrp.repl,rem]=strtok(rem,'_');
[volcanodata.sgrp.day,~]=strtok(rem,'_');

volcanodata.sgrp.day=strrep(volcanodata.sgrp.day,'D3','D03');
volcanodata.sgrp.day=strrep(volcanodata.sgrp.day,'D6','D06');
volcanodata.sgrp.day=strrep(volcanodata.sgrp.day,'2','D02');
volcanodata.sgrp.day=strrep(volcanodata.sgrp.day,'D1D02','D12');
volcanodata.sgrp.day=strrep(volcanodata.sgrp.day,'16','D16');

volcanodata.samples.samples=strcat(volcanodata.sgrp.mouse,'_',volcanodata.sgrp.tryp,'_',volcanodata.sgrp.day);
volcanodata.samples.legend=strcat(volcanodata.sgrp.mouse,'/',volcanodata.sgrp.tryp,' -',volcanodata.sgrp.day);



%% volcano 19 groups

[volcanodata.samples.uid,~,volcanodata.samples.uidx]=unique(volcanodata.samples.legend);
[volcanodata.samples.usid,~,volcanodata.samples.usidx]=unique(volcanodata.samples.samples);


for i=1:19
   volcanodata.(volcanodata.samples.usid{1,i}).data=[];
   volcanodata.(volcanodata.samples.usid{1,i}).samples={};
   volcanodata.(volcanodata.samples.usid{1,i}).legend={};
   for j=1:89
       if volcanodata.samples.usidx(j,1)==i
       volcanodata.(volcanodata.samples.usid{1,i}).data=[volcanodata.(volcanodata.samples.usid{1,i}).data,volcanodata.data(:,j)];
       volcanodata.(volcanodata.samples.usid{1,i}).samples=[volcanodata.(volcanodata.samples.usid{1,i}).samples,volcanodata.samples.samples(:,j)];
       volcanodata.(volcanodata.samples.usid{1,i}).legend=[volcanodata.(volcanodata.samples.usid{1,i}).legend,volcanodata.samples.legend(:,j)];
       end
   end
end
%% mattest Infected vs Uninfected [Early(D03) and Late(D10)]

%[PValues, TScores] = mattest(DataX, DataY)

[mymattest.C57_uninfected_927D3.P,mymattest.C57_uninfected_927D3.T]=mattest(volcanodata.C57_control_D02.data,volcanodata.C57_927_D03.data);
[mymattest.C57_uninfected_247D3.P,mymattest.C57_uninfected_247D3.T]=mattest(volcanodata.C57_control_D02.data,volcanodata.C57_247_D03.data);
[mymattest.C57_uninfected_927D10.P,mymattest.C57_uninfected_927D10.T]=mattest(volcanodata.C57_control_D02.data,volcanodata.C57_927_D10.data);
[mymattest.C57_uninfected_247D10.P,mymattest.C57_uninfected_247D10.T]=mattest(volcanodata.C57_control_D02.data,volcanodata.C57_247_D10.data);


[mymattest.BALBC_uninfected_927D3.P,mymattest.BALBC_uninfected_927D3.T]=mattest(volcanodata.BALBC_control_D02.data,volcanodata.BALBC_927_D03.data);
[mymattest.BALBC_uninfected_247D3.P,mymattest.BALBC_uninfected_247D3.T]=mattest(volcanodata.BALBC_control_D02.data,volcanodata.BALBC_247_D03.data);
[mymattest.BALBC_uninfected_927D10.P,mymattest.BALBC_uninfected_927D10.T]=mattest(volcanodata.BALBC_control_D02.data,volcanodata.BALBC_927_D10.data);
[mymattest.BALBC_uninfected_247D10.P,mymattest.BALBC_uninfected_247D10.T]=mattest(volcanodata.BALBC_control_D02.data,volcanodata.BALBC_247_D10.data);

%% mattest C57 vs BALBC and 927 vs 247 [Early(D03) and Late(D10)]

[mymattest.C57_BALBC_D3.P,mymattest.C57_BALBC_D3.T]=mattest([volcanodata.C57_927_D03.data,volcanodata.C57_247_D03.data],[volcanodata.BALBC_927_D03.data,volcanodata.BALBC_247_D03.data]);
[mymattest.C57_BALBC_D10.P,mymattest.C57_BALBC_D10.T]=mattest([volcanodata.C57_927_D10.data,volcanodata.C57_247_D10.data],[volcanodata.BALBC_927_D10.data,volcanodata.BALBC_247_D10.data]);

[mymattest.D3_927_247.P,mymattest.D3_927_247.T]=mattest([volcanodata.C57_927_D03.data,volcanodata.BALBC_927_D03.data],[volcanodata.C57_247_D03.data,volcanodata.BALBC_247_D03.data]);
[mymattest.D10_927_247.P,mymattest.D10_927_247.T]=mattest([volcanodata.C57_927_D10.data,volcanodata.BALBC_927_D10.data],[volcanodata.C57_247_D10.data,volcanodata.BALBC_247_D10.data]);

%% volcano plot infec vs uninfec

volcanodata.C57_uninfected_927D3=mavolcanoplot(volcanodata.C57_control_D02.data,volcanodata.C57_927_D03.data,mymattest.C57_uninfected_927D3.P,'PlotOnly', true,'Labels', fpkm.genes);
title('C57/Uninfected vs C57/927(Day-3)' ,'fontsize',12);
xlim([-3,3]);
ylim([0,5]);
print(gcf, '-dpng', '-r300', 'C57_uninfected_927D3.png');

close gcf

volcanodata.C57_uninfected_927D10=mavolcanoplot(volcanodata.C57_control_D02.data,volcanodata.C57_927_D10.data,mymattest.C57_uninfected_927D10.P,'PlotOnly', true,'Labels', fpkm.genes);
title('C57/Uninfected vs C57/927(Day-10)' ,'fontsize',12);
xlim([-3,3]);
ylim([0,5]);
print(gcf, '-dpng', '-r300', 'C57_uninfected_927D10.png');
close gcf

volcanodata.BALBC_uninfected_927D3=mavolcanoplot(volcanodata.BALBC_control_D02.data,volcanodata.BALBC_927_D03.data,mymattest.BALBC_uninfected_927D3.P,'PlotOnly', true,'Labels', fpkm.genes);
title('BALBc/Uninfected vs BALBc/927(Day-3)' ,'fontsize',12);
xlim([-3,3]);
ylim([0,5]);
print(gcf, '-dpng', '-r300', 'BALBC_uninfected_927D3.png');
close gcf

volcanodata.BALBC_uninfected_927D10=mavolcanoplot(volcanodata.BALBC_control_D02.data,volcanodata.BALBC_927_D10.data,mymattest.BALBC_uninfected_927D10.P,'PlotOnly', true,'Labels', fpkm.genes);
title('BALBc/Uninfected vs BALBc/927(Day-10)' ,'fontsize',12);
xlim([-3,3]);
ylim([0,5]);
print(gcf, '-dpng', '-r300', 'BALBC_uninfected_927D10.png');
close gcf

volcanodata.C57_uninfected_247D3=mavolcanoplot(volcanodata.C57_control_D02.data,volcanodata.C57_247_D03.data,mymattest.C57_uninfected_247D3.P,'PlotOnly', true,'Labels', fpkm.genes);
title('C57/Uninfected vs C57/247(Day-3)' ,'fontsize',12);
xlim([-3,3]);
ylim([0,5]);
print(gcf, '-dpng', '-r300', 'C57_uninfected_247D3.png');

close gcf

volcanodata.C57_uninfected_247D10=mavolcanoplot(volcanodata.C57_control_D02.data,volcanodata.C57_247_D10.data,mymattest.C57_uninfected_247D10.P,'PlotOnly', true,'Labels', fpkm.genes);
title('C57/Uninfected vs C57/247(Day-10)' ,'fontsize',12);
xlim([-3,3]);
ylim([0,5]);
print(gcf, '-dpng', '-r300', 'C57_uninfected_247D10.png');
close gcf

volcanodata.BALBC_uninfected_247D3=mavolcanoplot(volcanodata.BALBC_control_D02.data,volcanodata.BALBC_247_D03.data,mymattest.BALBC_uninfected_247D3.P,'PlotOnly', true,'Labels', fpkm.genes);
title('BALBc/Uninfected vs BALBc/247(Day-3)' ,'fontsize',12);
xlim([-3,3]);
ylim([0,5]);
print(gcf, '-dpng', '-r300', 'BALBC_uninfected_247D3.png');
close gcf

volcanodata.BALBC_uninfected_247D10=mavolcanoplot(volcanodata.BALBC_control_D02.data,volcanodata.BALBC_247_D10.data,mymattest.BALBC_uninfected_247D10.P,'PlotOnly', true,'Labels', fpkm.genes);
title('BALBc/Uninfected vs BALBc/247(Day-10)' ,'fontsize',12);
xlim([-3,3]);
ylim([0,5]);
print(gcf, '-dpng', '-r300', 'BALBC_uninfected_247D10.png');
close gcf

%% Volcano plots

% [mymattest.C57_BALBC_D3.P,mymattest.C57_BALBC_D3.T]=mattest([volcanodata.C57_927_D03.data,volcanodata.C57_247_D03.data],[volcanodata.BALBC_927_D03.data,volcanodata.BALBC_247_D03.data]);
% volcanodata.BALBC_uninfected_247D10=mavolcanoplot(volcanodata.BALBC_control_D02.data,volcanodata.BALBC_247_D10.data,mymattest.BALBC_uninfected_247D10.P,'PlotOnly', true,'Labels', fpkm.genes);

volcanodata.C57_BALBC_D3=mavolcanoplot([volcanodata.C57_927_D03.data,volcanodata.C57_247_D03.data],[volcanodata.BALBC_927_D03.data,volcanodata.BALBC_247_D03.data],mymattest.C57_BALBC_D3.P,'PlotOnly', true,'Labels', fpkm.genes);
title('C57 vs BALBc(Day-3)' ,'fontsize',12);
xlim([-3,3]);
ylim([0,5]);
print(gcf, '-dpng', '-r300', 'C57_BALBC_D3.png');
close gcf


volcanodata.C57_BALBC_D10=mavolcanoplot([volcanodata.C57_927_D10.data,volcanodata.C57_247_D10.data],[volcanodata.BALBC_927_D10.data,volcanodata.BALBC_247_D10.data],mymattest.C57_BALBC_D10.P,'PlotOnly', true,'Labels', fpkm.genes);
title('C57 vs BALBc(Day-10)' ,'fontsize',12);
xlim([-3,3]);
ylim([0,5]);
print(gcf, '-dpng', '-r300', 'C57_BALBC_D10.png');
close gcf

%[mymattest.D3_927_247.P,mymattest.D3_927_247.T]=mattest([volcanodata.C57_927_D03.data,volcanodata.BALBC_927_D03.data],[volcanodata.C57_247_D03.data,volcanodata.BALBC_247_D03.data]);
volcanodata.D3_927_247=mavolcanoplot([volcanodata.C57_927_D03.data,volcanodata.BALBC_927_D03.data],[volcanodata.C57_247_D03.data,volcanodata.BALBC_247_D03.data],mymattest.D3_927_247.P,'PlotOnly', true,'Labels', fpkm.genes);
title('Tryp927 vs Tryp247 (Day-3)' ,'fontsize',12);
xlim([-3,3]);
ylim([0,5]);
print(gcf, '-dpng', '-r300', 'D3_927_247.png');
close gcf

volcanodata.D10_927_247=mavolcanoplot([volcanodata.C57_927_D10.data,volcanodata.BALBC_927_D10.data],[volcanodata.C57_247_D10.data,volcanodata.BALBC_247_D10.data],mymattest.D10_927_247.P,'PlotOnly', true,'Labels', fpkm.genes);
title('Tryp927 vs Tryp247 (Day-10)' ,'fontsize',12);
xlim([-3,3]);
ylim([0,5]);
print(gcf, '-dpng', '-r300', 'D10_927_247.png');
close gcf
