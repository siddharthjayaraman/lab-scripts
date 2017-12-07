%% Load Histogram data

% pacbio.reads lists all the read headers with their length in pacbio.readlength

% blastresult.reads list all the reads with length between 1400-2000bp and
% blastresult.alignment column 1: Read length / column 2: Alignment length
% / column 3: Alignmnet percentage

load('PacBio_manuscript_histogram_data.mat');

%% Used matlab histogram command for figures
%figure
%title('Unfiltered read length distribution')
rdl = pacbio.readlength;
rdl(rdl>2500) = 2500;
h = histogram(rdl,100,'linewidth',.7);
%line([1400 2000],[6.7e4 6.7e4],'color','k','linewidth',2)
rectangle('position',[1400 0 600 8.3e4],'linewidth',2.5)
box off
set(gca,'fontsize',14)
xlabel('Read length (bp)','fontsize',16)
ylabel('Number of reads','fontsize',16)

%%
figure(3)
%title('Alignment-length-coverage');
algn = blastresult.alignment(:,3);
algn(algn>100) = 100;
h = histogram(algn,'linewidth',.7);
rectangle('position',[60 0 40 9.5e4],'linewidth',2.5)
box off
set(gca,'fontsize',14)
xlabel('Alignment (%)','fontsize',16)
ylabel('Number of reads','fontsize',16)
