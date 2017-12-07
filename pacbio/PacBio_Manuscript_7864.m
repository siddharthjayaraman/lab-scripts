%% Pacbio data analysis and figures script


%%Import reads and number of passes

[reads.header,reads.sequence]=fastaread('reads_of_insert_job016520.fasta');
reads.header=reads.header';
reads.sequence=reads.sequence';
nop=importdata('number_of_passes_output_job016520.txt');

% Read length
for i=1:length(reads.header)
    reads.length(i,1)=length(reads.sequence{i});
end

% Read headers
reads.header=strrep(reads.header,'m140507_195301_42215_c100646022550000001823121309101432_s1_p0','balbc_10_0');
reads.header=strrep(reads.header,'m141217_202019_42215_c100735622550000001823161406051586_s1_p0','balbc_10_1');
reads.header=strrep(reads.header,'m141217_233935_42215_c100735622550000001823161406051587_s1_p0','balbc_10_2');
reads.header=strrep(reads.header,'m141218_025845_42215_c100735732550000001823161406051540_s1_p0','balbc_10_4');
reads.header=strrep(reads.header,'m141218_061758_42215_c100735732550000001823161406051541_s1_p0','balbc_10_5');
reads.header=strrep(reads.header,'m141218_093713_42215_c100735732550000001823161406051542_s1_p0','balbc_12_1');
reads.header=strrep(reads.header,'m141218_130029_42215_c100735732550000001823161406051543_s1_p0','balbc_12_2');
reads.header=strrep(reads.header,'m141218_161838_42215_c100735732550000001823161406051544_s1_p0','balbc_12_3');
reads.header=strrep(reads.header,'m150105_171452_42215_c100721762550000001823136804301585_s1_p0','balbc_12_4');
reads.header=strrep(reads.header,'m150105_203358_42215_c100721762550000001823136804301586_s1_p0','balbc_12_5');
reads.header=strrep(reads.header,'m140507_131439_42215_c100646022550000001823121309101430_s1_p0','balbc_3_0');
reads.header=strrep(reads.header,'m141213_100224_42215_c100735072550000001823161406051517_s1_p0','balbc_3_1');
reads.header=strrep(reads.header,'m141213_132137_42215_c100735282550000001823161406051540_s1_p0','balbc_3_2');
reads.header=strrep(reads.header,'m141213_164050_42215_c100735282550000001823161406051541_s1_p0','balbc_3_3');
reads.header=strrep(reads.header,'m141213_200003_42215_c100735282550000001823161406051542_s1_p0','balbc_3_4');
reads.header=strrep(reads.header,'m140507_163348_42215_c100646022550000001823121309101431_s1_p0','balbc_6_0');
reads.header=strrep(reads.header,'m141213_231916_42215_c100735282550000001823161406051543_s1_p0','balbc_6_1');
reads.header=strrep(reads.header,'m141214_023829_42215_c100735282550000001823161406051544_s1_p0','balbc_6_2');
reads.header=strrep(reads.header,'m141214_060128_42215_c100735282550000001823161406051545_s1_p0','balbc_6_4');
reads.header=strrep(reads.header,'m141217_170335_42215_c100735622550000001823161406051585_s1_p0','balbc_6_5');

nop.rowheaders=strrep(nop.rowheaders,'m140507_195301_42215_c100646022550000001823121309101432_s1_p0','balbc_10_0');
nop.rowheaders=strrep(nop.rowheaders,'m141217_202019_42215_c100735622550000001823161406051586_s1_p0','balbc_10_1');
nop.rowheaders=strrep(nop.rowheaders,'m141217_233935_42215_c100735622550000001823161406051587_s1_p0','balbc_10_2');
nop.rowheaders=strrep(nop.rowheaders,'m141218_025845_42215_c100735732550000001823161406051540_s1_p0','balbc_10_4');
nop.rowheaders=strrep(nop.rowheaders,'m141218_061758_42215_c100735732550000001823161406051541_s1_p0','balbc_10_5');
nop.rowheaders=strrep(nop.rowheaders,'m141218_093713_42215_c100735732550000001823161406051542_s1_p0','balbc_12_1');
nop.rowheaders=strrep(nop.rowheaders,'m141218_130029_42215_c100735732550000001823161406051543_s1_p0','balbc_12_2');
nop.rowheaders=strrep(nop.rowheaders,'m141218_161838_42215_c100735732550000001823161406051544_s1_p0','balbc_12_3');
nop.rowheaders=strrep(nop.rowheaders,'m150105_171452_42215_c100721762550000001823136804301585_s1_p0','balbc_12_4');
nop.rowheaders=strrep(nop.rowheaders,'m150105_203358_42215_c100721762550000001823136804301586_s1_p0','balbc_12_5');
nop.rowheaders=strrep(nop.rowheaders,'m140507_131439_42215_c100646022550000001823121309101430_s1_p0','balbc_3_0');
nop.rowheaders=strrep(nop.rowheaders,'m141213_100224_42215_c100735072550000001823161406051517_s1_p0','balbc_3_1');
nop.rowheaders=strrep(nop.rowheaders,'m141213_132137_42215_c100735282550000001823161406051540_s1_p0','balbc_3_2');
nop.rowheaders=strrep(nop.rowheaders,'m141213_164050_42215_c100735282550000001823161406051541_s1_p0','balbc_3_3');
nop.rowheaders=strrep(nop.rowheaders,'m141213_200003_42215_c100735282550000001823161406051542_s1_p0','balbc_3_4');
nop.rowheaders=strrep(nop.rowheaders,'m140507_163348_42215_c100646022550000001823121309101431_s1_p0','balbc_6_0');
nop.rowheaders=strrep(nop.rowheaders,'m141213_231916_42215_c100735282550000001823161406051543_s1_p0','balbc_6_1');
nop.rowheaders=strrep(nop.rowheaders,'m141214_023829_42215_c100735282550000001823161406051544_s1_p0','balbc_6_2');
nop.rowheaders=strrep(nop.rowheaders,'m141214_060128_42215_c100735282550000001823161406051545_s1_p0','balbc_6_4');
nop.rowheaders=strrep(nop.rowheaders,'m141217_170335_42215_c100735622550000001823161406051585_s1_p0','balbc_6_5');
nop.rowheaders=strrep(nop.rowheaders,'-','_');
