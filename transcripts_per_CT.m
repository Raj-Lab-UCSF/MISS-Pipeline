matdir = '/Users/christophermezias/Documents/MISS_General/MatFiles'; 
load([matdir filesep 'PresetInputs.mat']);
transcript_sums = sum(genevct);
ctgen_binary = genevct;
ctgen_binary(ctgen_binary>0) = 1;
unique_genes = sum(ctgen_binary);

figure;
subplot(1,2,1);
bar(transcript_sums,'LineWidth',1.5)
xticks(1:25);
yticks([ceil(max(transcript_sums)/3) ceil((max(transcript_sums)/3)*2) ceil(max(transcript_sums))]);
xticklabels(classkey);
xtickangle(90);
title('Average \Sigma_{Transcripts} per Cell Type','FontSize',18);
set(gca,'FontSize',16);

subplot(1,2,2);
bar(unique_genes,'LineWidth',1.5)
xticks(1:25);
yticks([ceil(max(unique_genes)/3) ceil((max(unique_genes)/3)*2) ceil(max(unique_genes))]);
xticklabels(classkey);
xtickangle(90);
title('Average Unique Transcripts per Cell Type','FontSize',18);
set(gca,'FontSize',16);

fname = 'lambda250_seed0_superfine.mat';
load([matdir filesep fname]);

idx = 235;

Bmeans = outstruct(idx).Bmeans;
Bsums = outstruct(idx).Bsums;

Bmeanmean = mean(Bmeans,1);
Bsumsum = sum(Bsums,1);
rval_meantot = corr(Bmeanmean.',transcript_sums.');
rval_meanuni = corr(Bmeanmean.',unique_genes.');
rval_sumtot = corr(Bsumsum.',transcript_sums.');
rval_sumuni = corr(Bsumsum.',unique_genes.');

mrx3_gennames = {}; %load copied MRx3 nG = 529 gene names here
selinds = ismember(gene_names,mrx3_gennames);
selgenes = genevct(selinds,:);
transcript_sums = sum(selgenes);
ctgen_binary = selgenes;
ctgen_binary(ctgen_binary>0) = 1;
unique_genes = sum(ctgen_binary);

rval_meantot_mrx3 = corr(Bmeanmean.',transcript_sums.');
rval_meanuni_mrx3 = corr(Bmeanmean.',unique_genes.');
rval_sumtot_mrx3 = corr(Bsumsum.',transcript_sums.');
rval_sumuni_mrx3 = corr(Bsumsum.',unique_genes.');


