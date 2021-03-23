function tx_per_CT(genevct,classkey)
%This function creates the Number of Unique and Number of Total Transcripts
%per Major Cell Type from Tasic, et al., 2018 Supplementary Panel.

%Calculating unique and total transcripts per major cell type label
transcript_sums = sum(genevct);
ctgen_binary = genevct;
ctgen_binary(ctgen_binary>0) = 1;
unique_genes = sum(ctgen_binary);

%Plotting total transcripts per major cell type
figure;
subplot(1,2,1);
bar(transcript_sums,'LineWidth',1.5)
xticks(1:25);
yticks([ceil(max(transcript_sums)/3) ceil((max(transcript_sums)/3)*2) ceil(max(transcript_sums))]);
xticklabels(classkey);
xtickangle(90);
title('Average \Sigma_{Transcripts} per Cell Type','FontSize',18);
set(gca,'FontSize',16);

%Plotting number of unique transcripts per cell type
subplot(1,2,2);
bar(unique_genes,'LineWidth',1.5)
xticks(1:25);
yticks([ceil(max(unique_genes)/3) ceil((max(unique_genes)/3)*2) ceil(max(unique_genes))]);
xticklabels(classkey);
xtickangle(90);
title('Average Unique Transcripts per Cell Type','FontSize',18);
set(gca,'FontSize',16);

end