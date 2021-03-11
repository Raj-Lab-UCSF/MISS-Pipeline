matdir = '/Users/christophermezias/Documents/MISS_General/MatFiles';
lambda = 250;

load([matdir filesep 'PresetInputs.mat'],'regvgene','genevct','gene_names');
voxvgene = regvgene;
preloadinds_t = MRx3_Selector(genevct,voxvgene,529,lambda);
tasic_genes = gene_names(preloadinds_t);

%Gene Overlap Visualization
mrx3genes = genevct(preloadinds_t,:);
figure('Position',[0 0 1200 700]);
imagesc(mrx3genes.'); colormap(bone);
title('Gene Overlap Among MRx3 Chosen Genes, Tasic, et al., 2018','FontSize',20);
xticks([]);
yticks([]);
ylabel('Cell Types','FontSize',20);
xlabel('scRNAseq Gene Expression','FontSize',20);
set(gca,'FontSize',20);
colorbar;

load([matdir filesep 'Zeisel_Inputs.mat'],'regvgene','entrez_names','meanexprmat');
voxvgene = regvgene;
genevct = meanexprmat.';
preloadinds_z = MRx3_Selector(genevct,voxvgene,1168,lambda);
zeisel_genes = gene_names(preloadinds_z);

mrx3genes = genevct(preloadinds_z,:);
figure('Position',[0 0 1200 700]);
clims = [min(genevct(:)) prctile(genevct(:),99)];
imagesc(mrx3genes.',clims); colormap(bone);
title('Gene Overlap Among MRx3 Chosen Genes, Zeisel, et al., 2018','FontSize',20);
xticks([]);
yticks([]);
ylabel('Cell Types','FontSize',20);
xlabel('scRNAseq Gene Expression','FontSize',20);
set(gca,'FontSize',20);
colorbar;
