function genelist = GeneList_Generator(study,elbowind,geneinds,directory)

if nargin < 1
    study = 'Tasic';
end

if nargin < 4
    directory = '/Users/christophermezias/Documents/MISS-MatFiles';
    if nargin < 3
        if strcmp(study,'Tasic')
            load([directory filesep 'MRx3_L90_inds.mat'],'geneinds');
        elseif strcmp(study,'Zeisel')
            load([directory filesep 'Zeisel_MRx3Inds.mat'],'geneinds');
        end
        if nargin < 2
            if strcmp(study,'Tasic')
                elbowind = 606;
            elseif strcmp(study,'Zeisel')
                elbowind = 1360;
            end
        end
    end
end

if strcmp(study,'Tasic')
    load([directory filesep 'Tasic_Inputs.mat'],'voxvgene','genevct','gene_names');
elseif strcmp(study,'Zeisel')
    load([directory filesep 'Zeisel_Inputs.mat'],'voxvgene','gene_names','genevct');
end        

%Generating gene list
elbow_t = elbowind;
preloadinds_t = geneinds(1:elbow_t);
genelist = gene_names(preloadinds_t);

%Gene Overlap Visualization
mrx3genes = genevct(preloadinds_t,:);
figure('Position',[0 0 1200 700]);
imagesc(mrx3genes.'); colormap(bone);
title(['Gene Overlap Among MRx3 Chosen Genes, ' study ', et al., 2018'],'FontSize',20);
xticks([]);
yticks([]);
ylabel('Cell Types','FontSize',20);
xlabel('scRNAseq Gene Expression','FontSize',20);
set(gca,'FontSize',20);
colorbar;

end
        