%This wrapper script runs through the figures and results from the MISS
%manuscript, and regenerates them in order, beginning with the basic
%pipeline results using the AIBS human scRANseq dataset

%FILEPATH
matdir = '/Users/justintorok/Documents/MATLAB/MISS/Human-related/RawDataFiles'; %define directory to draw from and save data to

%LOADING INITIAL INPUT DATA
load([matdir filesep 'AIBS_Human_scRNAseq_datastruct_subclass.mat'],'subclass_struct');
load([matdir filesep 'abagen426_1mm_datastruct.mat'],'regexpr_struct');
genevct = subclass_struct.averages.genevct;
classkey = subclass_struct.cell_types;
gene_names = subclass_struct.gene_names;
regvgene = regexpr_struct.regvgene;
regvgene_sums = sum(regvgene);
regvgene_sums = repmat(regvgene_sums,size(regvgene,1),1);
regvgene = regvgene ./ regvgene_sums;
clear regexpr_struct subclass_struct

%MRx3 GENE RANKING USING scRNAseq FROM AIBS
lambda = 90; %percentile of genes to exclude from MRx3 ranking based on projection error added
geneinds =  MRx3_Selector_Prefilter(genevct,regvgene,size(regvgene,2),lambda,0); %generating MRx3 gene indices
save([matdir filesep 'MRx3_l90_human_geneinds.mat'],'geneinds');

%DEFINING MAPPING PARAMETER INPUTS
ng_param_list = [100:100:1000, size(regvgene,2)]; %values of nG to test and map, going through genes in MRx3 ranked order
% ng_param_list = [50:10:1000, 1050:50:2000, 2100:100:3900, 4000:500:14500, size(regvgene,2)];
missmethod = 'MRx3'; %gene ranking/subsetting method as a label, options are 'MRx3' and 'none'
infmethod = 'inv+res'; %inversion method between E and C*D, options are 'inversion', 'inv+res' to also get residuals, and 'corr' for correlation mapping
savename = 'CellDensity_Human_allnG.mat'; %set name of file to be saved

%GENERATING MAPS, RESIDUALS, METADATA, & ELBOW INDEX/CURVE
outstruct = Cell_Density_Outstruct_Human(genevct,regvgene,... 
    gene_names,ng_param_list,lambda,missmethod,infmethod,...
    [],matdir);  %getting outstruct of maps and residuals across nG            

%DEFINING ELBOW
makefig = 1; %binary flag whether or not to make elbow curve plot
elbowind = elbow_selector_human(outstruct,makefig); %getting elbow index value and generating elbow curve

save([matdir filesep savename],'outstruct','ng_param_list','lambda',... %saving cell mapping output
        'missmethod','infmethod','geneinds','elbowind'); 