%This wrapper script runs through the basics of how to use the MISS
%pipeline to generate maps and ultimately select an elbow index and then
%using that choose an optimal mapm of cell types. At top, the user must
%define their filepath and then the dataset name. We only include 2
%datasets, 'tasic' & 'zeisel', so the user must initially choose from these
%options. If the user wants to add another option, they can trivially do so
%by adding an extra elseif statement to the if statements below. They
%should then change 'rnaseqdataname' to match whatever moniker they apply
%to their chosen dataset.

%FILEPATH
matdir = '/Users/christophermezias/Documents/MISS-MatFiles'; %define directory to draw from and save data to

%LOADING INPUT DATA
rnaseqdataname = 'tasic'; %define input scRNAseq dataset to use
if strcmp(rnaseqdataname,'tasic')
    load([matdir filesep 'Tasic_Inputs.mat'],'voxvgene','classkey','gene_names','genevct')
elseif strcmp(rnaseqdataname,'zeisel')
    load([matdir filesep 'Zeisel_Inputs.mat'],'voxvgene','classkey','gene_names','genevct')
end

%MRx3 GENE RANKING
lambda = 90; %percentile of genes to exclude from MRx3 ranking based on projection error added
makenew = 0; %binary flag for loading in already calculated MRx3 gene indices (0) or creating them anew (1)
if makenew
    geneinds =  MRx3_Selector_Prefilter(genevct,voxvgene,size(voxvgene,2),lambda,0); %generating MRx3 gene indices
elseif strcmp(rnaseqdataname,'tasic')
    load([matdir filesep 'MRx3_L90_inds'],'geneinds'); %Tasic MRx3 gene indices
elseif strcmp(rnaseqdataname,'zeisel')
    load([matdir filesep 'Zeisel_MRx3Inds.mat'],'geneinds'); %Zeisel MRx3 gene indices
end

%DEFINING MAPPING PARAMETER INPUTS
ng_param_list = 26:3855; %values of nG to test and map, going through genes in MRx3 ranked order
missmethod = 'MRx3'; %gene ranking/subsetting method as a label, options are 'MRx3' and 'none'
infmethod = 'inv+res'; %inversion method between E and C*D, options are 'inversion', 'inv+res' to also get residuals, and 'corr' for correlation mapping
savename = 'MISSMaps_2BSaved.mat'; %set name of file to be saved

%GENERATING MAPS, RESIDUALS, METADATA, & ELBOW INDEX/CURVE
makenew = 0; %binary flag for loading in already calculated outstruct (0) or creating it anew (1)
if makenew  
    outstruct = Cell_Density_Outstruct(genevct,voxvgene,... %getting outstruct of maps and residuals across nG
        gene_names,ng_param_list,lambda,missmethod,infmethod,...
        geneinds,matdir);               
elseif strcmp(rnaseqdataname,'tasic')
    load([matdir filesep 'Tasic_outstruct.mat'],'outstruct'); %Tasic, et al. outstruct
elseif strcmp(rnaseqdataname,'zeisel')
    load([matdir filesep 'Zeisel_outstruct.mat'],'outstruct'); %Zeisel, et al. outstruct
end

%DEFINING ELBOW
makefig = 1; %binary flag whether or not to make elbow curve plot
elbowind = elbow_selector(outstruct,makefig); %getting elbow index value and generating elbow curve

%SAVING OUTPUT
save([matdir filesep savename],'outstruct','ng_param_list','lambda',... %saving cell mapping & elbow output
    'missmethod','infmethod','geneinds','elbowind','-v7.3');