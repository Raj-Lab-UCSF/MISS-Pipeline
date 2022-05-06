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
matdir = '/Users/justintorok/Documents/MATLAB/MISS/MISS-MatFiles'; %define directory to draw from and save data to

%LOADING INPUT DATA
rnaseqdataname = 'tasic'; %define input scRNAseq dataset to use
if strcmp(rnaseqdataname,'tasic')
    load([matdir filesep 'Tasic_Inputs.mat'],'voxvgene','classkey','gene_names','genevct')
elseif strcmp(rnaseqdataname,'zeisel')
    load([matdir filesep 'Zeisel_Inputs.mat'],'voxvgene','classkey','gene_names','genevct')
end

%MRx3 GENE RANKING
lambda = 90; %percentile of genes to exclude from MRx3 ranking based on projection error added
resrate = 0.8;
nits = 100;
makenew = 0;
if makenew
    geneinds_bootstrap =  MRx3_Selector_Prefilter_Bootstrap(genevct,voxvgene,...
        size(voxvgene,2),lambda,resrate,nits,0); %generating MRx3 gene indices
    save([matdir filesep 'mrx3inds_bootstrapping.mat'],'geneinds_bootstrap');
else
    load([matdir filesep 'mrx3inds_bootstrapping.mat'],'geneinds_bootstrap');
end

%DEFINING MAPPING PARAMETER INPUTS
for i = 1:nits
    tic
    fprintf('Iteration %d/%d\n',i,nits);
    geneinds = geneinds_bootstrap(i,:);
%     ngen_param = 100:100:800;
    ngen_param = [26:20:806,850:50:1500,1600:100:2200,3849]; %coarser range
    missmethod = 'MRx3'; %gene ranking/subsetting method as a label, options are 'MRx3' and 'none'
    infmethod = 'inv+res';
    outstruct_temp = Cell_Density_Outstruct(genevct,voxvgene,gene_names,...
                        ngen_param,lambda,missmethod,infmethod,geneinds,matdir);
    elbowind = elbow_selector(outstruct_temp,0);
    outstruct_bootstrap(i) = outstruct_temp(elbowind);
    toc
end

%SAVING OUTPUT
save([matdir filesep 'outstruct_bootstrap.mat'],'outstruct_bootstrap','lambda',... %saving cell mapping & elbow output
    'missmethod','infmethod','geneinds_bootstrap','-v7.3');