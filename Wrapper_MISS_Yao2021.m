%This wrapper script runs through the MISS pipeline using the Yao, et al., 
% 2021 scRNAseq dataset, which includes multiple neocortical and
% hippocampal regions

%FILEPATH
matdir = '/Users/justintorok/Documents/MATLAB/MISS/MISS-MatFiles'; %define directory to draw from and save data to
addpath('/Users/justintorok/Documents/MATLAB/MISS/MISS-Pipeline/');
% matdir = '/data/rajlab1/user_data/justin/MatFiles'; %define directory to draw from and save data to
% addpath('/Home/jtorok/MISS-Pipeline/');

%LOADING INITIAL INPUT DATA
load([matdir filesep 'Yao_Inputs.mat'],'voxvgene','classkey','gene_names','genevct')

%MRx3 GENE RANKING USING scRNAseq FROM YAO, ET AL., 2021
lambda = 90;
load([matdir filesep 'Yao_MRx3_inds'],'geneinds'); %Tasic MRx3 gene indices

%DEFINING MAPPING PARAMETER INPUTS
ng_param_list = 400; %values of nG to test and map, going through genes in MRx3 ranked order
missmethod = 'MRx3'; %gene ranking/subsetting method as a label, options are 'MRx3' and 'none'
infmethod = 'inv+res'; %inversion method between E and C*D, options are 'inversion', 'inv+res' to also get residuals, and 'corr' for correlation mapping
savename = 'CellDensity_Yao_elbow.mat'; %set name of file to be saved

%GENERATING MAPS, RESIDUALS, METADATA, & ELBOW INDEX/CURVE 
outstruct = Cell_Density_Outstruct(genevct,voxvgene,... %getting outstruct of maps and residuals across nG
    gene_names,ng_param_list,lambda,missmethod,infmethod,...
    geneinds,matdir);              
save([matdir filesep savename],'outstruct','ng_param_list','lambda',... %saving cell mapping output
    'missmethod','infmethod','geneinds','classkey','-v7.3'); 

% %DEFINING ELBOW
% makefig = 1; %binary flag whether or not to make elbow curve plot
% elbowind = elbow_selector(outstruct,makefig); %getting elbow index value and generating elbow curve

% %MAIN TEXT FIGURES & RESULTS
% 
% %FIGURE 1 (PIPELINE SCHEMATIC)
% slicelocs = 34; %define what slice number to plot
% savenclose = 0; %binary flag for whether to save and close current figure window
% gnames = {'Slc17a6','Vip'}; %gene names to map (placeholder)
% Gene_Expression_Slice_Maps(gnames,slicelocs,savenclose,matdir); %panel STEP 1
% Cell_Type_tSNE('tasic',savenclose,matdir); %panel STEP 2
% %Panels STEP3-5 were made by hand in Keynote and Powerpoint using illustrations and simulated data
% 
% %LOAD COMAPRISON MAPPING
% nosubout = load([matdir filesep 'allnG_MRx3Prefilter_l90.mat'],'outstruct'); %loading no subsetting outstruct
% nosubout = nosubout.outstruct;
% corrout = load([matdir filesep 'CellDensity_corr_Tasic.mat'],'outstruct'); %loading correlation mapping outstruct
% corrout = corrout.outstruct;
% 
% %FIGURE 2 (GABAERGIC INTERNEURONS)
% savenclose = 0; %binary flag for whether to save and close current figure window
% %generates panel A
% region = 'neo'; %defining area over brain over which to make comparisons with data from Kim, et al., 2017
% Kim_Study_Comparison(outstruct,elbowind,region,savenclose,matdir) %generates panel B
% Kim_Study_Comparison(nosubout,1,region,savenclose,matdir) %generates panel C
% Kim_Study_Comparison(corrout,1,region,savenclose,matdir) %generates panel D
% Kim_Study_Comparison(corrout,2,region,savenclose,matdir) %generates panel E
% 
% %FIGURE 3 (NEOCORTICAL LAMINAR GLUTAMATERGIC NEURONS)
% mapmeths = {'MISS','No Subset Inv.','Subset Corr.'}; %giving titles to each generate panel
% slicelocs = [24 32]; %defining which slices to plot by number
% savenclose = 0; %binary flag for whether to save and close current figure window
% Figure3_taulayerslice(outstruct,elbowind,mapmeths{1},slicelocs,savenclose,matdir); %generates panel A, col 1
% % Figure3_taulayerslice(nosubout,1,mapmeths{2},slicelocs,savenclose,matdir); %generates panel A, col 2
% Figure3_taulayerslice(corrout,1,mapmeths{3},slicelocs,savenclose,matdir); %generates panel A, col 3
% % makenew = 1; %binary flag for loading in already calculated tau values (0) or calculating them anew (1)
% % if makenew %this if statement is for the calculation or loading of tau values
% %     [Rval_pv,Rval_sst,Rval_vip,tauvec] = rvalNtau_calc(outstruct,matdir); %calculating rvals and tau values
% % else
% %     load([matdir filesep 'rval_wholerange_Tasic.mat']); %loading in precalculated r values
% %     load([matdir filesep 'Tasic_MRx3Prefilter_fulltau.mat'],'tauvec'); %loading in precalculated tau values
% % end
% % rvalNtau_V_fronorm(Rval_pv,Rval_sst,Rval_vip,tauvec,outstruct); %generates panel B
% 
% %FIGURE 4 (GLIAL & ENDOTHELIAL CELLS)
% types = 22:25; %plot glia & vascular cells, 22 = astro, 23 = micro/macro, 24 = oligo, 25 = endo
% cmap_range = {[1 1 1;0.5 1 0.62],[1 0.75 0;1 0.25 0],[1 1 1;1 0.5 0.5],[0.75 0.25 0.4;1 0.25 0.1]}; %glia & endo colormaps
% xfac = 0.75; %Lines 87-88 are parameters for 3D brain illustrations
% voxthresh = 0.65;
% savenclose = 1; %binary flag for saving & closing (1) or opening (0) fig GUI window
% img_name = 'Tasic_G&E_'; %image name snippet for 3D renderings 
% %glia & endo 3D brain renderings (glia = Panels A-C, endo = Panel E)
% MISS_Brainframe(outstruct,types,'Tasic',elbowind,xfac,savenclose,voxthresh,cmap_range,img_name,matdir); 
% savenames = {'Astro','Micro','Oligo','Endo'};
% colorbar_creator(cmap_range,savenames) %making colorbars for glia (Panels A-C) & endo (Panel E) 3D renderings
% Figure_4d_glia(outstruct,elbowind,savenclose,matdir); %Panel D
% 
% %LOADING INPUT DATA USING scRNAseq FROM ZEISEL, ET AL., 2018
% load([matdir filesep 'Zeisel_Inputs.mat'],'voxvgene','classkey','gene_names','genevct')
% 
% %MRx3 GENE RANKING USING scRNAseq FROM ZEISEL, ET AL., 2018
% lambda = 90; %percentile of genes to exclude from MRx3 ranking based on projection error added
% makenew = 0; %binary flag for loading in already calculated MRx3 gene indices (0) or creating them anew (1)
% if makenew
%     geneinds =  MRx3_Selector_Prefilter(genevct,voxvgene,size(voxvgene,2),lambda,0); %generating MRx3 gene indices
% else
%     load([matdir filesep 'Zeisel_MRx3Inds.mat'],'geneinds'); %Tasic MRx3 gene indices
% end
% 
% %DEFINING MAPPING PARAMETER INPUTS
% ng_param_list = 201; %values of nG to test and map, going through genes in MRx3 ranked order
% missmethod = 'MRx3'; %gene ranking/subsetting method as a label, options are 'MRx3' and 'none'
% infmethod = 'inv+res'; %inversion method between E and C*D, options are 'inversion', 'inv+res' to also get residuals, and 'corr' for correlation mapping
% savename = 'CellDensity_Zeisel.mat'; %set name of file to be saved
% 
% %GENERATING MAPS, RESIDUALS, METADATA, & ELBOW INDEX/CURVE
% makenew = 0; %binary flag for loading in already calculated outstruct (0) or creating it anew (1)
% if makenew  
%     outstruct = Cell_Density_Outstruct(genevct,voxvgene,... %getting outstruct of maps and residuals across nG
%         gene_names,ng_param_list,lambda,missmethod,infmethod,...
%         geneinds,matdir);              
% %     save([matdir filesep savename],'outstruct','ng_param_list','lambda',... %saving cell mapping output
% %         'missmethod','infmethod','geneinds','-v7.3'); 
% else
%     load([matdir filesep 'Zeisel_outstruct_mod.mat'],'outstruct'); %Tasic, et al. outstruct
% end
% 
% %DEFINING ELBOW
% makefig = 1; %binary flag whether or not to make elbow curve plot
% elbowind = elbow_selector(outstruct,makefig); %getting elbow index value and generating elbow curve
% 
% %FIGURE 5 (MAPPING USING scRNAseq FROM ZEISEL, ET AL., 2018)
% typeinds = [153 27 83 117]; %cell type indices for slice mapping (minus one, find which one!!!)
% slicelocs = {[35 45 60];[35 40 65];[40 45 65];[30 50 65]}; %slice locations per cell type (same here!!!)
% savenclose = 0; %binary flag for saving as above
% Figure_5abcd_zeisel_typemaps(outstruct,elbowind,typeinds,slicelocs,savenclose,matdir); %Panels A-D
% Figure_5e_MISS_correlations_comp(outstruct,elbowind,savenclose,matdir); %Panel E
% types = [153 27 83 117]; %types to 3D render using scRNAseq from Zeisel, et al., 2018, Panel F
% cmap_range = {[1 1 0;0 1 0],[1 1 0;0 1 0],[1 1 0;0 1 0],[1 1 0;0 1 0]}; %3D rendering colormaps
% xfac = 0.75; %Lines 139-40 are parameters for 3D brain illustrations
% voxthresh = 0.65;
% savenclose = 1; %binary flag for saving as above
% img_name = 'Zeisel_'; %image name snippet for saving 3D renderings
% %3D brain renderings using scRNAseq data from Zeisel, et al., 2018, Panel F
% MISS_Brainframe(outstruct,types,'Zeisel',elbowind,xfac,savenclose,voxthresh,cmap_range,img_name,matdir);
% savenames = {'Zt153','Zt27','Zt83','Zt117'};
% colorbar_creator(cmap_range,savenames); %making colorbars for Zeisel, et al., 2018 maps (Panel F)
% 
