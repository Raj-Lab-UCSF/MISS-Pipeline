%This wrapper script runs through the supplemental material from the MISS
%manuscript, and regenerates them in order, beginning with the basic
%pipeline results  and then proceeding to the supplementary figures,
%tables, and datasets. S. Table 1, however, was created as an explainer and
%was not generated using code, and so is not recreated here. 

%FILEPATH
matdir = '/Users/christophermezias/Documents/MISS-MatFiles'; %define directory to draw from and save data to

%S. TABLE 2
load([matdir filesep 'listB.mat'],'listB')
regnames = listB(:,1); %This is the region names list in S. Table 2

%LOADING TASIC, ET AL., 2018 INPUT DATA
load([matdir filesep 'Tasic_Inputs.mat'],'voxvgene','classkey','gene_names','genevct')

%MRx3 GENE RANKING USING scRNAseq FROM TASIC, ET AL., 2018
lambda = 90; %percentile of genes to exclude from MRx3 ranking based on projection error added
makenew = 0; %binary flag for loading in already calculated MRx3 gene indices (0) or creating them anew (1)
if makenew
    geneinds =  MRx3_Selector_Prefilter(genevct,voxvgene,size(voxvgene,2),lambda,0); %generating MRx3 gene indices
else
    load([matdir filesep 'MRx3_L90_inds'],'geneinds'); %Tasic MRx3 gene indices
end

%DEFINING MAPPING PARAMETER INPUTS
ng_param_list = 26:3855; %values of nG to test and map, going through genes in MRx3 ranked order
missmethod = 'MRx3'; %gene ranking/subsetting method as a label, options are 'MRx3' and 'none'
infmethod = 'inv+res'; %inversion method between E and C*D, options are 'inversion', 'inv+res' to also get residuals, and 'corr' for correlation mapping
savename = 'CellDensity_Tasic.mat'; %set name of file to be saved

%GENERATING MAPS, RESIDUALS, METADATA, & ELBOW INDEX/CURVE
makenew = 0; %binary flag for loading in already calculated outstruct (0) or creating it anew (1)
if makenew  
    outstruct = Cell_Density_Outstruct(genevct,voxvgene,... %getting outstruct of maps and residuals across nG
        gene_names,ng_param_list,lambda,missmethod,infmethod,...
        geneinds,matdir);              
    save([matdir filesep savename],'outstruct','ng_param_list','lambda',... %saving cell mapping output
        'missmethod','infmethod','geneinds','-v7.3'); 
else
    load([matdir filesep 'Tasic_outstruct.mat'],'outstruct'); %Tasic, et al. outstruct
end

%S. FIGURE 1, RESULTS USING scRNAseq FROM TASIC, ET AL., 2018
makefig = 1; %binary flag whether or not to make elbow curve plot
%S. Figure 1a was created as an illustration, partially by hand, and so is not easily recreatable via code here
elbowind = elbow_selector(outstruct,makefig); %getting elbow index value and elbow curve generated is S.Figures 1b

%S. TABLE 4
STable4_CT = outstruct(elbowind).Bmeans; %This output is S.Table 4 Cell Type per Region values, scRNAseq from Tasic, et al., 2018

%S. FIGURE 2
savenclose = 1; %binary flag for whether to save and close or open figure GUI
nG = outstruct(elbowind).nGen;
S_Figure_Residuals(nG,lambda,geneinds,savenclose,matdir); %Generates S. Figures 2a-b
naround = 100; %parameter for how many genes around elbow to include in S. Figure 2c
SFigure2c_Generator(elbowind,outstruct,ng_param_list,naround); %Generates S. Figure 2c

%S FIGURE & TABLE 3, RESULTS USING scRNAseq FROM TASIC, ET AL., 2018
genelist = GeneList_Generator('Tasic',elbowind,geneinds,matdir); %Generates S. Figure 3a & S.Table 3, Sheet 1
tx_per_CT(genevct,classkey); %Generates S. Figures 3c & d

%LOADING ZEISEL, ET AL., 2018 INPUT DATA
load([matdir filesep 'Zeisel_Inputs.mat'],'voxvgene','classkey','gene_names','genevct')

%MRx3 GENE RANKING USING scRNAseq FROM ZEISEL, ET AL., 2018
lambda = 90; %percentile of genes to exclude from MRx3 ranking based on projection error added
makenew = 0; %binary flag for loading in already calculated MRx3 gene indices (0) or creating them anew (1)
if makenew
    geneinds =  MRx3_Selector_Prefilter(genevct,voxvgene,size(voxvgene,2),lambda,0); %generating MRx3 gene indices
else
    load([matdir filesep 'Zeisel_MRx3Inds.mat'],'geneinds'); %Tasic MRx3 gene indices
end

%DEFINING MAPPING PARAMETER INPUTS
ng_param_list = 201:3803; %values of nG to test and map, going through genes in MRx3 ranked order
missmethod = 'MRx3'; %gene ranking/subsetting method as a label, options are 'MRx3' and 'none'
infmethod = 'inv+res'; %inversion method between E and C*D, options are 'inversion', 'inv+res' to also get residuals, and 'corr' for correlation mapping
savename = 'CellDensity_Zeisel.mat'; %set name of file to be saved

%GENERATING MAPS, RESIDUALS, METADATA, & ELBOW INDEX/CURVE
makenew = 0; %binary flag for loading in already calculated outstruct (0) or creating it anew (1)
if makenew  
    outstruct = Cell_Density_Outstruct(genevct,voxvgene,... %getting outstruct of maps and residuals across nG
        gene_names,ng_param_list,lambda,missmethod,infmethod,...
        geneinds,matdir);              
    save([matdir filesep savename],'outstruct','ng_param_list','lambda',... %saving cell mapping output
        'missmethod','infmethod','geneinds','-v7.3'); 
else
    load([matdir filesep 'Zeisel_outstruct.mat'],'outstruct'); %Tasic, et al. outstruct
end

%S. FIGURE 1, RESULTS USING scRNAseq FROM ZEISEL, ET AL., 2018
makefig = 1; %binary flag whether or not to make elbow curve plot
elbowind = elbow_selector(outstruct,makefig); %getting elbow index value and elbow curve generated is S.Figures 1c

%S. TABLE 5
STable5_CT = outstruct(elbowind).Bmeans; %This output is S.Table 5 Cell Type per Region values, scRNAseq from Zeisel, et al., 2018

%S FIGURE & TABLE 3, RESULTS USING scRNAseq FROM ZEISEL, ET AL., 2018
genelist = GeneList_Generator('Zeisel',outstruct(elbowind).nGen,geneinds,matdir); %Generates S. Figure 3b & S.Table 3, Sheet 2



