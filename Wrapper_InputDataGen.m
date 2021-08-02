%FILEPATH
matdir = '/Users/christophermezias/Documents/MISS-MatFiles';

%TASIC, ET AL., 2018 PREPROCESSING AND PROCESSED DATA GENERATION
classstruct = scRNAseq_Data_Extract_Tasic(matdir); %generating struct object for all cells, hierarchically clustered, with reads appropriately normalized
excl_names = 'VLMC'; %cell classes to exclude due to paucity of data
[genevct, genevct_subt,classkey,classkey_subt, ...
    C_indivcells,ct_labvec,subt_labvec,ct_group,gene_names] = ... 
    Cell_Type_Data_Extract_Tasic(classstruct,excl_names,matdir); %generating processed scRNAseq input data
voxvgene = ISH_Data_Extract_Tasic(classstruct,matdir); %processing and generating ISH AGEA input data
load([matdir filesep 'ISH_input_data.mat'],'GENGDmod',...
    'structList','structIndex','nonzerovox','listBmap'); %loading in ISH input data variables

%SAVING TASIC, ET AL., 2018 INPUT DATA
save([matdir filesep 'Tasic_Inputs.mat'],'genevct','genevct_subt',...
    'classkey','classkey_subt','C_indivcells','ct_labvec','subt_labvec',...
    'ct_group','gene_names','voxvgene','GENGDmod','listBmap',...
    'structList','structIndex','nonzerovox','-v7.3'); %saving Tasic, et al., 2018 input data

%CLEAR!
clearvars -except matdir GENGDmod listBmap structList structIndex nonzerovox

%ZEISEL, ET AL., 2018 PREPROCESSING AND PROCESSED DATA GENERATION
[genevct,C_indivcells,classkey,gene_names,...
    ct_labvec,ct_namevec] = Cell_Type_Data_Extract_Zeisel(matdir); %Zeisel scRNAseq input data processing and generation
voxvgene = ISH_Data_Extract_Zeisel(matdir); %Zeisel ISH
[Zeisel_gene_names,repcells,repcellinds] = Zeisel_geneset(matdir); %Zeisel mapping unity geneset with coronal AGEA ISH atlas
load([matdir filesep 'ISH_input_data.mat'],'GENGDmod',...
    'structList','structIndex','nonzerovox','listBmap'); %loading in ISH input data variables

%SAVING ZEISEL, ET AL., 2018 INPUT DATA
save([matdir filesep 'Zeisel_Inputs.mat'],'genevct','C_indivcells',...
    'classkey','gene_names','ct_labvec','ct_namevec',...
    'voxvgene','GENGDmod','listBmap','structList','structIndex',...
    'nonzerovox','-v7.3'); %saving Zeisel, et al., 2018 input data
save([filedir filesep 'Zeisel_coronal_geneset.mat'],...
    'Zeisel_gene_names','repcells','repcellinds'); %saving Zeisel coronal geneset data for Figure 5