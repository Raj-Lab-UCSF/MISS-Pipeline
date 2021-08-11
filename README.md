# MISS-Pipeline
[Chris Mezias](https://github.com/chm2062) & [Justin Torok](https://github.com/justin-torok)

The set of code necessary to run the MISS (Matrix Inversion with Subset Selection) pipeline to extract cell counts from scRNAseq and ISH transcriptomic data. Full methodological details can be found in our [preprint](https://www.biorxiv.org/content/10.1101/833566v1).  

## 1. Setup
All code is written in MATLAB and requires version 2018b or later.

Step 1: Clone this repository into a local directory of your choice.<br>
Step 2: Download the following [Zip File](https://drive.google.com/file/d/1nAS0Oqlrm2PYGermmpc0Fx2Umx0AtcYt/view?usp=sharing) and unpack its contents in the local copy of the repository, which contains all of the data dependencies necessary to run the code.
Step 3: Some figures require the use of a 3D brain rendering visualization package, also by the authors, which can be found at the following [GitHub link](https://github.com/Raj-Lab-UCSF/Brainframe.git). Clone this repository and add it to your Matlab startup path.
Step 4: To recreate input data from raw datasets, using the Wrapper_InputDataGen pipeline, please download .csv and .loom raw data files for scRNAseq data from [Tasic, et al., 2018]() and [Zeisel, et al., 2018]() at the links given. Input data cannot be regenerated without this raw data, which should come from source authors, but be warned these files are quite large (5-10GB each). 
Step 5: To recreate cell type maps and analyses and figures from the manuscript, use the Wrapper_MISSManuscript pipeline script. To recreate only the cell type maps, use Wrapper_CellDensityMaps script. To recreate supplementary analyses and datasets, use the Wrapper_SuppMaterial pipeline script.

## 2. Code Files
Below is a short description of each of the code files contained in the MISS-Pipeline folder, grouped by general functionality. Pipeline wrappers are described in their own section first. Functions, grouped by functionality, have their inputs and outputs described, with required inputs in boldface text and optional inputs with their default setting in parentheses. The sections code files are grouped into, in order, are:
- Pipeline Wrappers
- Input Data Generation
- Cell Type Density Map Generation
- Main Text Analyses & Figure Generation
- Supplementary Analyses & Figure Generation

### Pipeline Wrappers
- 'Wrapper_InputDataGen.m': This wrapper recreates the cleaned, formatted, processed, and normalized input data necessary for running the MISS pipeline to create cell type maps. This wrapper takes raw scRNAseq and ISH data, as well as associated metadata, and puts it into the formats necessary for running the MISS pipeline. Users must set "matdir" to the file path the .mat data files are stored in on their local machine.
- 'Wrapper_MISSManuscript.m': This wrapper recreates cell type map outputs for both scRNAseq data from Tasic, et al., 2018 (lines 10-42) and Zeisel, et al., 2018 (lines 105-137) and also recreates the main text analyses, results, and figures, in the order they appear in the MISS manuscript, found at the BioRxiv link above. In multiple places, users have the option to change a binary flag variable called "makenew." Setting this to 0 loads in the output data used in the MISS manuscript, and setting this to 1 causes the MRx3 gene reordering and cell type density maps to be regenerated. Users must set "matdir" to the file path the .mat data files are stored in on their local machine.
- 'Wrapper_SuppMaterial.m': This wrapper recreates the supplementary figures and analyses found in the MISS manuscript, split into sections using the scRNAseq data from Tasic, et al., 2018 and Zeisel, et al., 2018. Within these sections, code is ordered based on the supplementary figure, table, or dataset order from the MISS manuscript. Users must set "matdir" to the file path the .mat data files are stored in on their local machine.

### Input Data Generation
- 'scRNAseq_Data_Extract_Tasic.m': Data preprocessing function that output a structure of per cell gene expression, class ID, and correlation with cell class centroid, organized hierarchically by cell class and subclass divisions. Specific to scRNAseq data from Tasic, et al., 2018, which can be downloaded at the original source following the links above.
	- ***Inputs***:
		- **directory**: File path specified by the user directing data file loads and save calls to the MISS-MatFiles folder.
	- ***Outputs***: 
		- **classstruct**: A structure of per cell gene expression, class ID, and correlation with cell class centroid, organized hierarchically by cell class and subclass divisions.
- 'Cell_Type_Data_Extract_Tasic.m': Data preprocessing function that outputs average expression profiles per cell type or class, along with individual cell expression arrays and associated metadata.
	- ***Inputs***:
		- **classstruct**: A structure of per cell gene expression, class ID, and correlation with cell class centroid, organized hierarchically by cell class and subclass divisions. This forms the basis data for the eventual data arrays output.
		- **excl_names**: Cell classes to exclude from output aggregated data and associated metadata due to paucity or instability.
		- **directory**: File path specified by the user directing data file loads and save calls to the MISS-MatFiles folder.
	- ***Outputs***: 
		- **meanexprmat_ct**: An array of aggregate expression scores across genes per cell class/type, with averaging weighted by each member cell's correlation with the cell type or class centroid expression vector across genes.
		- **meanexprmat_subt**: An array of aggregate expression scores across genes per cell subclass/subtype, with averaging weighted by each member cell's correlation with the cell subtype or subclass centroid expression vector across genes.
		- **classkey**: Character array of cell class/type names.
		- **classkey_subt**: Character array of cell subclass/subtype names.
		- **C_indivcells**: Matrix of individual cell gene expression across the unity set of genes described in gene_names.
		- **ct_labvec**: Vector of cell class/type numeric ID per individual cell.
		- **subt_labvec**: Vector of cell subclass/subtype numeric ID per individual cell.
		- **ct_group**: Vector of major cell class/type numeric ID associated with each subclass/subtype numeric ID.
		- **entrez_names**: Character array of gene names in the unity set between scRNAseq data and the coronal AGEA ISH atlas.
- 'ISH_Data_Extract_Tasic.m': Data preprocessing function that outputs ISH expression values per voxel, across the set of unity genes between Tasic, et al., 2018 scRNAseq data and AGEA ISH expression scores from the coronal atlas.
	- ***Inputs***: 
		- **classstruct**: A structure of per cell gene expression, class ID, and correlation with cell class centroid, organized hierarchically by cell class and subclass divisions.
		- **directory**: File path specified by the user directing data file loads and save calls to the MISS-MatFiles folder.
	- ***Outputs***: 
		- **voxvgene**: ISH expression per voxel using the unity gene set between Tasic, et al., 2018 scRNAseq data and ISH expression scores from the coronal AGEA. This is the same set of genes described in gene_names.
- 'Cell_Type_Data_Extract_Zeisel.m': Data preprocessing function that outputs average expression profiles per cell type or class, along with individual cell expression arrays and associated metadata. Specific to scRNAseq data from Zeisel, et al., 2018. This requires download of loom input files from Zeisel, et al., 2018, described in the links above.
	- ***Inputs***: 
		- **directory**: File path specified by the user directing data file loads and save calls to the MISS-MatFiles folder.
	- ***Outputs***: 
		- **genevct**: An array of aggregate expression scores across genes per cell class/type, with averaging weighted by each member cell's correlation with the cell type or class centroid expression vector across genes.
		- **C_indivcells**: Matrix of individual cell gene expression across the unity set of genes described in gene_names.
		- **classkey**: Character array of cell class/type names.
		- **entrez_names**: Character array of gene names in the unity set between scRNAseq data and the coronal AGEA ISH atlas.
		- **ct_labvec**: Vector of cell class/type numeric ID per individual cell.
		- **ct_namevec**: Vector of numeric IDs per cell class/type.
- 'ISH_Data_Extract_Zeisel.m': Data preprocessing function that outputs ISH expression values per voxel, across the set of unity genes between Zeisel, et al., 2018 scRNAseq data and AGEA ISH expression scores from the coronal atlas.
	- ***Inputs***: 
		- **directory**: File path specified by the user directing data file loads and save calls to the MISS-MatFiles folder.
	- ***Outputs***: 
		- **voxvgene**: ISH expression per voxel using the unity gene set between Zeisel, et al., 2018 scRNAseq data and ISH expression scores from the coronal AGEA. This is the same set of genes described in gene_names.
- 'Zeisel_geneset.m': This function extract the unity set of genes between those chosen for mapping in Zeisel et al., 2018 and the coronal AGEA, the cell classes/types with at least one gene in this unity set, and the index of the cell classes/types with at least one representative gene in this set. This data is used in Figure 5 for comparison purposes.
	- ***Inputs***: 
		- **filedir**: File path specified by the user directing data file loads and save calls to the MISS-MatFiles folder.
	- ***Outputs***: 
		- **Zeisel_gene_names**: Gene names in the unity set of genes between those selected by Zeisel, et al., 2018 for their mapping and the coronal AGEA.
		- **repcells**: Cell classes/types with at least one gene representing them in the unity set given by Zeisel_gene_names.
		- **repcellinds**: Indices of the cell classes/types in classkey and genevct (columns) represented in repcells.

### Cell Type Density Map Generation
- 'MRx3_Selector_Prefilter': This is our gene reordering function, which ranks genes according to their effective entropy across cell types, divided by their redundancy with genes already in the set, with genes that contribute the most to error in E - C*D being automatically penalized to rank towards the end of the list. This algorithms is initialized at the gene with the lowest entropy that does not strongly contribute to projection error in E - C*D. 
	- ***Inputs***: 
		- **C_raw**: An array of aggregate expression scores across genes per cell class/type, with averaging weighted by each member cell's correlation with the cell type or class centroid expression vector across genes.
		- **E**: ISH expression per voxel using the unity gene set between scRNAseq data and ISH expression scores from the coronal AGEA. This is the same set of genes described in gene_names.
		- **n**: Number of cell types/classes.
		- **lambda**: Value governing E - C*D projection error penalty.
		- **useParallel**: Binary flag determine whether multiple cores should be used for parallel processing.
	- ***Outputs***:
		- **geneinds**: Vector of gene ranks, in order, where each entry is the index of genes in the original arbitrary order.
- 'Cell_Density_Outstruct.m': This function wraps together all of the other functions necessary to produce the per voxel cell type maps, including testing across different nG (number of genes) parameter values, performing the linear inversion or correlation mapping, and calculating the E - C*D residual and its Frobenius norm.
	- ***Inputs***:
		- **genevct_**: An array of aggregate expression scores across genes per cell class/type, with averaging weighted by each member cell's correlation with the cell type or class centroid expression vector across genes.
		- **voxvgene_**: ISH expression per voxel using the unity gene set between scRNAseq data and ISH expression scores from the coronal AGEA. This is the same set of genes described in gene_names.
		- **gene_names_**: Character array of gene names in the unity set between scRNAseq data and the coronal AGEA ISH atlas.
		- **ng_param_list_**: A vector of all nG values to be tested.
		- **lambda_**: Value governing E - C*D projection error penalty in the MRx3 algorithm.
		- **missmethod_**: Gene ranking/subsetting algorithm to choose. Default is MRx3, but classic mRMR, colAMD, DBSCAN, DiffExp (Differential Expression) can also be chosen.
		- **infmethod_**: Mapping method using nonnegative linear inversion (inversion), correlation based mapping (corr), or nonnegative linear inversion plus residual and Frobenius norm calculation (inv+res). Default is inv+res.
		- **preloadinds_**: Vector of ranked gene indices produced by the MRx3 (or another) gene ranking algorithm. If this ranking has nor been performed already, preloadinds can be set to an empty array, which is the default.
	- ***Outputs***: 
		- **outstruct**: An output structure for cell type mapping results, including the following fields across indices given by each nG parameter tested: corrB (the per voxel per cell type density matrix), nGen (the nG parameter associated with each corrB map produced), lambda (the lambda value used in the MRx3 algorithm), mnresnorm (mean normalized residual per nG parameter for E - C*D), fronorm (Frobenius norm of the residual for E - C*D), Bsums (the summed value across voxels per region for each cell type), and Bmeans (the average value across voxels per region for each cell type). This data structure object will be the basis for most subsequent analyses. 
- 'GeneSelector.m': This function takes the gene subset, based on the gene ranking algorithm and current nG parameter, and outputs the gene index reduced E and C matrices.
	- ***Inputs***: 
		- **genevct**: An array of aggregate expression scores across genes per cell class/type, with averaging weighted by each member cell's correlation with the cell type or class centroid expression vector across genes.
		- **voxvgene**: ISH expression per voxel using the unity gene set between scRNAseq data and ISH expression scores from the coronal AGEA. This is the same set of genes described in gene_names.
		- **gene_names_**: Character array of gene names in the unity set between scRNAseq data and the coronal AGEA ISH atlas.
		- **ngen_param**: Current nG (number of genes) parameter value.
		- **lambda**: Value governing E - C*D projection error penalty in the MRx3 algorithm.
		- **method**: Gene ranking/subsetting algorithm to choose. Default is MRx3, but classic mRMR, colAMD, DBSCAN, DiffExp (Differential Expression) can also be chosen.
		- **preloadinds_**: Vector of ranked gene indices produced by the MRx3 (or another) gene ranking algorithm. If this ranking has nor been performed already, preloadinds can be set to an empty array, which is the default.
	- ***Outputs***:
		- **E_red**: Matrix of ISH gene expression per gene per voxel using the reduced set of genes based on gene ranking from MRx3 (or another algorithm) and the current nG parameter value.
		- **C_red**: Matrix of RNAseq gene expression across cell types using the reduced set of genes based on gene ranking from MRx3 (or another algorithm) and the current nG parameter value.
		- **nGen**: Current nG parameter value.
		- **reduced_gene_names**: Character array of gene names in the reduced set.
		- **mrmrinds**: Vector of indices of the genes in the reduced set from their original arbitrary order, with current ordering based on MRx3 (or another algorithm's) ranking.
- 'CellDensityInference.m': Function for performing linear inversion C*D = E to create MISS maps, using E_red and C_red.
	- ***Inputs***: 
		- **E_red**: Matrix of ISH gene expression per gene per voxel using the reduced set of genes based on gene ranking from MRx3 (or another algorithm) and the current nG parameter value.
		- **C_red**: Matrix of RNAseq gene expression across cell types using the reduced set of genes based on gene ranking from MRx3 (or another algorithm) and the current nG parameter value.
	- ***Outputs***: 
		- **B**: Matrix of voxels versus cell types where each entry is each cell type's arbitrarily scaled density in that given voxel.
- 'CellDensityInference_resnorm.m': Same as CellDensityInference.m, except that it also calculates the mean normalized and Frobenius norm of the residual.
	- ***Inputs***: 
		- **E_red**: Matrix of ISH gene expression per gene per voxel using the reduced set of genes based on gene ranking from MRx3 (or another algorithm) and the current nG parameter value.
		- **C_red**: Matrix of RNAseq gene expression across cell types using the reduced set of genes based on gene ranking from MRx3 (or another algorithm) and the current nG parameter value.
	- ***Outputs***: 
		- **B**: Matrix of voxels versus cell types where each entry is each cell type's arbitrarily scaled density in that given voxel.
		- **mnresnorm**: Mean normalized residual of E_red - C_red * B.
		- **fronorm**: Frobenius norm of the residual of E_red - C_red * B.
- 'CellDensity_Corr.m': Function for performing correlation mapping to create maps, using E_red and C_red.
	- ***Inputs***: 
		- **E_red**: Matrix of ISH gene expression per gene per voxel using the reduced set of genes based on gene ranking from MRx3 (or another algorithm) and the current nG parameter value.
		- **C_red**: Matrix of RNAseq gene expression across cell types using the reduced set of genes based on gene ranking from MRx3 (or another algorithm) and the current nG parameter value.
	- ***Outputs***: 
		- **B**: Matrix of voxels versus cell types where each entry is each cell type's Pearson correlation (R) value in that given voxel.
- 'Voxel_To_Region.m': This function takes per voxel cell density or correlation values and aggregates them into per region values, with regions given according to the AGEA and mouse CCF. 
	- ***Inputs***: 
		- **D**: Matrix of voxels versus cell types where each entry is each cell type's arbitrarily scaled density or Pearson correlation (R) value in that given voxel.
		- **directory**: File path specified by the user directing data file loads and save calls to the MISS-MatFiles folder.
	- ***Outputs***: 
		- **cell_types_sum**: The summed value across voxels per region for each cell type.
		- **cell_types_mean**: The average value across voxels per region for each cell type.
- 'elbow_selector.m': Function that uses the Frobenius norm of the residual across tested nG parameter values to give an elbow in the resulting curve, where the elbow index selected is the nG value corresponding to the x coordinate of the point in the curve closest to the origin [0,0]. 
	- ***Inputs***: 
		- **outstruct**: An output structure with fields containing all MISS maps at per voxel and per region resolutions, as well as attendant metadata, such as the mean normalized residual, the Frobenius norm of the residual, with the structure having indices across all tested nG parameter values. See above for a more detailed description.
		- **makefig**: Binary flag for indicating whether the elbow curve figure should be generated (1) or not (0).
	- ***Outputs***:
		- **elbowind**: Index of the chosen elbow value in the ng_param_list vector of all tested nG parameter values. To get the nG parameter value associated with the elbow index, simply do either ng_param_list(elbowind) or outstruct(elbowind).nGen.

### Main Text Analyses & Figure Generation
- 'MISS_Brainframe.m': Function that does 3D brain renderings for the MISS manuscript. This calls 'brainframe.m', the 3D brain rendering function in the 'Brainframe' package [linked here](https://github.com/Raj-Lab-UCSF/Brainframe.git). Download the linked package and add it to your Matlab path to recreate these renderings. This is used to produce many figure panels in the manuscript.
	***Inputs***:
		- **outstruct**: An output structure with fields containing all MISS maps at per voxel and per region resolutions, as well as attendant metadata, such as the mean normalized residual, the Frobenius norm of the residual, with the structure having indices across all tested nG parameter values. See above for a more detailed description.
		- **types**: Indices of cell types you want to visualize, i.e. column index of outstruct(i).corrB.
		- **torz**: 'Tasic' or 'Zeisel', whichever scRNAseq study is to be used.
		- **elbowind**: Index of the chosen elbow value in the ng_param_list vector of all tested nG parameter values. To get the nG parameter value associated with the elbow index, simply do either ng_param_list(elbowind) or outstruct(elbowind).nGen. Default is 582.
		- **xfac**: Multiplier on the number of points visualized per voxel relative to input data scale. Default is 1.
		- **savenclose**: logical flag that, when true, saves axial, coronal, sagittal, and/or custom views as low-compression .tiff files and then closes the MATLAB figure. Default is 0.
		- **voxthresh**: Percentile of voxels thresholded as above desired minimum signal for visualization. Written in decimal. Default is 1.
		- **cmap_range**: The colormap for each cell type indicated in types. Should be a cell array of 2x3 vectors n types long. Default is for 1 type and is [0 0 0; 1 1 1].
		- **img_name**: Unique filename string. Default is 'invNsubs'.
		- **matdir**: File path specified by the user directing data file loads and save calls to the MISS-MatFiles folder.
	- ***Outputs***: Outputs are the 3D brain renderings, no variables.
- 'colorbar_creator.m': A function to create colorers based on two color poles per map needing a color scale.
	- ***Inputs***: 
		- **cmap_range**: The colormap for each cell type indicated in types. Should be a cell array of 2x3 vectors n types long.
		- **savenames**: Unique filename strings saved into a cell array.
		- **savenclose**: logical flag that, when true, saves axial, coronal, sagittal, and/or custom views as low-compression .tiff files and then closes the MATLAB figure. Default is 0.
	- ***Outputs***: A colorbar or colorers between the color poles specified.
- 'Kim_Study_Comparison.m': A function that produces correlation plots, with R & $\rho values, and stars indicating significance level included. This is used to produce panels 2b-e in the manuscript.
	- ***Inputs***: 
		- **outstruct**: An output structure with fields containing all MISS maps at per voxel and per region resolutions, as well as attendant metadata, such as the mean normalized residual, the Frobenius norm of the residual, with the structure having indices across all tested nG parameter values. See above for a more detailed description.
		- **idx**: index in outstruct to use. We use elbowind.
		- **region**: Area of brain over which to make comparisons. Options are 'neo' (neocortex), 'fore' (forebrain outside of neocortex), 'neo+fore' (whole forebrain), and 'whole' (whole brain). Default is 'neo'.
		- **savenclose**: logical flag that, when true, saves axial, coronal, sagittal, and/or custom views as low-compression .tiff files and then closes the MATLAB figure. Default is 0.
		- **directory**: File path specified by the user directing data file loads and save calls to the MISS-MatFiles folder.
	- ***Outputs***: Correlation plots, with R & $\rho values, and stars indicating significance level included. No variables produced.
- 'Figure3_taulayerslice.m': This function generates the slice maps of neocortical laminar glutamatergic cells sequenced in Tasic, et al., 2018, visualized at specified slice indices, using the cell type map data in outstruct(i).corrB at the index specified by i. T$ values are included. This is used to produce panel 3a in the manuscript.
	- ***Inputs***: 
		- **outstruct**: An output structure with fields containing all MISS maps at per voxel and per region resolutions, as well as attendant metadata, such as the mean normalized residual, the Frobenius norm of the residual, with the structure having indices across all tested nG parameter values. See above for a more detailed description.
		- **idx**: index in outstruct to use. We use elbowind.
		- **mapmethod**: String name of mapping method used to produce cell type maps.
		- **slicelocs**: slice indices to visualize in AGEA (Allen Gene Expression Atlas) z-axis space.
		- **savenclose**: logical flag that, when true, saves axial, coronal, sagittal, and/or custom views as low-compression .tiff files and then closes the MATLAB figure. Default is 0.
		- **directory**: File path specified by the user directing data file loads and save calls to the MISS-MatFiles folder.
	- ***Outputs***: slice maps of neocortical laminar glutamatergic cells sequenced in Tasic, et al., 2018, visualized at specified slice indices, using the cell type map data in outstruct(i).corrB at the index specified by i. No variables are produced.
- 'rvalNtau_calc.m': This function calculates the r-values between our Pv+, Sst+, and Vip+ interneuron cell type maps and the Kim, et al., 2017 data and the T$ values for the laminar glutamatergic cells. All maps are produced using major class types from Tasic, et al., 2018.
	- ***Inputs***: 
		- **outstruct**: An output structure with fields containing all MISS maps at per voxel and per region resolutions, as well as attendant metadata, such as the mean normalized residual, the Frobenius norm of the residual, with the structure having indices across all tested nG parameter values. See above for a more detailed description.
		- **matdir**: File path specified by the user directing data file loads and save calls to the MISS-MatFiles folder.
	- ***Outputs***:
		- **Rval_pv**: R-values between our Pv+ interneuron maps and the data from Kim, et al., 2017. This is a length(outstruct) X 3 matrix, where the first column is only considering neocortex, the second only forebrain, and the third whole brain. 
		- **Rval_sst**: R-values between our Sst+ interneuron maps and the data from Kim, et al., 2017. This is a length(outstruct) X 3 matrix, where the first column is only considering neocortex, the second only forebrain, and the third whole brain. 
		- **Rval_vip**: R-values between our Vip+ interneuron maps and the data from Kim, et al., 2017. This is a length(outstruct) X 3 matrix, where the first column is only considering neocortex, the second only forebrain, and the third whole brain. 
		- **tauvec**: T$ values for correlating the empirical ranking of the neocortical laminar glutamatergic cells in order from the cortical surface with the one produced using our cell type maps, calculated per slice and then averaged. This is a length(outstruct) X 1 vector.
- 'rvalNtau_V_fronorm.m': This function generates a plot of the Frobenius norm of the residual and the R values and T$ values generated in the 'rvalNtau_calc.m' function above. These values are plotted versus nG, the number of genes used per MRx3 ranked subset. This generates panel 3b in the manuscript.
	- ***Inputs***: 
		- **Rval_pv**: R-values between our Pv+ interneuron maps and the data from Kim, et al., 2017. This is a length(outstruct) X 3 matrix, where the first column is only considering neocortex, the second only forebrain, and the third whole brain. 
		- **Rval_sst**: R-values between our Sst+ interneuron maps and the data from Kim, et al., 2017. This is a length(outstruct) X 3 matrix, where the first column is only considering neocortex, the second only forebrain, and the third whole brain. 
		- **Rval_vip**: R-values between our Vip+ interneuron maps and the data from Kim, et al., 2017. This is a length(outstruct) X 3 matrix, where the first column is only considering neocortex, the second only forebrain, and the third whole brain. 
		- **tauvec**: T$ values for correlating the empirical ranking of the neocortical laminar glutamatergic cells in order from the cortical surface with the one produced using our cell type maps, calculated per slice and then averaged. This is a length(outstruct) X 1 vector.
		- **outstruct**: An output structure with fields containing all MISS maps at per voxel and per region resolutions, as well as attendant metadata, such as the mean normalized residual, the Frobenius norm of the residual, with the structure having indices across all tested nG parameter values. See above for a more detailed description.
	- ***Outputs***: a plot of the Frobenius norm of the residual and the R values and T$ values generated in the 'rvalNtau_calc.m' function above. These values are plotted versus nG, the number of genes used per MRx3 ranked subset. No variables are output.
- 'Figure_4d_glia.m': This function generates a bar plot of cell density per major region groupings of interest across glial cell types in our maps, using data from Tasic, et al., 2018. This generates panel 4d in the manuscript.
	- ***Inputs***: 
		- **outstruct**: An output structure with fields containing all MISS maps at per voxel and per region resolutions, as well as attendant metadata, such as the mean normalized residual, the Frobenius norm of the residual, with the structure having indices across all tested nG parameter values. See above for a more detailed description.
		- **idx**: index in outstruct to use. We use elbowind.
		- **savenclose**: logical flag that, when true, saves axial, coronal, sagittal, and/or custom views as low-compression .tiff files and then closes the MATLAB figure. Default is 0.
		- **directory**: File path specified by the user directing data file loads and save calls to the MISS-MatFiles folder.
	- ***Outputs***: A bar plot of cell density per major region groupings of interest across glial cell types in our maps. No variables are output.
-'Figure_5abcd_zeisel_typemaps.m': This function generates slice maps using the geneset and correlation mapping procedure from Zeisel, et al., 2018, along with our maps produced using MRx3 derived gene subsets on the scRNAseq data from Zeisel, et al., 2018. Scatterplots and correlations between these two mappings are also included. This is used to generate panels 5a-d in the manuscript.
	- ***Inputs***: 
		- **outstruct**: An output structure with fields containing all MISS maps at per voxel and per region resolutions, as well as attendant metadata, such as the mean normalized residual, the Frobenius norm of the residual, with the structure having indices across all tested nG parameter values. See above for a more detailed description.
		- **idx**: index in outstruct to use. We use elbowind.
		- **typeinds**: Indices of cell types you want to visualize, i.e. column index of outstruct(i).corrB.
		- **slicelocs**: slice indices to visualize in AGEA (Allen Gene Expression Atlas) z-axis space.
		- **savenclose**: logical flag that, when true, saves axial, coronal, sagittal, and/or custom views as low-compression .tiff files and then closes the MATLAB figure. Default is 0.
		- **directory**: File path specified by the user directing data file loads and save calls to the MISS-MatFiles folder.
	- ***Outputs***: Slice maps using the geneset and correlation mapping procedure from Zeisel, et al., 2018, along with our maps produced using MRx3 derived gene subsets on the scRNAseq data from Zeisel, et al., 2018. Scatterplots and correlations between these two mappings are also included. No variables are output.
- 'Figure_5e_MISS_correlations_comp.m': This function generates a box and whisker plot, with points included, of the R value between all MISS and Zeisel, et al., 2018 maps across all cell types, grouped into major classes of types. This is used to generate panel 5e in the manuscript.
	- ***Inputs***: 
		- **outstruct**: An output structure with fields containing all MISS maps at per voxel and per region resolutions, as well as attendant metadata, such as the mean normalized residual, the Frobenius norm of the residual, with the structure having indices across all tested nG parameter values. See above for a more detailed description.
		- **idx**: index in outstruct to use. We use elbowind.
		- **slicelocs**: slice indices to visualize in AGEA (Allen Gene Expression Atlas) z-axis space.
		- **savenclose**: logical flag that, when true, saves axial, coronal, sagittal, and/or custom views as low-compression .tiff files and then closes the MATLAB figure. Default is 0.
		- **directory**: File path specified by the user directing data file loads and save calls to the MISS-MatFiles folder.
	- ***Outputs***: a box and whisker plot, with points included, of the R value between all MISS and Zeisel, et al., 2018 maps across all cell types, grouped into major classes of types.

### Supplementary Analyses & Figure Generation
- 'S_Figure_Residuals.m': This function generates a box and whisker plot and 3D brain rendering using spheres of the per-region residuals, at the elbow index cell type maps, in voxels that are from the scRNAseq sampled regions in Tasic, et al., 2018, and those regions that were not sampled. This generates supplementary panels 2a-b in the manuscript.
	- ***Inputs***: 
		- **nG**: Current nG parameter value. Default is 582, the elbow index.
		- **lambda**: Scalar indicating the lambda value for MRx3-based subset selection. Default is 90.
		- **preloadinds**: Vector of gene ranks, in order, where each entry is the index of genes in the original arbitrary order.
		- **savenclose**: logical flag that, when true, saves axial, coronal, sagittal, and/or custom views as low-compression .tiff files and then closes the MATLAB figure. Default is 0.
		- **directory**: File path specified by the user directing data file loads and save calls to the MISS-MatFiles folder.
	- ***Outputs***: A box and whisker plot and 3D brain rendering using spheres of the per-region residuals, at the elbow index cell type maps, in voxels that are from the scRNAseq sampled regions in Tasic, et al., 2018, and those regions that were not sampled.
- 'SFigure2c_Generator.m': This function creates a line plot of the r-values between cell type maps, per cell type, between different nG parameters, within a specified range. This generates supplementary panel 2c in the manuscript.
	- ***Inputs***: 
		- **elbowind**: Index of the chosen elbow value in the ng_param_list vector of all tested nG parameter values. To get the nG parameter value associated with the elbow index, simply do either ng_param_list(elbowind) or outstruct(elbowind).nGen.
		- **outstruct**: An output structure with fields containing all MISS maps at per voxel and per region resolutions, as well as attendant metadata, such as the mean normalized residual, the Frobenius norm of the residual, with the structure having indices across all tested nG parameter values. See above for a more detailed description.
		- **ng_param_list_**: A vector of all nG values mapped.
		- **naround**: The distance around the elbow index in number of genes over which to perform the correlations as denoted in the function description.
	- ***Outputs***: A line plot of the r-values between cell type maps, per cell type, between different nG parameters, within a specified range. No variables are output.
- 'tx_per_CT.m': This function generates bar plots of the number of total and unique transcripts per cell type using the MISS maps generated with the scRNAseq data from Tasic, et al., 2018. This is used to generates supplementary panels 3 c and d.
	- ***Inputs***: 
		- **genevct**: An array of aggregate expression scores across genes per cell class/type, with averaging weighted by each member cell's correlation with the cell type or class centroid expression vector across genes.
		- **classkey**: Character array of cell class/type names.
	- ***Outputs***: Bar plots of the number of total and unique transcripts per cell type using the MISS maps generated with the scRNAseq data from Tasic, et al., 2018
- 'GeneList_Generator.m': This function generates the MRx3 ranked subset list of gene names, in rank order, using the data from the specified scRNAseq dataset and the specified elbow index, and the specified MRx3 gene ranks. This function also generates an expression intensity map of all genes in the MRx3 ranked subset across all cell types. This function is used to generate supplementary panel 3a and supplementary table 3.
	- ***Inputs***:
		- **study**: A string specifying which scRNAseq study, Tasic, et al., 2018 ('Tasic') or Zeisel, et al., 2018 ('Zeisel') is to be used.
		- **elbowind**: Index of the chosen elbow value in the ng_param_list vector of all tested nG parameter values. To get the nG parameter value associated with the elbow index, simply do either ng_param_list(elbowind) or outstruct(elbowind).nGen.
		- **geneinds**: Vector of gene ranks, in order, where each entry is the index of genes in the original arbitrary order.
		- **directory**: File path specified by the user directing data file loads and save calls to the MISS-MatFiles folder.
	- ***Outputs***: An expression intensity map of all genes in the MRx3 ranked subset across all cell types. Variable output is described below.
		- **genelist**: the MRx3 ranked subset list of gene names, in rank order.		
### Supplementary Code
- 'Wrapper_CellDensityMaps.m': This wrapper script only performs the cell type mapping procedures sections from the 'Wrapper_MISSManuscript.m' script, and excludes all analysis and figure plotting code. This can be used to recreate the cell type maps using the two scRNAseq data sets employed in the current pipeline, or it can be modified to create cell type maps using other scRNAseq datasets.
- 'mRMR_Selector.m': This is the mRMR gene ranking algorithm used in MRx3 prefilter, which ranks genes according to their effective entropy across cell types, divided by their redundancy with genes already in the set. Residual error is not considered in this algorithm, differentiating it from MRx3.

## Data Files
- **allnG_MRx3Prefilter_l90.mat**: Outstruct of MISS cell type maps produced using all the genes in common between Tasic, et al., 2018 scRNAseq data and the AGEA ISH atlas.
- **CellDensity_corr_Tasic.mat**: Outstruct of elbow index and all gene correlation maps off all cell types from the Tasic, et al., 2018 scRNAseq data.
- **default_mouse.mat**: A struct object containing default mouse brain mapping parameters for the 3D brain renderings.
- **input_struct_voxelrender.mat**: The metadata and data necessary for creating the slice map visualizations in this manuscript.
- **ISH_gene_names.mat**: Gene names from the AGEA ISH atlas.
- **ISH_input_data.mat**: All numeric input data and metadata relating to the AGEA ISH atlas necessary for producing the MISS cell type maps.
- **kim_density_listB_order.mat**: Pv+, Sst+, and Vip+ interneuron densities per region.
- **kim_totals_reorder_m.mat**: Total Pv+, Sst+, and Vip+ interneuron counts per region.
- **listB.mat**: Brain region names and major region groupings.
- **mouse_ALM_2018-06-14_samples-columns.mat**: Anterolateral motor area sampled scRNAseq data from Tasic, et al., 2018.
- **mouse_LGd_2018-06-14_samples-columns.mat**: Dorsal lateral geniculate complex sampled scRNAseq data from Tasic, et al., 2018.
- **mouse_VISp_2018-06-14_genes-rows.mat**: Gene names and scRNAseq metadata from Basic, et al., 2018.
- **mouse_VISp_2018-06-14_samples-columns.mat**: Primary visual cortex (V1) sampled scRNAseq data from Tasic, et al., 2018.
- **MRx3_L90_inds.mat**: MRx3 ranked genes, listed by original arbitrary ordering index, using the scRNAseq data from Tasic, et al., 2018.
- **ReadData.mat**: 
- **regionlabs.mat**: Numeric region labels applied per voxel in the mouse brain.
- **rval_wholerange_Tasic.mat**: A nG X 3 matrix of r-values between MISS generated cell type maps and cell density data from Kim, et al., 2017, across Pv+, Sst+, and Vip+ interneurons, using the scRNAseq data from Tasic, et al., 2018 and the AGEA to create the cell type maps.
- **Tasic_Inputs.mat**: All of the input data and metadata from Tasic, et al., 2018 necessary to create the cell type maps, in conjunction with the AGEA ISH atlas.
- **Tasic_MRx3Prefilter_fulltau.mat**: An nG x 1 vector of tau values between the empirical and MISS cell type map derived ranks of neocortical laminar glutamatergic neurons by distance from the cortical surface, with MISS maps generated using the scRNAseq data from Tasic, et al., 2018.
- **Tasic_outstruct.mat**: The outstruct MISS cell type map output for all tested nG parameters, using the scRNAseq data from Tasic, et al., 2018.
- **tau_calc_dependencies.mat**: Metadata necessary for generating empirical tau ranks and slice maps before calculating tau values.
- **Zeisel_cellIDs.mat**: Numeric cell type IDs and whether each type meets inclusion or exclusion criteria after QC analysis, all derived from Zeisel, et al., 2018.
- **Zeisel_coronal_geneset.mat**: The unity set between the gene subset from Zeisel, et al., 2018 and the coronal AGEA ISH atlas.
- **Zeisel_Inputs**: All of the input data and metadata from Zeisel, et al., 2018 necessary to create the cell type maps, in conjunction with the AGEA ISH atlas.
- **Zeisel_MRx3Inds.mat**: MRx3 ranked genes, listed by original arbitrary ordering index, using the scRNAseq data from Zeisel, et al., 2018.
- **Zeisel_outstruct.mat**: The outstruct MISS cell type map output for all tested nG parameters, using the scRNAseq data from Zeisel, et al., 2018.