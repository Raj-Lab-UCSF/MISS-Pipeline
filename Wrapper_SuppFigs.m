%Supplemental Figures Wrapper
matdir = '/Users/christophermezias/Documents/MISS-MatFiles';

%% S.Figure 1
%S.Figure 1a: This panel was largely hand constructed as an illustration
%and cannot be readily reconstructed with code

%S.Figure 1b & 1c: (b) depicts the elbow curve used to choose the nG
%parameter and (c) shows how correlated nG values surrounding the elbow are
%with the chosen elbow.

load([matdir filesep 'inversion_aroundelbow_MRx3Prefilter_l90'],...
    'outstruct');
load([matdir filesep 'majtypes_elbows_dist2origin_prefilter_l90'],...
    'fitstruct','ng_param_list');
fitstruct = fitstruct.l90;
load([matdir filesep 'Tasic_Inputs'],'classkey','genevct');
naround = 100;

SFigure1bc_Generator(fitstruct,outstruct,ng_param_list,naround);

%S.Figure 1d & 1e: These demonstrate the residual between E - C*B in sample
%region voxels versus unsampled region voxels, demonstrating no bias in
%residual between sampled and unsampled regions.

load([matdir filesep 'MRx3_L90_inds.mat'],'geneinds');
nG = outstruct.nGen;
preloadinds = geneinds;
lambda = 90;
savenclose = 0;
S_Figure_Residuals(nG,lambda,preloadinds,savenclose,matdir)

%S.Figure 1f: This is run in its own separate script as it requires a core
%or server or powerful desktop computer to run. DO NOT ATTEMPT TO RUN ON A
%LAPTOP!

%S.Figure 1g: Total and n-Unique transcripts per major cell class in
%scRNAseq data from Tasic, et al., 2018

tx_per_CT(genevct,classkey);

