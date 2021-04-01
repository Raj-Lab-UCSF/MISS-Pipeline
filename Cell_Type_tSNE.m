function Cell_Type_tSNE(study, savenclose, directory)

if nargin < 3
    directory = [cd filesep 'MatFiles'];
    if nargin < 2
        savenclose = 0;
        if nargin < 1
            study = 'tasic';
        end
    end
end
rng(0);
if strcmp(study,'zeisel')
    load([directory filesep 'Zeisel_Inputs.mat'],'genevct');
elseif strcmp(study,'tasic')
    load([directory filesep 'Tasic_Inputs.mat'],'genevct');
end

C = genevct;
ctmean = mean(C,1);
ctmean = repmat(ctmean,size(C,1),1);
C = C ./ ctmean;
Y = tsne(C);
[idx] = kmeans(C,size(C,2));
figure;
gscatter(Y(:,1), Y(:,2),idx,[],[],10,0)
% xlabel('tSNE Dimension 1'); ylabel('tSNE Dimension 2');
xticks([]); yticks([]);
set(gca,'XTickLabel',[],'YTickLabel',[]);
set(gca,'FontSize',20);
if savenclose
    print('tSNE_C','tiffn');
    close
end
end