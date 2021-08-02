function Figure_4d_glia(outstruct,idx,savenclose,directory)
% Designed to handle an outstruct of dimension length(ngenelist) that has
% field "corrB." idx is the index of interest that is within the range of
% ngenelist

if nargin < 4
    directory = [cd filesep 'MatFiles'];
    if nargin < 3
        savenclose = 0;
    end
end

load([directory filesep 'Tasic_Inputs.mat'],'listBmap');
listBmap = listBmap(:);
for i = 1:212
    voxels(i) = sum(listBmap==i);
end
reglabs = {'Amg','Cer','Sub','Hip','Hyp','Neo','Med','Mid','Olf','Pal','Pns','Str','Tha';...
            1:11,12:22,23:25,26:36,37:56,57:94,95:119,120:140,141:148,149:156,157:169,170:177,178:212};
astrosum = (sum(outstruct(idx).Bmeans(1:213,22),2) + sum(outstruct(idx).Bmeans(214:end,22),2)) / 2;
microsum = (sum(outstruct(idx).Bmeans(1:213,23),2) + sum(outstruct(idx).Bmeans(214:end,23),2)) / 2;
oligosum = (sum(outstruct(idx).Bmeans(1:213,24),2) + sum(outstruct(idx).Bmeans(214:end,24),2)) / 2;
for i = 1:size(reglabs,2)
    curinds = reglabs{2,i};
    astroplot(i) = mean(astrosum(curinds)/sum(astrosum(:)));
    microplot(i) = mean(microsum(curinds)/sum(microsum(:)));
    oligoplot(i) = mean(oligosum(curinds)/sum(oligosum(:)));
end
gliaplot = [astroplot; microplot; oligoplot]; gliaplot = gliaplot.';

f1 = figure('Units','inch','Position',[0 0 12 4]); hold on;

figure(f1);
bar(gliaplot);
set(gca,'XTick',1:size(reglabs,2));
set(gca,'XTickLabel',reglabs(1,:));
set(gca,'YTick',[0 0.005 0.01]);
ylabel('Normalized Cell Density');
set(gca,'TickLength',[0 0]);
legend('Astro','Micro','Oligo','Location','northeast');
set(gca,'FontSize',18);
title('Mean Normalized Glial Density Per Major Region','FontSize',20);

if savenclose
    print('Figure_3de_gliaplots','-dtiffn');
    close
end
end