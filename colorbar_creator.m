function colorbar_creator(cmap_range,savenames)

if nargin < 2
    savenames = {'Astro','Micro','Oligo','Endo'}; %glia filenames
    if nargin < 1
        cmap_range = {[1 1 1;0.5 1 0.62],[1 0.75 0;1 0.25 0],[1 1 1;1 0.5 0.5],[0.75 0.25 0.4;1 0.25 0.1]}; %glia colormaps
    end
end

for i = 1:length(cmap_range)
    color1 = cmap_range{i}(1,:);
    color2 = cmap_range{i}(2,:);
    cmap = twocolor(color1,color2,500);
    figure('Units','inches','Position',[0 0 1.5 10])
    imagesc(fliplr(1:500).'); colormap(cmap);
    xticks([]); yticks([]);
    saveas(gcf,savenames{i},'tiff')
end
