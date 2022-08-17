function elbowind = elbow_selector_human(outstruct,makefig)

%This function uses the outstruct output from Cell_Density_Outstruct to
%create an albow curve and choose an elbow index. This elbow index is then
%used to select what nG value is used to create the MISS cell type maps.
%
%   INPUT
%       outstruct: this is the output structure of Cell_Density_Outstruct
%       and includes all the cell type maps for all nG, performed in the
%       MRx3 listed order of genes. This structure contains a fronorm field
%       which is the Frobenius norm of E - C*D. This field is turned into a
%       vector and is used to generate the elbow curve and select the elbow
%       index value.
%       makefig: a binary flag for indicating whether an elbow curve figure
%       should be generated.
%
%   OUTPUT
%       elbowind: this is the elbow index for selecting what nG value, and
%       therefore what MRx3 geneset, will be used to make MISS cell type
%       maps.
%       Figure: a figure of the elbow curve may also be generated using
%       this function is makefig == 1.

fronorm = zeros(length(outstruct),1);
ng_param_list = fronorm;
for i = 1:length(outstruct)
    fronorm(i) = outstruct(i).fronorm;
    ng_param_list(i) = outstruct(i).nGen;
end

% if ng_param_list(end) == 3855 || ng_param_list(end) == 3803
%     ngmax = ng_param_list(end);
% else
%     if size(outstruct(end).corrB,2) == 25
%         ngmax = 3855;
%     else
%         ngmax = 3803;
%     end
% end
ngmax = 14766;

fronorm = fronorm ./ ng_param_list;
error = fronorm;
norm_error = (error - min(error)) / (max(error) - min(error));
normnG = (ng_param_list / (ngmax));
dist2origin = sqrt((normnG-0).^2 + (norm_error-0).^2);
[~,elbowind] = min(dist2origin);

if makefig
    figure;
    plot(normnG,norm_error,'LineWidth',2.5); hold on;
    plot([normnG(elbowind) normnG(elbowind)],[0 1],'LineWidth',2,'LineStyle','--'); hold on;
    yticks([]);
    xlim([0 1]);
    xticks([0 1]);
    xticklabels({num2str(1),'{\it N}_G'});
    set(gca,'FontSize',18);
end

end
    
    