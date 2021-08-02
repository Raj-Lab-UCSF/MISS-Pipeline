function rvalNtau_V_fronorm(rval_pv,rval_sst,rval_vip,tauvec,outstruct)

curvenames = {'||Residual||_F','R-Pv+','R-Sst+','R-Vip+','\tau_A_d_j','Elbow'};

fronorm = zeros(length(outstruct),1);
ng_param_list = zeros(length(outstruct),1);
for i = 1:length(outstruct)
    fronorm(i) = outstruct(i).fronorm;
    ng_param_list(i) = outstruct(i).nGen;
end
fronorm = fronorm ./ ng_param_list;
error = fronorm;

norm_error = (error - min(error)) / (max(error) - min(error));

figure('Position',[0 0 1000 400]); hold on;
plot(ng_param_list,norm_error,'k-','LineWidth',2); hold on;
set(gca,'FontSize',18)
yticks([0 0.2 0.4 0.6 0.8 1]);
xticks([0 0.25*3855 0.5*3855 0.75*3855 3855]);
xticklabels({'0',num2str(floor(3855*0.25)),num2str(floor(3855*0.5)),num2str(floor(3855*0.75)),'3855'});
xlabel('\itnG','FontSize',20);
ylabel('Residual','FontSize',20);
title('Resdiual Vs. MISS QC','FontSize',20);

yyaxis right
yticks([0 0.2 0.4 0.6 0.8 1]);
ylabel('R & \tau Values','FontSize',18);
plot(ng_param_list,rval_pv(:,1),'r-','LineWidth',2); hold on;
plot(ng_param_list,rval_sst(:,1),'g-','LineWidth',2); hold on;
plot(ng_param_list,rval_vip(:,1),'b-','LineWidth',2); hold on;
plot(ng_param_list,tauvec,'y-','LineWidth',2); hold on;
plot([606 606],[0 1],'m--','LineWidth',1.5); hold on;

legend(curvenames);

end