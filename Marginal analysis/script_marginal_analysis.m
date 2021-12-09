

%%% This script details how to estimate a GP regression model using a
%%% matrice of predictors xs (for the scale) and xg (for the shape), and a vector of exceedances y.
%%% Select a subset of the exceedances assumed to have homogeneous extreme
%%% value regression models (e.g. from a single investment style). 

% Set the predictor in matrices called xs and xg. Standardize the
% predictors with sample means and standard deviations
xs=(xs-repmat(mean(xs),length(xs),1))./repmat(std(xs),length(xs),1);
xg=(xg-repmat(mean(xg),length(xg),1))./repmat(std(xg),length(xg),1);
% Define the dimension of the matrices
ds=size(xs,2);
dg=size(xg,2);
% p0: starting values. Replace '2' and '.8' by unconditional estimates of the
% scale (in log scale) and the shape parameter, computed e.g. with the
% gpfit function on y.
% Rem: here we use a unit link function (nonlog=1) for the shape. If an exponential link function is used for the
% shape (nonlog=0), use a log scale for starting values of the shape.
% Remind that we always use an exponential link function for the scale
% parameter since strict positivity is necessary.
nonlog=1;
p0=[2;zeros(ds,1)+0.001;.8;zeros(dg,1)+0.001]; % assume effects close to 0 for all covariates.
% Estimate the model. See gpd_regression for exact definitions of the
% output. sigma_estimated and gamma_estimated give the estimated shape and
% scale parameter for each observation. hessian can be used for inference.
[params_estimated,~,~,~,grad0,hessian0,sigma_estimated,gamma_estimated]=gpd_regression(p0,y,xs,xg,nonlog);

% Repeat with a regularization procedure, for a grid of regularization
% (\lambda) parametr, as in Hambuckers et al. (2018), Journal of Applied
% Econometrics.
% Build the grid of regularization parameters
ngrid=40; % number of points. Overall, we need to perform ngrid^2 estimations
% Bounds of the grid
gr=[0.001 0.3];
% HS: grid of parameter for the scale. HG: same for the shape
HS=linspace(gr(1),gr(2),ngrid); % LASSO grid
HG=linspace(gr(1),gr(2),ngrid);

% Estimate the model via penalized likelihood for each combination of the
% grid.
l=1; % must be one if we use the LASSO (L1 norm)
for j=1:ngrid
    for k=1:ngrid
        [params_estimated_penalized{j}(:,k),sbic2(k,j)]=gpd_regression_auto(p0,y,xs,xg,[HS(k) HG(j)],l,[],[],1);
        % We store the estimated parameter, and the BIC criterion (computed
        % with the active set as approximate degrees of freedom)
    end
end

% Select the penalization parameter with the smallest BIC
[~,idx]=min(sbic2(:));
[ids,idg]=ind2sub(size(sbic2),idx);
% Selected parameters ?
nu_g=HG(idg);
nu_s=HS(ids);
% Isolate the biased penalized likelihood estimator 
params_fin=params_estimated_penalized{idg}(:,ids);
% Compute the post-LASSO estimator
idsfin=find(pfin(2:ds+1)~=0);
idgfin=find(pfin(ds+3:end)~=0);
% Starting values for estimating the post-LASSO estimator. Should be
% updated accordingly for the constants.
p00=[2;zeros(length(idsfin),1)+0.001;.7;zeros(length(idgfin),1)+0.005];
[params_ult,~,~,~,~,hessian,sigma_ultimate,gamma_ultimate]=gpd_regression(p00,y,xs(:,idsfin),xg(:,idgfin),1);

%Graphical goodness-of-fit test based on quantile residuals
% PIT quantile residuals (should be gaussian)
rr=norminv(gpcdf(y,gamma_ultimate,sigma_ultimate,0));
% Theoretical quantiles
qqn=norminv(0.001:.001:0.999);
% Empirical quantiles
qqe=quantile(rr,0.001:.001:.999);
% QQ-plot
close all
figure
plot(qqn,qqn)
hold on
scatter(qqn,qqe,'.r')
title('my strategy')