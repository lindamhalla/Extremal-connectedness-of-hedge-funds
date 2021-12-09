function [LL,sigma,gamma] = gpdregress_app_full2_auto(params,y,xs,xg,hs,l,id_unpenS,id_unpenG,nonlog)

%%% GPD regression, with both parameters that are functions of covariates, and l-norm penalty for automatic variable selection 

% params: column vector of size (d(sigma)+d(gamma)+2)x1. The parameters of
% sigma are given first, then the parameters of gamma (the cst is the first
% parameter).
% y: vector of observed response variable (n X 1)
% xs: matrix of explanatory variable for sigma (n x d(sigma))
% xg: matrix of explanatory variables for gamma (n x d(gamma))
% hs: vector of smoothing parameter (2 x 1)
% l: power of the norm penalty (l=1 is the LASSO case)
% b0s: estimated parameters for sigma, at the previous iteration, used for quadratic
% approximation. Vector of size d(sigma)x1
% b0g: estimated parameters for gamma, at the previous iteration, used for quadratic
% approximation. Vector of size d(gamma)x1
% id_unpenS: position in the vector of parameters, of the ones that are
% not penalized, for sigma.
% id_unpenG: position in the vector of parameters, of the ones that are
% not  penalized, for gamma.

% Parametric gamlss model for GPD distribution with scale sigma and shape gamma - probability function

%%% Define the sizes
ds=size(xs,2);


%%% Reshape the parameters
cs=params(1);
betas=params(2:(ds+1));
if isempty(betas)
    sigma=exp(cs);
else
    sigma=exp(cs+xs*betas); % strict positivity
end
cg=params(ds+2);
betag=params(ds+3:end);
if nonlog==0
    if isempty(betag)
        gamma=exp(cg);
    else
        gamma=exp(cg+xg*betag); % strict positivity
    end
else
    if isempty(betag)
        gamma=(cg);
    else
        gamma=(cg+xg*betag); % strict positivity
    end
end

%%% Penalties

% penalty for sigma

penS=sqrt(betas.^2+10^-7);
%0.5*l(1)*(b0s.^2+10^-7).^(l(1)/2-1).*(betas.^2); 
% dpenfinS=length(y)*hs(1)*l(1)*(b0s.^2+10^-7).^(l(1)/2-1).*betas;
% d2penfinS=length(y)*hs(1)*l(1)*(b0s.^2+10^-7).^(l(1)/2-1);

% Zero penalty for un-regularized variables
penS(id_unpenS)=0;
% dpenfinS(id_unpenS)=0;
% d2penfinS(id_unpenS)=0;
% penalty for gamma

penG=sqrt(betag.^2+10^-7);
%5=0.5*l(2)*(b0g.^2+10^-7).^(l(2)/2-1).*(betag.^2); 
% dpenfinG=length(y)*hs(2)*l(2)*(b0g.^2+10^-7).^(l(2)/2-1).*betag;
% d2penfinG=length(y)*hs(2)*l(2)*(b0g.^2+10^-7).^(l(2)/2-1);

% Zero penalty for un-regularized variables
penG(id_unpenG)=0;
% dpenfinG(id_unpenG)=0;
% d2penfinG(id_unpenG)=0;

%%% Negative log-likelihood

lnopen = sum(log(gppdf(y,gamma,sigma,0)));
LL=-lnopen+length(y)*hs(2)*sum(penG)+length(y)*hs(1)*sum(penS);


%Gradient of the negative LL% NOT VALID FOR nonlog=1. Use numerical
%approximation instead.

% dlds=-sum(repmat(((1./gamma+1).*y.*gamma)./(sigma.*(1+y.*gamma./sigma))-1,1,ds+1).*[ones(n,1) xs]);
% if nonlog==0
%     dldg=-sum(repmat(log(1+gamma.*y./sigma)./(gamma.^2)-(1./gamma+1).*(y./sigma)./(1+gamma.*y./sigma),1,dg+1).*repmat(gamma,1,dg+1).*[ones(n,1) xg]);
% else
%     dd=repmat(log(1+gamma.*y./sigma)./(gamma.^2)-(1./gamma+1).*(y./sigma)./(1+gamma.*y./sigma),1,dg+1).*[ones(n,1) xg];
% %     dd(y>-sigma./gamma,:)=100;
%     dldg=-sum(dd);
% end

%grad=[dlds(:);dldg(:)]+[0;dpenfinS;0;dpenfinG];

% Hessian

% H=hessian_gpd_reg(y,sigma,gamma,xg,xs,ds,dg,n,nonlog);
% penH=[0;d2penfinS;0;d2penfinG];
% H=H+diag(penH);


end

