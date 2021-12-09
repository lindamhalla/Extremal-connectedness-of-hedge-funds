function [LL,grad,H,sigma,gamma] = gpdregress_full2(params,y,xs,xg,nonlog)

%%% GPD regression, with both parameters that are functions of covariates. 

% params: column vector of size (d(sigma)+d(gamma)+2)x1. The parameters of
% sigma are given first, then the parameters of gamma (the cst is the first
% parameter).
% y: vector of observed response variable (n X 1)
% xs: matrix of explanatory variable for sigma (n x d(sigma))
% xg: matrix of explanatory variables for gamma (n x d(gamma))


% Parametric gamlss model for GPD distribution with scale sigma and shape gamma - probability function

%%% Define the sizes
ds=size(xs,2);
n=size(y,1);
dg=size(xg,2);

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

%%% Negative log-likelihood

LL = -sum(log(gppdf(y,gamma,sigma,0)));


%Gradient of the negative LL% 

dlds=-sum(repmat(((1./gamma+1).*y.*gamma)./(sigma.*(1+y.*gamma./sigma))-1,1,ds+1).*[ones(n,1) xs]);
if nonlog==0
    dldg=-sum(repmat(log(1+gamma.*y./sigma)./(gamma.^2)-(1./gamma+1).*(y./sigma)./(1+gamma.*y./sigma),1,dg+1).*repmat(gamma,1,dg+1).*[ones(n,1) xg]);
else
    dldg=-sum(repmat(log(1+gamma.*y./sigma)./(gamma.^2)-(1./gamma+1).*(y./sigma)./(1+gamma.*y./sigma),1,dg+1).*[ones(n,1) xg]);
end

grad=[dlds(:);dldg(:)];

% Hessian % Analytical Hessian gives similar values than the ones obtained numerically.

% hessian for sigma parameters

ccs=(-sigma.*y.*gamma).*(1./gamma+1)./(sigma+y.*gamma).^2;
XS=[ones(n,1) xs];
% Hs=nan(ds+1);
% for i=1:ds+1
%     for j=1:ds+1
%         Hs(i,j)=sum(XS(:,i).*XS(:,j).*ccs);
%     end
% end
Hs=XS'*(XS.*repmat(ccs,1,ds+1));
Hs=-Hs;

% hessian for gamma parameters
if nonlog==0
    ccg=(y./sigma)./(1+gamma.*y./sigma)-log(1+gamma.*y./sigma)./gamma-(y./sigma).*((1+gamma.*y./sigma).*gamma-(1+gamma).*(y./sigma).*gamma)./(1+gamma.*y./sigma).^2;
    XG=[ones(n,1) xg];
    % Hg=nan(dg+1);
    % for i=1:dg+1
    %     for j=1:dg+1
    %         Hg(i,j)=sum(XG(:,i).*XG(:,j).*ccg);
    %     end
    % end
    Hg=XG'*(XG.*repmat(ccg,1,dg+1));
    Hg=-Hg;
else
    ccg=y./(1+gamma.*y./sigma)./(sigma.*gamma.^2)-2*log(1+gamma.*y./sigma)./(gamma.^3)+((y./sigma).^2).*(1./gamma+1)./((1+gamma.*y./sigma).^2)+(y./sigma)./(1+gamma.*y./sigma)./(gamma.^2);
    XG=[ones(n,1) xg];
    % Hg=nan(dg+1);
    % for i=1:dg+1
    %     for j=1:dg+1
    %         Hg(i,j)=sum(XG(:,i).*XG(:,j).*ccg);
    %     end
    % end
    Hg=XG'*(XG.*repmat(ccg,1,dg+1));
    Hg=-Hg;
end

%%% cross-derivatives (between sigma and gamma parameters)
if nonlog==0
    ccsg=y.*gamma.*(sigma-y)./((sigma+y.*gamma).^2);
    % Cd=nan(ds+1,dg+1);
    % for i=1:ds+1
    %     for j=1:dg+1
    %         Cd(i,j)=sum(XS(:,i).*XG(:,j).*ccsg);
    %     end
    % end
    Cd=XS'*(XG.*repmat(ccsg,1,dg+1));
    Cd=-Cd;
else
    %ccs=(-sigma.*y.*gamma).*(1./gamma+1)./(sigma+y.*gamma).^2;
    %dlds=-sum(repmat(((1./gamma+1).*y.*gamma)./(sigma.*(1+y.*gamma./sigma))-1,1,ds+1).*[ones(n,1) xs]);
    %ccsg=y./sigma./(1+y.*gamma./sigma)-((y./sigma).^2).*(1+gamma)./((1+gamma.*y./sigma).^2);
    ccsg=y.*(sigma-y)./((sigma+y.*gamma).^2);
    % Cd=nan(ds+1,dg+1);
    % for i=1:ds+1
    %     for j=1:dg+1
    %         Cd(i,j)=sum(XS(:,i).*XG(:,j).*ccsg);
    %     end
    % end
    Cd=XS'*(XG.*repmat(ccsg,1,dg+1));
    Cd=-Cd;
end

% Final hessian

H=zeros(ds+dg+2);
H(1:ds+1,1:ds+1)=Hs;
H(ds+2:end,ds+2:end)=Hg;
H(1:ds+1,ds+2:end)=Cd;
H(ds+2:end,1:ds+1)=Cd';


end

