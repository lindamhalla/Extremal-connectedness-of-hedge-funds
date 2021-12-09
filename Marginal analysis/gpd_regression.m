

function [pest,fval,exitflag,output,grad,hessian,sigma,gamma]=gpd_regression(p0,y,xs,xg,nonlog)

% p0: starting value for the estimation procedure.
% y: vecteur of exceedances (n x 1)
% xs: matrix of J explanatory variables for sigma (n x J) (si pas de
% covariable - i.e. only one constant: give an empty vector [])
% xg: matrix of K explanatory variables for gamma (n x K)(si pas de
% covariable - i.e. only one constant: give an empty vector [])
% nonlog: 1 si pas de transformation exponential, 0 sinon (0 si
% transformation exponentiel)

% %Example:
% n=1000;
% xs=randn(n,3);
% xg=randn(n,3);
% betag=[.5 0.1 0.1 0.1];
% betas=[.1 0.1 0.1 0.1];
% gamma=exp([ones(n,1) xg]*betag');
% sigma=exp([ones(n,1) xs]*betas');
% nonlog=0;
% y=gprnd(gamma,sigma,0);
% %solution initiale ?
% [cst]=gpfit(y);
% p0 = [log(cst(2));zeros(3,1)+.01;log(cst(1));zeros(3,1)+.01];
% [pest,fval,exitflag,output,grad,hessian,sigma_est,gamma_est]=gpd_regression(p0,y,xs,xg,nonlog);

[pest,fval,exitflag,output,grad,hessian]=fminunc('gpdregress_full2',p0,optimset('hessian','on','algorithm','trust-region','largescale','on','gradobj','on','display','iter','DerivativeCheck','off'),y,xs,xg,nonlog); % initial unpenalized estimation
[~,~,~,sigma,gamma] = gpdregress_full2(pest,y,xs,xg,nonlog);


end