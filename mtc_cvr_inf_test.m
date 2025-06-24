function [mtc_diff, mtc_pval, cvr_diff, cvr_pval] = mtc_cvr_inf_test(b, B, rt1, rt0, alpha)

%==========================================================================
%==========================================================================
% mtc_cvr_inf_test test
%                        H0: MTC(rt1) - MTC(rt0) = 0
%                        H1: MTC(rt1) - MTC(rt0) > 0
% and
%                        H0: CVaR(rt1) - CVaR(rt0) = 0
%                        H1: CVaR(rt1) - CVaR(rt0) > 0
% and
%                        H0: SHR(rt1) - SHR(rt0) = 0
%                        H1: SHR(rt1) - SHR(rt0) > 0


% and report the p-values along with differences of MTCs and CVaRs


% Please cite the associated paper as in README file:



% Inputs:

% b: Length of block
% B: Number of bootstrap iterations
% rt1, rt0: Return (or excess return) time series with length of S from
% which we resample blocks of length b with replacement
% alpha: Confidence level for calculating CVaR and MTC



% Outputs:

% mtc_diff: MTC(rt1) - MTC(rt0)
% cvr_diff: CVaR(rt1) - CVaR(rt0)
% mtc_pval: The p-value of testing MTC differences
% cvr_pval: The p-value of testing CVaR differences
%==========================================================================
%==========================================================================
ret_matrx = [rt1 rt0];
S = size(ret_matrx,1);% lenght of time series rt1, rt0

%==========================================================================
%            Calculate CVaR(rt1), MTC(rt1), CVaR(rt0), MTC(rt0)           %
%==========================================================================
options = optimset('linprog');
options.Display = 'off';

% CVaR calculation based on Rockafellar-Uryasev formulation
[xcv1, fv1] = linprog([(1/S)*ones(S,1); (1-alpha)], ...
    [-eye(S) -ones(S,1)],rt1,[],[],[zeros(S,1);-inf], [inf*ones(S,1);inf],options);
[xcv0, fv0] = linprog([(1/S)*ones(S,1); (1-alpha)], ...
    [-eye(S) -ones(S,1)],rt0,[],[],[zeros(S,1);-inf], [inf*ones(S,1);inf],options);
cvr_1 = (fv1/(1-alpha));
cvr_0 = (fv0/(1-alpha));

mtc_1 = mean(rt1)/cvr_1;% this is MTC(rt1)
mtc_0 = mean(rt0)/cvr_0;% this is MTC(rt0)

%==========================================================================
%==========================================================================
%                       H0: MTC(rt1) - MTC(rt0) = 0                       %
%                       H1: MTC(rt1) - MTC(rt0) > 0                       %
%==========================================================================
%==========================================================================

%==========================================================================
% Steps 1 and 2 of inference test algorithm: construct bootstraped sample

% generate randomized starting points for overlapped blocks (rsp_k)
rsb0 = S-b+1;% total number of block
for m = 1:B
    rsp_all(:,m) = randi(rsb0,rsb0,1);% generate rsb0 randomized starting points
end

k = ceil(S/b);% required number of blocks to cover a bootstraped sample of length S
for m = 1:B
    rsp_k(:,m) = rsp_all(randi(size(rsp_all,1),k,1),m);% select randomly k starting points from rsp_all
end

% construct bootstraped sample using randomized starting points in rsp_k
for m = 1:B
    rt1boot = [];% bootstrapped sample of rt1
    rt0boot = [];% bootstrapped sample of rt0
    
    for i=1:size(rsp_k,1)
        
        % extract the blocks of length b using block starting points in rsp_k
        minv = rsp_k(i,m);
        maxv = rsp_k(i,m)+b-1;
        sample = linspace(minv,maxv,b)';
        
        % allocate these elements to the bootstrapped series
        rt1boot = [rt1boot; ret_matrx(sample,1)];
        rt0boot = [rt0boot; ret_matrx(sample,2)];
        
    end
    
    rt1star(:,m) = rt1boot(1:S,1);% this is bootstraped sample for rt1
    rt0star(:,m) = rt0boot(1:S,1);% this is bootstraped sample for rt0
    
end

%==========================================================================
% Step 3 of inference test algorithm: calculating \bar{Y}^*

for m = 1:B
    
    ret_1 = rt1star(:,m);
    ret_0 = rt0star(:,m);
    
    [xcv0, fv0] = linprog([(1/S)*ones(S,1);(1-alpha)], ...
        [-eye(S) -ones(S,1)],ret_0,[],[],[zeros(S,1);-inf], [inf*ones(S,1);inf],options);
    [xcv1, fv1] = linprog([(1/S)*ones(S,1);(1-alpha)], ...
        [-eye(S) -ones(S,1)],ret_1,[],[],[zeros(S,1);-inf], [inf*ones(S,1);inf],options);
    CVaR_0 = (1/(1-alpha))*fv0;
    CVaR_1 = (1/(1-alpha))*fv1;
    
    YbStar(m,:) = [mean(ret_1)  CVaR_1  mean(ret_0)  CVaR_0];% calculate \bar{Y}^*
    
end

%==========================================================================
% Step 4 of inference test algorithm: calculate \hat{\tau}_0^2

mtc_diff = mtc_1 - mtc_0;% MTC differences
mtcdiff_s = sqrt(S)*(mtc_1 - mtc_0);% test statistics

Mbar = (1/B)*sum(YbStar);
SigmaStar_r = zeros(4,4);
for m = 1:B
    SigmaStar_r = SigmaStar_r + (1/B)*((YbStar(m,:) - Mbar)'*(YbStar(m,:) - Mbar));
end

c_und = [1/cvr_1;-0.5*((mtc_1+mtc_0)/cvr_1);-1/cvr_0;0.5*((mtc_1+mtc_0)/cvr_0)];
tau2 = S*c_und'*SigmaStar_r*c_und;

%==========================================================================
% Step 5 of inference test algorithm: calculate p-value

% p-value calculation
if abs(mtc_1 - mtc_0) <= 1e-6
    mtc_pval = 1;
else
    mtc_pval = 1-normcdf(mtcdiff_s, 0, sqrt(tau2));%
end

%==========================================================================
%==========================================================================
%                     H0: CVaR(rt1) - CVaR(rt0) = 0                       %
%                     H1: CVaR(rt1) - CVaR(rt0) > 0                       %
%==========================================================================
%==========================================================================

%==========================================================================
% Steps 1-3 are the same as MTC differences test
%==========================================================================
% Step 4 of inference test algorithm: calculate \nu{\tau}_0^2

cvr_diff = cvr_1 - cvr_0;% CVaR differences
cvrdiff_s = sqrt(S)*(cvr_1 - cvr_0);% test statistics

e_und = [0;1;0;-1];
nu2 =  S*e_und'*SigmaStar_r*e_und;

%==========================================================================
% Step 5 of inference test algorithm:  calculate p-value

% p-value calculation
if abs(cvr_1 - cvr_0) <= 1e-6
    cvr_pval = 1;
else
    cvr_pval = 1-normcdf(cvrdiff_s, 0, sqrt(nu2));
end
%==========================================================================

end
