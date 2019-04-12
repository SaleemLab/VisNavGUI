function pFisher = getFisherpvalue(pval)
    chi22 = -2*sum(log(pval));
    pFisher = chi2cdf(chi22,2*numel(pval),'upper');
end