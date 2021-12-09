
L. Mhalla, J. Hambuckers and M. Lambert, "Extremal connectedness of hedge fund". Journal of Applied econometrics, forthcoming.

The data used in the empirical analyses have been provided by Hedge Funds Research (HFR) and are available at https://www.hfr.com/hfr-database under proper licensing and registration.

Marginal analysis
-----------------

Fund-specific variables have been provided by HFR. It includes the returns (i.e. our target variable, taken as the adjusted performance calculated by HFR), the asset under management (scaled by a factor of 1000), the level of incentive and management fees, and the existence of a leverage or not.  
HFR provided us also with the sub-strategy classification between the 12 investment styles considered in the the paper, as well as the country of registration. A detailed description of which funds are excluded from this initial sample (e.g. because of missing data, early inception date, etc.)
is given in Section 1 of the Supplementary material available at https://github.com/lindamhalla/Extremal-connectedness-of-hedge-funds.

To associate macroeconomic indicators to each fund, we use the country of registration. Macroeconomic variables (i.e. unemployment rate, short- and long-term interest rates, industrial production) have been obtained from the OECD website (see http://stats.oecd.org - where detailed definitions
of the time series are provided). Stock volatility has been calculated as monthly average of 360-day volatility of the national stock market index. Data for these calculations have been obtained from Bloomberg.

Finally, we use the publicly available 7 risk factors of David Hsieh available at https://faculty.fuqua.duke.edu/~dah7/HFData.htm, and the 5 factors available on Kenneth French website at https://mba.tuck.dartmouth.edu/pages/faculty/ken.french/data_library.html

Time series of the macro-financial variables used in this paper are made available in the file 'macro_finance.xlsx', Funds returns for which macro-financial data are missing are excluded (this concerns mostly small countries and a limited number of funds, since overall 95% of our observations are tied to funds registered in the US, the UK, Switzerland, France, and Canada, for which no macro-financial data are missing).
Factors data can be found in the files 'fung_hsieh.xlsx' and 'fama_french.xlsx'. 

Remark that all returns have been multiplied by -1 and we focused on the right tail, since our interest is on losses.

The code in 'script_marginal_analysis.m' allows estimating an extreme value regression model, either with a likelihood or a penalized-likelihood approach, as in the paper.


Multivariate analysis
--------------------

After estimating the marginal tail distributions of the different investment styles of thehedge funds, the negative returns are transformed to the unit-Frechet scale using eq. (10), as described at the beginning of Section 3.3 of the paper.

Then, (bivariate) vectors of transformed returns are built for each pair of investment styles. Finally, only observations with a radial component exceeding a high threshold are kept for inference, where a Husler--Reiss spectral density is fitted.

In the Rdata file 'df.multivariate.Rdata', we provide the negative returns transformed to the unit Frechet scale (for anonymity reasons), as well as the set of covariates used to model the dynamic tail connectedness measures, i.e., the VIX (retrieved from the Federal Reserve of Economic Data), the negative returns of the MSCI (retrieved from the MSCI World page), the EPU (retrieved from the Federal Reserve of Economic Data), and the FSI (retrieved from the Federal Reserve of Economic Data).

The R script used to estimate the extremal connectedness between investment styles of hedge funds can be found in the file 'modeling_extremal_connectedness_funds.R'. This R script allows also the computations of the ECoVaR for all pairs of strategies.
The functions needed to fit the Husler--Reiss spectral density can be found in 'spectral_HR.R'.
