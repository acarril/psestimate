{smcl}
{* *! version 1.0.1 May 2016}{...}
{cmd:help psestimate}{...}
{hline}

{title:Title}

{pstd}
{hi:psestimate} {hline 2} Estimate the propensity score proposed by Imbens and Rubin (2015)

{title:Syntax}

{p 8 16 2}
{cmd:psestimate} {depvar} [{indepvars}]
[, {opt t:otry(varlist)} {opt cl:inear(real)} {opt cq:uadratic(real)}]
{p_end}

{synoptset 20 tabbed}{...}
{synopthdr}
{synoptline}
{synoptset 20 tabbed}{...}
{synopt:{opt t:otry}} specify list of covariates to try; default is all{p_end}
{synopt:{opt cl:inear}} threshold value for likelihood ratio test of first order covariates; default is 1{p_end}
{synopt:{opt cq:uadratic}} threshold value for likelihood ratio test of second order covariates; default is 2.71{p_end}
{p 4 6 2}

{title:Description}

{pstd}
The {cmd:psestimate} command estimates the propensity score proposed by {help psestimate##imbens_rubin_2015:Imbens and Rubin (2015)}.
In particular, it implements the algorithm outlined by {help psestimate##imbens_2015:Imbens (2015)}, which estimates the propensity score for a binary dependent variable indicating treatment status.
A subset of covariates may be explicitly included in the linear part of the specification as independent variables.

{pstd}
The algorithm chooses first order terms from all remaining variables of the dataset, unless a subset of variables is specified with the {opt t:otry(varlist)} option.
The selection of first order terms is performed in a stepwise fashion, comparing the base (nested) model to a model with one single additional covariate.
A likelihood ratio test (LRT; see {help lrtest}) to test the null hypothesis of the non-significance of the additional coefficient is performed.
A comparison is made for all remaining covariates, selecting the covariate associated with the highest LRT statistic.
This covariate is then included in the model.

{pstd}
The process of selecting and including additional linear terms is carried out until the highest LRT statistic is less than {opt cl:inear(real)} or there are no remaining covariates to add.
Then the algorithm chooses second order terms only from covariates selected for the linear specification, performing analogue tests for selection among all interactions and quadratic terms.
This second process of selecting and including additional quadratic terms is carried out until the highest LRT statistic is less than {opt cq:uadratic(real)} or there are no remaining covariates to add.

{title:Options}

{phang}
{opt t:otry(varlist)} specifies the vector of covariates to try into the first and second order terms.
The default is to include all variables in the dataset, exluding the {depvar} and other base model covariates indicated in {indepvars} (if any).

{phang}
{opt cl:inear(real)} specifies the threshold value used for the addition of first order (linear) terms.
The decision is based on a likelihood ratio test statistic for the null hypothesis that the coefficient of the additional first order term is equal to zero.
Default value is 1.

{phang}
{opt cq:uadratic(real)} specifies the threshold value used for the addition of second order (quadratic) terms.
The decision is based on a likelihood ratio test statistic for the null hypothesis that the coefficient of the additional second order term is equal to zero.
Default value is 2.71.

{title:Examples}

{title:Author}

{pstd}
Alvaro Carril{break}
Research Analyst at J-PAL LAC{break}
acarril@fen.uchile.cl

{title:Acknowledgements}

{pstd}
This program started as a refinement of Juan Ignacio Elorrieta's work for Bustos, Pomeranz and Zucman (forthcoming).
Juan Ignacio's code provided a significant headstart, while comments from Dina Pomeranz and Sebastián Bustos helped to fine tune the program. All remaining errors are my own.

{title:References}
{marker imbens_rubin_2015}
{phang}Imbens, Guido W. and Donald B. Rubin. 2015. {it: Causal Inference in Statistics, Social, and Biomedical Sciences}. New York: Cambridge University Press.
{marker imbens_2015}
{phang}Imbens, Guido W. 2015. "Matching Methods in Practice: Three Examples". {it:Journal of Human Resources} 50(2): 373-419.
