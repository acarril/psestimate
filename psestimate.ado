*! 1.0.1 Alvaro Carril 07may2016
program define psestimate
	version 12
	
syntax varlist(min=1) [, Totry(varlist) CLinear(real 1) CQuadratic(real 2.71)]

*-------------------------------------------------------------------------------
* Inputs
*-------------------------------------------------------------------------------

* Extract treatment variable and base covariates from varlist
local treatvar :	word 1 of `varlist'
local K_b :			list varlist - treatvar

* Try all variables not defined in varlist
if missing("`totry'") {
	qui ds
	local totry `r(varlist)'
}
local totry :		list totry - varlist

* Thresholds:
local C_lin			`clinear'
local C_qua			`cquadratic'

*-------------------------------------------------------------------------------
* Initial setup
*-------------------------------------------------------------------------------
local h `K_b' `K_l' `K_q' // generic vector of functions
local llrt_max = `C_lin' // set equal to linear threshold to start while loop

*-------------------------------------------------------------------------------
* Select first order covariates (steps 1-5)
*-------------------------------------------------------------------------------
* Estimate base model:
qui logit `treatvar' `h'
estimates store null

* Indicate progress of first order covaraites loop:
local N_foc : list sizeof totry
nois _dots 0, reps(`N_foc') title(Selecting first order covariates...)
local rep 1

while `llrt_max' >= `C_lin' {
	local llrt_max = `C_lin'
	foreach v of varlist `totry' {
		qui logit `treatvar' `h' `v'
		estimates store alternative_`v'
		qui lrtest null alternative_`v', force
		if (`r(chi2)' >= `llrt_max') {
			local v_max `v' // store covariate with max llrt stats
			local llrt_max = `r(chi2)' // update maximum llrt stat
		}
	nois _dots `rep++' 0
	}
	if "`v_max'" != "" {
		qui estimates restore alternative_`v_max' // restore computed estimates for selected covariate
		estimates clear // clear all other estimates
		estimates store null // update null model estimates with the selected covariate
		local K_l `K_l' `v_max'
		local h `K_b' `K_l' `K_q'
		local totry: list totry - v_max
		local v_max
		local success = -1 // update success for progress bar
		nois _dots `rep++' `success' // update progress bar
	}
	else {
		di as text _newline "Selected first order covariates are: " as result "`K_l'"
		estimates drop alternative_*
		continue, break
	}
}

*-------------------------------------------------------------------------------
* Select second order covariates (steps 6-10)
*-------------------------------------------------------------------------------

* Generate interactive variables from linear model
*-------------------------------------------------------------------------------
local num_h : word count `h'
local totry // clear totry varlist
forval i = 1/`num_h' {
  forval j = 1/`=`i'-1' {
    local x : word `i' of `h'
    local y : word `j' of `h'
    generate `x'X`y' = `x' * `y'
	local labx : variable label `x'
	local laby : variable label `y'
	if (!missing("`labx'") & !missing("`laby'")) label var `x'X`y' `"`labx' X `laby'"'
	local totry `totry' `x'X`y'
  }
}

* Generate quadratic varaibles from linear model
*-------------------------------------------------------------------------------
* Collect dummies
qui ds `h', has(type numeric)
local h_numeric `r(varlist)'
foreach v of local h_numeric {
	capture assert missing(`v') | inlist(`v', 0, 1)
	if _rc != 0 local nondummy `nondummy' `v'
}

foreach z of local nondummy {
	gen `z'_2 = `z'^2
	local labz : variable label `z'
	if !missing("`labz'") label var `z'_2 `"`labz' squared"'
	local totry `totry' `z'_2
}

local quadvars `totry' // preserve list of all quadratic terms to try

* Select second order terms
*-------------------------------------------------------------------------------
* Indicate progress of first order covaraites loop:
local N_foc : list sizeof totry
nois _dots 0, reps(`N_foc') title(Selecting first order covariates...)
local rep 1

local llrt_max = `C_qua'
while `llrt_max' >= `C_qua' {
	local llrt_max = `C_qua'
	foreach v of varlist `totry' {
		qui logit `treatvar' `h' `v', vce(robust)
		estimates store alternative_`v'
		qui lrtest null alternative_`v', force
		if (`r(chi2)' >= `llrt_max') {
			local v_max `v' // store covariate with max llrt stats
			local llrt_max = `r(chi2)' // update maximum llrt stat
		}
	nois _dots `rep++' 0
	}
	if "`v_max'" != "" {
		qui estimates restore alternative_`v_max' // restore computed estimates for selected covariate
		estimates clear // clear all other estimates
		estimates store null // update null model estimates with the selected covariate
		local K_q `K_q' `v_max'
		local h `K_b' `K_l' `K_q'
		local totry: list totry - v_max
		local v_max
		local success = -1 // update success for progress bar
		nois _dots `rep++' `success' // update progress bar
	}
	else {
		di as text _newline "Selected second order covariates are: " as result "`K_q'"
		predict ps if e(sample) == 1, pr // creates PS hat
		gen log_odds = ln(ps / (1 - ps)) // creates log odds ratio
		lab var log_odds "Log odds ratio"
		local droplist: list quadvars - K_q
		drop `droplist'
		estimates drop _all
		continue, break
	}
}
di as text "Final model is: " as result "`h'"

end
