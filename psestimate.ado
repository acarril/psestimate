*! 1.4 Alvaro Carril 20nov2016
program define psestimate, rclass
	version 11
	
syntax varlist(min=1 fv) [if] [in] [, ///
	Totry(varlist fv) ///
	NOTry(varlist) ///
	CLinear(real 1) ///
	CQuadratic(real 2.71) ///
	ITERate(passthru) ///
	GENPScore(name) ///
	GENLor(name) ///
	noLin ///
	noQuad ///
	]
	
marksample touse

*-------------------------------------------------------------------------------
* Inputs
*-------------------------------------------------------------------------------
* Factor variables:
local fvops = "`s(fvops)'" == "true" | _caller() >= 11
if `fvops' {
	gettoken lhs rest : varlist // extract depvar
	_fv_check_depvar `lhs' // check depvar is not a factor variable
}

* Checks:
foreach g in `genpscore' `genlor' {
	if "`g'" != "" confirm new var `g'
}

if ("`lin'" == "nolin" & "`quad'" == "noquad") {
	display as error "options nolin and noquad may not be combined"
	exit 198
}

* Try all variables not defined in varlist
if missing("`totry'") {
	qui ds `varlist' `notry' __00*, not
	local totry `r(varlist)'
}
else {
	local totry : list totry - varlist
	local totry : list totry - notry
}

* Extract treatment variable and base covariates from varlist
local treatvar :	word 1 of `varlist'
local K_b :			list varlist - treatvar

* Thresholds:
local C_lin			`clinear'
local C_qua			`cquadratic'

*-------------------------------------------------------------------------------
* Initial setup
*-------------------------------------------------------------------------------
local h `K_b' `K_l' `K_q' // generic vector of functions
local llrt_max = `C_lin' // set equal to linear threshold to start while loop

* Estimate base model:
qui logit `treatvar' `h' if `touse', `iterate'
estimates store null

*-------------------------------------------------------------------------------
* Select first order covariates (steps 1-5)
*-------------------------------------------------------------------------------

if "`lin'" != "nolin" { 
	* Indicate progress of first order covaraites loop:
	local N_foc : list sizeof totry
	local N_foc = `N_foc'*(`N_foc'+1)/2
	nois _dots 0, reps(`N_foc') title(Selecting first order covariates...)
	local rep 1
	
	* Start first order covariates loop
	while `llrt_max' >= `C_lin' {
		local llrt_max = `C_lin'
		if !missing("`totry'") {
			foreach v of local totry {
				local estrep = `estrep'+1
				capture quietly logit `treatvar' `h' `v' if `touse', `iterate'
				if _rc == 0 {
					estimates store est`estrep'
					qui lrtest null est`estrep', force
					if (`r(chi2)' >= `llrt_max') {
						local v_max `v' // store covariate with max llrt stats
						local llrt_max = `r(chi2)' // update maximum llrt stat
						capture estimates drop v_max
						qui estimates restore est`estrep'
						estimates store v_max
						estimates drop est`estrep'
					}
				}
				local N_foc : list sizeof totry
				if `estrep' != `N_foc' {
					nois _dots `rep++' 0
				}
				else {
					if "`v_max'" != "" {
						nois _dots `rep++' -1
					}
				}
			}
		}
		if "`v_max'" != "" {
			qui estimates restore v_max // restore computed estimates for selected covariate
			estimates clear // clear all other estimates
			estimates store null // update null model estimates with the selected covariate
			local K_l `K_l' `v_max'
			local h `K_b' `K_l' `K_q'
			local totry: list totry - v_max
			local v_max
			local estrep 0
		}
		else {
			if !missing("`K_l'") {
				di as text _newline "Selected first order covariates are: " as result "`K_l'"
			}
			else di as text _newline "No first order covariates selected"
			estimates drop _all
			continue, break
		}
	}
}

*-------------------------------------------------------------------------------
* Select second order covariates (steps 6-10)
*-------------------------------------------------------------------------------
if "`quad'" != "noquad" { 

* Separate lists of vars and fvvars
*-------------------------------------------------------------------------------
	foreach v in `h' {
		capture _fv_check_depvar `v'
		if _rc {
			local h_cat `h_cat' `v'
			local h_fv `h_fv' `v'
		}
		else {
			local h_nocat `h_nocat' c.`v'
			local h_fv `h_fv' c.`v'
		}
	}

* Generate interactive variables from linear model
*-------------------------------------------------------------------------------
	local num_h : word count `h_fv'
	local totry // clear totry varlist
	forval i = 1/`num_h' {
		forval j = 1/`=`i'-1' {
			local x : word `i' of `h_fv'
			local y : word `j' of `h_fv'
			local totry `totry' `x'#`y'
		}
	}
	
	di "totry:  `totry'"
	
* Generate quadratic varaibles from linear model
*-------------------------------------------------------------------------------

	* Collect dummies
	qui ds `h', has(type numeric)
	local h_numeric `r(varlist)'
	foreach v of local h_numeric {
		capture assert missing(`v') | inlist(`v', 0, 1)
		if _rc != 0 local nondummy `nondummy' `v'
	}
	
	* Add quadratic terms of non-dummies
	foreach z of local nondummy {
		local totry `totry' c.`z'#c.`z'
	}

	local quadvars `totry' // preserve list of all quadratic terms to try

* Select second order terms
*-------------------------------------------------------------------------------
	* Estimate base model again:
	qui logit `treatvar' `h' if `touse', `iterate'
	estimates store null
	
	* Indicate progress of second order covaraites loop:
	local N_soc : list sizeof totry
	local N_soc = `N_soc'*(`N_soc'+1)/2
	nois _dots 0, reps(`N_soc') title(Selecting second order covariates...)
	local rep 1
	
	* Start second order covariates loop
	local llrt_max = `C_qua'
	local estrep = 0
	while `llrt_max' >= `C_qua' {
		local llrt_max = `C_qua'
		if !missing("`totry'") {
			foreach v in `totry' {
				local estrep = `estrep'+1
				capture quietly logit `treatvar' `h' `v' if `touse', `iterate'
				if _rc == 0 {
					estimates store est`estrep'
					qui lrtest null est`estrep', force
					if (`r(chi2)' >= `llrt_max') {
						local v_max `v' // store covariate with max llrt stats
						local llrt_max = `r(chi2)' // update maximum llrt stat
						capture estimates drop v_max
						qui estimates restore est`estrep'
						estimates store v_max
						estimates drop est`estrep'
					}
				}
				local N_soc : list sizeof totry
				if `estrep' < `N_soc' {
					nois _dots `rep++' 0
				}
				else {
					if "`v_max'" != "" {
						nois _dots `rep++' -1
					}
				}
			}
		}
		if "`v_max'" != "" {
			qui estimates restore v_max // restore computed estimates for selected covariate
			estimates clear // clear all other estimates
			estimates store null // update null model estimates with the selected covariate
			local K_q `K_q' `v_max'
			local h `K_b' `K_l' `K_q'
			local totry: list totry - v_max
			local v_max
			local estrep 0
		}
		else {
			if !missing("`K_q'") {
				di as text _newline "Selected second order covariates are: " as result "`K_q'"
			}
			else di as text _newline "No second order covariates selected"
			estimates drop _all
			continue, break
		}
	}
}
* Show final model
di as text "Final model is: " as result "`h'"

* Save return results
return local h `h'
return local K_q `K_q'
return local K_l `K_l'
return local K_b `K_b'
return local tvar `treatvar'
return scalar C_q = `C_qua'
return scalar C_l = `C_lin'

* Estimate final model to save eresults
qui logit `treatvar' `h' if `touse'
* Generate PS hat and generate log odds ratio
tempvar `genpscore' `genlor' ps
qui predict `ps' if e(sample) == 1, pr
if "`genlor'" != "" {
	qui gen `genlor' = ln(`ps'/(1-`ps')) if `touse'
	lab var `genlor' "Log odds ratio"
	if "`genpscore'" != "" {
		qui rename `ps' `genpscore'
		lab var `genpscore' "Propensity score"
		order `genpscore' `genlor', last
	}
}
else {
	if "`genpshat'" != "" {
		qui rename `ps' `genpscore'
		lab var `genpscore' "Propensity Score"
	}
}

end
