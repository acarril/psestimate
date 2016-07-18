*! 1.2 Alvaro Carril 18jul2016
program define psestimate, rclass
	version 11
	
syntax varlist(min=1) [if] [in] [, ///
	Totry(varlist) ///
	NOTtry(varlist) ///
	CLinear(real 1) ///
	CQuadratic(real 2.71) ///
	ITERate(passthru) ///
	GENPShat(name) ///
	GENLor(name) ///
	noLin ///
	noQuad ///
	]

marksample touse
*-------------------------------------------------------------------------------
* Inputs
*-------------------------------------------------------------------------------
* Checks:
foreach g in `genpshat' `genlor' {
	if "`g'" != "" confirm new var `g'
}

if ("`lin'" == "nolin" & "`quad'" == "noquad") {
	display as error "options nolin and noquad may not be combined"
	exit 198
}

* Extract treatment variable and base covariates from varlist
local treatvar :	word 1 of `varlist'
local K_b :			list varlist - treatvar

* Try all variables not defined in varlist
if missing("`totry'") {
	qui ds
	local totry `r(varlist)'
}
local totry :	list totry - varlist
local totry :	list totry - nottry

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
	nois _dots 0, reps(`N_foc') title(Selecting first order covariates...)
	local rep 1
	
	* Start first order covariates loop
	while `llrt_max' >= `C_lin' {
		local llrt_max = `C_lin'
		if !missing("`totry'") {
			foreach v of varlist `totry' {
				capture quietly logit `treatvar' `h' `v' if `touse', `iterate'
				if _rc == 0 {
					estimates store `v'
					qui lrtest null `v', force
					if (`r(chi2)' >= `llrt_max') {
						local v_max `v' // store covariate with max llrt stats
						local llrt_max = `r(chi2)' // update maximum llrt stat
					}
				}
				local N_soc : list sizeof totry
				if `estrep' != `N_soc' {
					nois _dots `rep++' 0
				}
				else {
					nois _dots `rep++' -1
				}
			}
		}
		if "`v_max'" != "" {
			qui estimates restore `v_max' // restore computed estimates for selected covariate
			estimates clear // clear all other estimates
			estimates store null // update null model estimates with the selected covariate
			local K_l `K_l' `v_max'
			local h `K_b' `K_l' `K_q'
			local totry: list totry - v_max
			local v_max
		}
		else {
			di as text _newline "Selected first order covariates are: " as result "`K_l'"
			estimates drop _all
			continue, break
		}
	}
}

*-------------------------------------------------------------------------------
* Select second order covariates (steps 6-10)
*-------------------------------------------------------------------------------
if "`quad'" != "noquad" { 
* Generate interactive variables from linear model
*-------------------------------------------------------------------------------
	local num_h : word count `h'
	local totry // clear totry varlist
	forval i = 1/`num_h' {
		forval j = 1/`=`i'-1' {
			local x : word `i' of `h'
			local y : word `j' of `h'
			local totry `totry' c.`x'#c.`y'
		}
	}

di "TOTRY: `totry'"
	
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
		qui gen `z'_2 = `z'^2
		local labz : variable label `z'
		if !missing("`labz'") label var `z'_2 `"`labz' squared"'
		local totry `totry' `z'_2
	}

	local quadvars `totry' // preserve list of all quadratic terms to try

* Select second order terms
*-------------------------------------------------------------------------------
	* Estimate base model again:
	qui logit `treatvar' `h' if `touse', `iterate'
	estimates store null
	
	* Indicate progress of second order covaraites loop:
	local N_soc : list sizeof totry
	nois _dots 0, reps(`N_soc') title(Selecting second order covariates...)
	local rep 1
	
	* Start second order covariates loop
	local llrt_max = `C_qua'
	while `llrt_max' >= `C_qua' {
		local llrt_max = `C_qua'
		if !missing("`totry'") {
			foreach v in `totry' {
				local estrep = `rep'+1
				capture quietly logit `treatvar' `h' `v' if `touse', `iterate'
				if _rc == 0 {
					estimates store est`estrep'
					qui lrtest null est`estrep', force
					if (`r(chi2)' >= `llrt_max') {
						local v_max `v' // store covariate with max llrt stats
						local llrt_max = `r(chi2)' // update maximum llrt stat
					}
				}
				local N_soc : list sizeof totry
				if `estrep' != `N_soc' {
					nois _dots `rep++' 0
				}
				else {
					nois _dots `rep++' -1
				}
			}
		}
		if "`v_max'" != "" {
			qui estimates restore est`estrep' // restore computed estimates for selected covariate
			estimates clear // clear all other estimates
			estimates store null // update null model estimates with the selected covariate
			local K_q `K_q' `v_max'
			local h `K_b' `K_l' `K_q'
			local totry: list totry - v_max
			local v_max
			local estrep
		}
		else {
			di as text _newline "Selected second order covariates are: " as result "`K_q'"
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
tempvar `genpshat' `genlor' ps
qui predict `ps' if e(sample) == 1, pr
if "`genlor'" != "" {
	qui gen `genlor' = ln(`ps'/(1-`ps')) if `touse'
	lab var `genlor' "Log odds ratio"
	if "`genpshat'" != "" {
		qui rename `ps' `genpshat'
		lab var `genpshat' "Propensity score"
		order `genpshat' `genlor', last
	}
}
else {
	if "`genpshat'" != "" {
		qui rename `ps' `genpshat'
		lab var `genpshat' "Propensity Score"
	}
}
end
