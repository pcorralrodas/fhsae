*! fayherriot December, 2017 - muestra INEI
* Paul Corral (World Bank Group - Poverty and Equity Global Practice)
* William Seitz (World Bank Group - Poverty and Equity Global Practice)

cap program drop fayherriot
cap set matastrict off
program define fayherriot, eclass
	version 11.2
	#delimit ;
	syntax varlist (min=1 numeric fv) [if] [in],
	REvar(varlist max=1 numeric)
	[FHpredict(string)
	FHSEpredict(string)
	FHCVpredict(string)
	DSEpredict(string)
	DCVpredict(string)
	AREApredict(string)
	GAMMApredict(string)
	OUTsample
	NONEGative];
#delimit cr		
set more off

marksample touse11
tokenize `varlist'
local depvar `1'
macro shift
local indeps `*'

//TEMP VARS AND NAMES
tempname beta vcov tomata s2u
tempvar lev yhat res touse out

//Create variable for out of sample prediction
if("`outsample'"!=""){
	qui:gen `out' = `depvar'==.
	foreach x of local indeps{
		qui:replace `out' = 0 if `x' == .
	}
}



//Remove collinear variables
_rmcoll `indeps', forcedrop
local indeps `r(varlist)'



//regression
noi: dis in green "=============================================================================="
noi: dis in yellow "Fay-Herriot Model" 
noi: dis in green "=============================================================================="

noi: regress	`depvar'  `indeps' if `touse11'==1
scalar er2a=e(r2_a)
scalar er2=e(r2)
scalar ef = e(F)
scalar enum = e(N)

	
//Predict vectors
qui{
gen `touse'=e(sample)
predict `lev', hat
predict `yhat', xb
predict `res', res

local `tomata' lev res revar region depvar indeps

	foreach x of local `tomata'{
		mata:st_view(`x'=.,.,"``x''", "`touse'")	
	}
	
	//Send to mata for FH
	noi{
	mata: _Fh = _fh(depvar, indeps, revar, res, lev)
	mata: fhse     = *(_Fh[1,1])
	mata: fhcv     = *(_Fh[1,2]) 
	mata: dse      = *(_Fh[1,3]) 
	mata: dcv      = *(_Fh[1,4]) 
	mata: st_matrix("`beta'",*(_Fh[1,5]))
	mata: st_matrix("`vcov'",*(_Fh[1,6]))
	mata: fh      = *(_Fh[1,7])
	mata: area    = *(_Fh[1,8]) 
	mata: st_matrix("`s2u'", *(_Fh[1,9]))
	mata: gamma   = *(_Fh[1,10])
	}
	 
	 noi{
	 mat rownames `beta' = `indeps' _cons
	 mat `beta' = `beta''
	 mat rownames `vcov' = `indeps' _cons
	 mat colnames `vcov' = `indeps' _cons
	 }
	
	if("`outsample'"!=""){
		mata: st_view(_x1=.,., "`indeps'","`out'")
		noi:mata: _out = _fhout(_x1,*(_Fh[1,5]),*(_Fh[1,6]),*(_Fh[1,9]))
		mata: fh_out      = *(_out[1,1])
		mata: fhse_out    = *(_out[1,2])
		mata: fhcv_out    = *(_out[1,3])
	 }
	 
	 local predicts fh fhse fhcv dse dcv area gamma
	 foreach x of local predicts{
		if("``x'predict'"!=""){
			local nm: list sizeof `x'predict
			if (`nm'>1){
				display as error "Only one variable name for option `x'predict"
				error 202020202
			}
			gen ``x'predict' = .
			mata: st_store(.,tokens("``x'predict'"),"`touse'",`x')

			if(("`x'"=="fh"|"`x'"=="fhse"|"`x'"=="fhcv") & "`outsample'"!=""){
				mata: st_store(.,tokens("``x'predict'"),"`out'",`x'_out)
			}	
		}	 
	 }
	 	
	 ereturn post `beta' `vcov', depname(`depvar') esample(`touse')
	 noi display in yellow "GLS coefficients"
	 noi ereturn display
	 ereturn scalar sigma2u = `s2u'[1,1]	 
	 ereturn scalar r2_a = er2a
	 ereturn scalar r2   = er2
	 ereturn scalar F    = ef
	 ereturn scalar N    = enum
	 
}
end

mata

function _fh(y, x, sigma2, eps, lv){
	pointer(real matrix) rowvector _fhval
	_fhval = J(1,10,NULL)
	x=x,J(rows(x),1,1)
	dof1 = cols(x)
	
	noneg = (st_local("nonegative")~="")

	//The shrinkage component
	delta = quadsum(eps:^2) - quadsum((sigma2:*(1:-lv)))
	sigma2u = delta/(rows(x)-dof1)
	if (sigma2u<0){
		sigma2u=0
		V=sigma2
	}
	else{
		V = sigma2:+sigma2u
	}
	b_gls    = _fheblup(V,x,y)
	_beta    = *b_gls[1,1]
	_varbeta = *b_gls[1,2]
	yhat     =  quadcross(x',_beta)
	if (noneg==1) yhat = yhat:*(yhat:>0)	

	//shrinkage estimator
	gamma = sigma2u:/(V)
	//EBLUP
	FH = (gamma:*y)+((1:-gamma):*yhat)
	//Random area effects
	a_eff = gamma:*(y-yhat)
	
	//MSE-estimate
	g1 = gamma:*sigma2
		
	g2 = ((1:-gamma):^2):*quadrowsum((quadcross(x',(_varbeta)):*x))
		vsig = 2*rows(x)/(quadsum(1:/V)^2)
	g3 = vsig:*((sigma2:^2):/(V:^3))
	
	mse_fh = g1+g2+(2:*g3)

	//FHSE
	_fhval[1,1] = &(sqrt(mse_fh))
	//FHCV
	_fhval[1,2] = &(100:*(sqrt(mse_fh):/FH))
	//DIRECT SE
	_fhval[1,3] = &(sqrt(sigma2))
	//DIRECT CV
	_fhval[1,4] = &(100:*(sqrt(sigma2):/y))
	//Betas
	_fhval[1,5] = &(_beta)
	//VCOV
	_fhval[1,6] = &(_varbeta)
	//FH
	_fhval[1,7] = &(FH)
	//Area effect
	_fhval[1,8] = &(a_eff)
	//Sigma2u
	_fhval[1,9] = &(sigma2u)
	//gamma
	_fhval[1,10] = &(gamma)
	
	return(_fhval)
}

function _fheblup(_v,_x,_y){
	pointer(real matrix) rowvector fhout
	fhout = J(1,2,NULL)

	_varbeta = invsym(quadcross(_x,(1:/_v),_x))
	_beta    = quadcross(_varbeta,quadcross(_x,(1:/_v),_y))
	
	fhout[1,1] = &(_beta)
	fhout[1,2] = &(_varbeta)
	return(fhout)
}

function _fhout(_x1,_bgls,_vcovgls,sigma2u){
	pointer(real matrix) rowvector _fhout1
	_fhout1 = J(1,3,NULL)
	_x1 = _x1,J(rows(_x1),1,1)
	noneg = (st_local("nonegative")~="")
	
	yhat = quadcross(_x1',_bgls)
	
	if (noneg==1) yhat = yhat:*(yhat:>0)
	
	_fhout1[1,1] = &(yhat)	
	_fhout1[1,2] = &(sqrt(quadrowsum((quadcross(_x1',(_vcovgls)):*_x1)):+ sigma2u))
	_fhout1[1,3] = &(100*((*_fhout1[1,2]):/(*_fhout1[1,1])))	
		
	return(_fhout1)
}
end

