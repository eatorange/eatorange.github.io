*	Applications of Cisse and Barrett (2018) approaches



	*	Data
	*	It is a part of Ethiopian household panel data with 4 outcomes (all exp, food exp, HDDS and TLU)
		*use	"", clear
	
		
	
	*	Variation in underlying distribution attumption of an outcome variable
	*	The following codes generate three Stata programs generating C&B measure under 3 different assumption (Normal, Gamma and Inverse Gaussian)
	
			
		*	Each program takes 4 arguments
			*	1.	Outcome variable
			*	2.	Lagged outcome variables, up to certain order
			*	3.	name of outcome variable (serves as a prefix of variables created. Useful when variable name is large
			*	4.	threshold value
			
		
		*	Set global macros for the arguments above
			
			*	All Expenditure
			global	depvar_aexp		lnrexpaeu_peryear
			global	ldepvar_aexp	cl.${depvar_aexp}##cl.${depvar_aexp}	// This example uses up to 2nd order, but researchers can choose which order to use.
			global	depvarname_aexp	aexp
			global	threshold_aexp	poverty_line
			
			*	Food Expenditure
			global	depvar_fexp		lnrexpaeu_peryear
			global	ldepvar_fexp	cl.${depvar_fexp}##cl.${depvar_fexp}	// This example uses up to 2nd order, but researchers can choose which order to use.
			global	depvarname_fexp	fexp
			global	threshold_fexp	poverty_line
			
			*	HDDS
			global	depvar_HDDS		HDDS
			global	ldepvar_HDDS	cl.${depvar_HDDS}##cl.${depvar_HDDS}	// This example uses up to 2nd order, but researchers can choose which order to use.
			global	depvarname_HDDS	HDDS
			global	threshold_HDDS	HDDS_threshold
			
			*	TLU (Inverse hyperbolic transformed)
			global	depvar_TLU		tlu_IHS
			global	ldepvar_TLU		cl.${depvar_TLU}##cl.${depvar_TLU}	// This example uses up to 2nd order, but researchers can choose which order to use.
			global	depvarname_TLU	TLU
			global	threshold_TLU	tlu_IHS_threshold
		
		*	(required) Specify Stata regression setting
		*	These include control variables, robust cluster standard error, etc.
			global	controls	lnhead_age	malehead	headnoed	occupation_non_farm	 i.year
			global	option		cluster(village_num)	//	Regression options (Specification that goes after comma. Leave it as blank if no option is used.
		

		*	Normal distrubtion
			
			cap	program	drop	resil_normal
			program	resil_normal
				args	depvar	ldepvar	depvarname	threshold	//	Outcome variable, name of outcome variable, and threshold value
						
				*	Step 1: Conditional mean
				reg	`depvar'	`ldepvar'	${controls},	${option}
				est	sto	m1_`depvarname'_normal
				
					*	Sample used to construct conditional mean
					cap	drop	`depvarname'_sample_normal
					gen	`depvarname'_sample_normal=1	if	e(sample)	
					
					*	Predict conditional mean
					predict	double	mean_`depvarname'_normal,	xb
					lab	var	mean_`depvarname'_normal	"Conditional mean (Normal)"
					
					*	Construct residual square
					predict	double	e_`depvarname'_normal,	residuals
					gen	e_`depvarname'_normal_sq	=	(e_`depvarname'_normal)^2	
					
				
				*	Step 2: Conditional variance
				reg	e_`depvarname'_normal_sq	`ldepvar'	${controls}	if	`depvarname'_sample_normal==1,	${option}
				est	store	m2_`depvarname'_normal
				
					*	Predict conditional variance
					predict	var_`depvarname'_normal, xb
					gen		var_`depvarname'_normal_abs	=	abs(var_`depvarname'_normal)	//	Absolute value of predicted variance (predicted value can be negative, so should use absolute value)	
					gen		sd_`depvarname'_normal	=	sqrt(var_`depvarname'_normal_abs)	//	standard deviation
				
				*	Step 3: Construct the resilience measure
					
					*	Construct Z-score of normal distribution
					gen thresh_`depvarname'_normal=(`threshold'-mean_`depvarname'_normal)/sd_`depvarname'_normal	// Z-score
					
					*	Construct probability variables
					gen prob_below_`depvarname'=normal(thresh_`depvarname'_normal)	//	normal CDF
					gen `depvarname'_resil_normal		=	1-prob_below_`depvarname'	//	C&B measure (1-CDF)
					gen	`depvarname'_resil_normal_scale	=	`depvarname'_resil_normal	*	100	// C&B measure, scaled from 0 to 100
				
					lab	var	`depvarname'_resil_normal		"`depvarname' Resilience (normal)"
					lab	var	`depvarname'_resil_normal_scale	"`depvarname' Resilience (normal scaled)"
					
					*	Regression of C&B on controls and cluster
					reg	`depvarname'_resil_normal	`ldepvar'	${controls},	${option}
					est	store	m3_`depvarname'_normal
				
				summ	mean_`depvarname'_normal	var_`depvarname'_normal	thresh_`depvarname'_normal	prob_below_`depvarname'	`depvarname'_resil_normal	`depvarname'_resil_normal_scale	
				
				
			end
				
		
			*	Gamma distribution
			*	Note: In this code, step 2 (regressing squared residual) is done with OLS, since glm might fail to converge			
		
			cap	program	drop	resil_gamma
			program	resil_gamma
				args	depvar	ldepvar	depvarname	threshold	//	Outcome variable, name of outcome variable, and thresholds.
				
				*	Step 1: generate conditional mean
				glm	`depvar'	`ldepvar'	${controls},	${option}	family(gamma)

					*	Sample used to construct conditional mean
					cap	drop	`depvarname'_sample_gamma
					gen	`depvarname'_sample_gamma=1	if	e(sample)	
					
					*	Predict conditional mean
					predict	double	mean_`depvarname'_gamma
					lab	var	mean_`depvarname'_gamma	"Conditional mean (Gamma)"
					
					*	Construct residual square
					predict	double	e_`depvarname'_gamma,	r
					gen	e_`depvarname'_gamma_sq	=	(e_`depvarname'_gamma)^2			
				
					*	Store marginal effects 
					eststo	m1_`depvarname'_gamma:	margins, dydx(*)	post
				
				*	Step 2: generate conditional variance
				reg	e_`depvarname'_gamma_sq	`ldepvar'	${controls}	if	`depvarname'_sample_gamma==1,	${option}
				est	store	m2_`depvarname'_gamma
				
					*	Predict conditional variance
					predict	var_`depvarname'_gamma, xb
					gen		var_`depvarname'_gamma_abs	=	abs(var_`depvarname'_gamma)		//	Absolute value of predicted variance (predicted value can be negative, so should use absolute value)
					lab	var	var_`depvarname'_gamma_abs	"Conditional variance (Gamma)"
					gen		sd_`depvarname'_gamma	=	sqrt(abs(var_`depvarname'_gamma))	//	Take square root of absolute value, since predicted value can be negative which does not have square root (thus will be missing)
				
				*	Step 3: Construct the resilience measure
					
					*	Contruct Gamma distribution parameters
					gen alpha1_`depvarname'_gamma	= (mean_`depvarname'_gamma)^2 / var_`depvarname'_gamma_abs	//	shape parameter of Gamma (alpha)
					gen beta1_`depvarname'_gamma	= var_`depvarname'_gamma_abs / mean_`depvarname'_gamma		//	scale parameter of Gamma (beta)

					*	Construct C&B measure
					gen `depvarname'_resil_gamma	= gammaptail(alpha1_`depvarname'_gamma, `threshold'/beta1_`depvarname'_gamma)	//	C&B measure (1-Gamma CDF)
					gen	`depvarname'_resil_gamma_scale	=	`depvarname'_resil_gamma*100	//	C&B scaled from 0 to 100
				
					lab	var	`depvarname'_resil_gamma			"`depvarname' Resilience (gamma)"
					lab	var	`depvarname'_resil_gamma_scale		"`depvarname' Resilience (gamma scaled)"
				
					*	Regression of C&B on controls and cluster
					reg	`depvarname'_resil_gamma	`ldepvar'	${controls},	${option}
					est	store	m3_`depvarname'_gamma
				
					*	Summary stats
					summ	mean_`depvarname'_gamma	var_`depvarname'_gamma_abs	alpha1_`depvarname'_gamma	beta1_`depvarname'_gamma	`depvarname'_resil_gamma	`depvarname'_resil_gamma_scale

			end
			
			*	Inverse Gaussian distribution
			*	Note: In this code, step 2 (regressing squared residual) is done with OLS, since glm might fail to converge
		
			cap	program	drop	resil_igaussian
			program	resil_igaussian
				args	depvar	ldepvar	depvarname	threshold	//	Outcome variable, name of outcome variable, and thresholds.
			
				*	Step 1: Conditional mean
				glm	`depvar'	`ldepvar'	${controls},	${option} family(igaussian)
					
					*	Sample used to construct conditional mean
					cap	drop	`depvarname'_sample_igaussian
					gen	`depvarname'_sample_igaussian=1	if	e(sample)	
					
					*	Predict conditional mean
					predict	double	mean_`depvarname'_igaussian
					lab	var	mean_`depvarname'_igaussian	"Conditional mean (Inv Gaussian)"
					
					*	Construct residual square
					predict	double	e_`depvarname'_igaussian,	r
					gen	e_`depvarname'_igaussian_sq	=	(e_`depvarname'_igaussian)^2			
				
					*	Store marginal effects 
					eststo	m1_`depvarname'_igaussian:	margins, dydx(*)	post
						
				*	Step 2: generate conditional variance
				reg	e_`depvarname'_igaussian_sq	`ldepvar'	${controls}	if	`depvarname'_sample_igaussian==1,	${option}
				est	store	m2_`depvarname'_igaussian
				
					*	Predict conditional variance
					predict	var_`depvarname'_igaussian, xb
					gen		var_`depvarname'_igaussian_abs	=	abs(var_`depvarname'_igaussian)		//	Absolute value of predicted variance (predicted value can be negative, so should use absolute value)
					lab	var	var_`depvarname'_igaussian_abs	"Conditional variance (Inv Gaussian)"
					
				*	Step 3: Construct the resilience measure
					
					*	Construct distribution parameters
					gen mu_`depvarname'_igaussian		= mean_`depvarname'_igaussian	//	mean parameter of igaussian
					gen lambda_`depvarname'_igaussian	= (mu_`depvarname'_igaussian)^3 / var_`depvarname'_igaussian_abs	//	shape parameter of igaussian (lambda). Variance is mu^3 / lambda, thus lambda = mu^3/variance

					*	Construct C&B measure
					gen `depvarname'_resil_igaussian			= 	igaussiantail(mu_`depvarname'_igaussian,lambda_`depvarname'_igauss,`threshold')	//	C&B measure (1 - Inv Guassian CDF)
					gen	`depvarname'_resil_igaussian_scale		=	`depvarname'_resil_igaussian*100	//	C&B measure scaled from 0 to 100
								
					lab	var	`depvarname'_resil_igaussian			"`depvarname' Resilience (igaussian)"
					lab	var	`depvarname'_resil_igaussian_scale		"`depvarname' Resilience (igaussian scaled)"
					
					*	Regress C&B on controls
					reg	`depvarname'_resil_igaussian	${controls},	${option}
					est	store	m3_`depvarname'_igaussian
			
			end
			
			
			
			*	It is always good to check summary stats and distribution.
			summ		${depvar_aexp},d
			kdensity	${depvar_aexp}
			
			*	Construct C&B measures of all exp, using three different assumptions
		
			resil_normal	${depvar_aexp}	${ldepvar_aexp}	${depvarname_aexp}	${threshold_aexp}	//	Normal
			resil_gamma		${depvar_aexp}	${ldepvar_aexp}	${depvarname_aexp}	${threshold_aexp}	//	Gamma
			resil_igaussian	${depvar_aexp}	${ldepvar_aexp}	${depvarname_aexp}	${threshold_aexp}	//	Inv Gaussian
				
			*	Summary stats of C&B measures
			summ	${depvarname_aexp}_resil_normal	${depvarname_aexp}_resil_gamma	${depvarname_aexp}_resil_igaussian
			
			*	Mean Square Error of different C&B measures
			*	Since we already computed sqaured residual, we just need to check their means.
			summ	e_${depvarname_aexp}_normal_sq	e_${depvarname_aexp}_gamma_sq	e_${depvarname_aexp}_igaussian_sq
			
		
			*	Similarly, we can construct resilience mesuare for other outcomes
			resil_normal	${depvar_HDDS}	${ldepvar_HDDS}	${depvarname_HDDS}	${threshold_HDDS}	//	HDDS
			resil_normal	${depvar_TLU}	${ldepvar_TLU}	${depvarname_TLU}	${threshold_TLU}	//	tlu
			
			
			*	Construct multivariate resilience measures
				
				*	Construct joint CDF - bivariate normal
				*	"binormal(z1,z2,rho)" function constructs binormal cdf. We use it to construct bivariate resilience measure
				*	Chris told that the multivariate measures should be based upon "predicted" values, not "realized" vlaues.
				*	Thus I use correlation among conditional means of each outcome
			
				pwcorr		mean_aexp_normal 	 mean_HDDS_normal 	mean_TLU_normal, sig star(0.05) // predicted outcome values
				mat	pearson_outcome_corr_coef	=	r(C)
				mat	pearson_outcome_corr_sig	=	r(sig)
				
				mat	list	pearson_outcome_corr_coef
				
				scalar	corr_allexp_HDDS		=	pearson_outcome_corr_coef[2,1]
				scalar	corr_allexp_TLU			=	pearson_outcome_corr_coef[3,1]
				scalar	corr_HDDS_TLU			=	pearson_outcome_corr_coef[3,2]
				
				
				*	Joint CDF	- bivariate normal
				loc	var	jcdf_ae_HDDS	
				cap	drop	`var'
				gen `var'	=	binormal(thresh_aexp_normal,thresh_HDDS_normal,corr_allexp_HDDS)
				lab	var	`var'	"Joint CDF (all exp and HDDS)"
				
				loc	var	jcdf_ae_TLU
				cap	drop	`var'
				gen `var'	=	binormal(thresh_aexp_normal,thresh_TLU_normal,corr_allexp_TLU)
				lab	var	`var'	"Joint CDF (all exp and TLU)"
				
				loc	var	jcdf_HDDS_TLU	
				cap	drop	`var'
				gen `var'	=	binormal(thresh_HDDS_normal,thresh_TLU_normal,corr_allexp_HDDS)
				lab	var	`var'	"Joint CDF (HDDS and TLU)"
				
				
				*	Trivariate normal CDF
				*	Requires "mvnormal()" which can be run only in Mata						
					
					*	(all exp, HDDS  and TLU(IHS)) 	  - outcome correlation
					pwcorr	mean_aexp_normal mean_HDDS_normal mean_TLU_normal	//	correlation matrix of three "predicted" outcomes
					*pwcorr	lnrexpaeu_peryear	HDDS	tlu_IHS	//	correlation matrix of three "realized" outcomes (wrong)
					mat corrmat	=	r(C)
					putmata	thresh_aexp_normal thresh_HDDS_normal thresh_TLU_normal, replace
					
					mata
					
						U = (thresh_aexp_normal, thresh_HDDS_normal, thresh_TLU_normal)
						W= st_matrix("corrmat")
						W
						R	=	vech(W)'
						R
						jcdf_ae_HDDS_TLU	=	mvnormal(U, R)	//	Trivariate normal cdf
						
					end
					
					*	Import Mata matrix back to main data, which is a joint trivariate CDF
					getmata	jcdf_ae_HDDS_TLU, replace
					lab	var	jcdf_ae_HDDS_TLU	"Joint CDF (all exp, HDDS, TLU)"
			
			
			
			*	Construct multivariate resilience measures
				
				*	Intersection measure (both thresholds should be satisfied)
				*	Bivariate: F(x1>x1_bar, x2>x2_bar) = 1 - F(x1<x1_bar) - F(x2<x2_bar) + F(x1<x1_bar,x2<x2_bar)
				*	Trivariate: F(x1>x_bar, x2>x2_bar, x3>x3_bar) = 1 - F(x1<x1_bar) - F(x2<x2_bar) - F(x3<x3_bar) + F(x1<x1_bar, x2<x2_bar) + F(x1<x1_bar, x3<x3_bar) + F(x2<x2_bar, x3<x3_bar) - F(x1<x1_bar, x2<x2_bar, x3<x3_bar) (Double-check!)
				
					*	All exp and HDDS
					loc	var	resil_int_ae_HDDS
					cap	drop	`var'
					gen	`var'	=	1	-	prob_below_aexp	-	prob_below_HDDS	+	jcdf_ae_HDDS
					lab	var	`var'		"Intersection resilience (all exp and HDDS)"
					
					*	All exp and tlu(IHS)
					loc	var	resil_int_ae_tlu_IHS
					cap	drop	`var'
					gen	`var'	=	1	-	prob_below_aexp	-	prob_below_TLU	+	jcdf_ae_TLU
					lab	var	`var'		"Intersection resilience (all exp and TLU(IHS))"
					
					*	HDDS and TLU (IHS)
					loc	var	resil_int_HDDS_TLU
					cap	drop	`var'
					gen	`var'	=	1	-	prob_below_HDDS	-	prob_below_TLU	+	jcdf_HDDS_TLU
					lab	var	`var'		"Intersection resilience (HDDS and TLU(IHS))"
				
					*	All exp, HDDS and TLU(IHS) (trivariate)
					loc	var	resil_int_ae_HDDS_TLU
					cap	drop	`var'
					gen	`var'	=	1	-	prob_below_aexp	-	prob_below_HDDS	-	prob_below_TLU	///
										+	jcdf_ae_HDDS	+	jcdf_ae_TLU	+	jcdf_HDDS_TLU	///
										-	jcdf_ae_HDDS_TLU
					lab	var	`var'	"Intersection resilience (all exp \& HDDS and TLU(IHS))"
						
				*	Union measure (at least one measure should be satisfied)
					*	Bivariate: F(x1>x1_bar, x2>x2_bar) = 1 - F(x1<x1_bar,x2<x2_bar)
					*	Trivariate: F(x1>x_bar, x2>x2_bar, x3>x3_bar) = 1 -  - F(x1<x1_bar, x2<x2_bar, x3<x3_bar) (Double-check!)
					
					*	All exp and HDDS
						loc	var	resil_int_ae_HDDS
						cap	drop	`var'
						gen	`var'	=	1	-	jcdf_ae_HDDS
						lab	var	`var'		"Union resilience (all exp and HDDS)"
						
						*	All exp and tlu(IHS)
						loc	var	resil_int_ae_tlu_IHS
						cap	drop	`var'
						gen	`var'	=	1	-	jcdf_ae_TLU
						lab	var	`var'		"Union resilience (all exp and TLU(IHS))"
						
						*	HDDS and TLU (IHS)
						loc	var	resil_int_HDDS_TLU
						cap	drop	`var'
						gen	`var'	=	1	-	jcdf_HDDS_TLU
						lab	var	`var'		"Union resilience (HDDS and TLU(IHS))"
					
						*	All exp, HDDS and TLU(IHS) (trivariate)
						loc	var	resil_int_ae_HDDS_TLU
						cap	drop	`var'
						gen	`var'	=	1	-	jcdf_ae_HDDS_TLU
						lab	var	`var'	"Union resilience (all exp \& HDDS and TLU(IHS))"
				