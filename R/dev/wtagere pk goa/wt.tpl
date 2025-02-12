//----------------------------
// Goal to predict fishery weights at age
// Data include survey and fishery data
//----------------------------
//----------------------------

DATA_SECTION
 // Retrospective years to test over (0 if current year) used to help w/ scoring which method
  init_int cur_yr ;// Retrospective years to test over (0 if current year)
  int retro;
 // Time frame for whole model (doesn't have to match data period)
  init_int styr;
  init_int endyr;
	// Number of datasets
  init_int ndat;
	// Scalar for surveys (always expects fishery wts as first dataset)...
  int nscale_parm;
  !! nscale_parm = ndat-1;
	// Number of years of observations in each data set
	init_ivector nyrs_data(1,ndat);
	!! cout <<"N yrs of data: "<<nyrs_data<<endl;
	// Actual years of observations in each data set (ragged array)
  init_imatrix yrs_data(1,ndat,1,nyrs_data);
	// First age and last, same for all datasets...
  init_int age_st;
  init_int age_end;
  int nages;
	// Data has all observations next to eachother then all sd's next to eachother
  init_3darray wt_obs(1,ndat,1,nyrs_data,age_st,age_end);
  init_3darray sd_obs(1,ndat,1,nyrs_data,age_st,age_end);
  int phase_d_scale;
  vector ages(age_st,age_end);
  // Kludge because ragged not working correctly...
  int max_nyrs_data;
 LOC_CALCS
   log_input(cur_yr);
   log_input(styr);
   log_input(endyr);
   log_input(ndat);
   log_input(nyrs_data);
   log_input(yrs_data);
   log_input(age_st);
   log_input(age_end);
   log_input(wt_obs);
   log_input(sd_obs);
   max_nyrs_data = max(nyrs_data);
  // Need to reset data-range for the years covered (for retrospective fitting)
  retro = endyr-cur_yr-2;
  // NOTE This may not reduce sparse year data appropriately (assumed to be annual)
  for (int i=1;i<=ndat;i++) 
    nyrs_data(i) -= retro;
  // here
  if (ndat>1) phase_d_scale = 1; else phase_d_scale = -1;
  nages = age_end - age_st + 1;
  for (int i=age_st;i<=age_end;i++) ages(i) = double(i);
  // COUT(sd_obs); 
  cout<< ndat <<" "<<styr <<" "<< endyr <<" "<< age_st <<" "<< age_end<<endl<<yrs_data<<endl; //exit(1);
 END_CALCS
  
INITIALIZATION_SECTION
  L1 27
  L2 46
  log_K -.13
  log_alpha -11.
  // log_t0 -3.4
  log_sd_coh -1.8
  log_sd_yr -.82

PARAMETER_SECTION
  init_bounded_number L1(5,50,1);
  init_bounded_number L2(10,110,2);
  init_number log_alpha(-1);
  init_number log_K(3);
  init_matrix d_scale(1,nscale_parm,age_st,age_end,phase_d_scale);
  // init_number log_t0(3);
  // Predicted weight matrix
  3darray wt_hat(1,ndat,1,max_nyrs_data,age_st,age_end);
  // 3darray wt_hat(1,ndat,1,nyrs_data,age_st,age_end);
  matrix  wt_pre(styr,endyr,age_st,age_end);
  number alpha;
  vector mnwt(age_st,age_end);
  vector wt_inc(age_st,age_end-1);
  sdreport_number sigma_coh;
  sdreport_number sigma_yr;
  sdreport_number K;
  // Current, next, and year-after next predictions here
  sdreport_matrix wt_hist(styr,endyr,age_st,age_end);
  sdreport_vector wt_cur(age_st,age_end);
  sdreport_vector wt_next(age_st,age_end);
  sdreport_vector wt_yraf(age_st,age_end);

  init_bounded_number log_sd_coh(-5,10,5);
  init_bounded_number log_sd_yr(-5,10,6);
  random_effects_vector coh_eff(styr,endyr,4);
  random_effects_vector yr_eff(styr,endyr,4);
  // Uncomment this (and comment above) for fixed effects estimation (boo)
  //init_vector coh_eff(styr,endyr,3);
  //init_vector yr_eff(styr,endyr,4);

  objective_function_value nll;

PROCEDURE_SECTION
  // Cohort and year effects variance parameters
  sigma_coh    = mfexp(log_sd_coh);
  sigma_yr     = mfexp(log_sd_yr );
  K            = mfexp(log_K);
  alpha        = mfexp(log_alpha);
  for (int j=age_st;j<=age_end;j++)
  {
    mnwt(j)    = alpha * pow(L1 + (L2-L1)*(1.-pow(K,double(j-age_st))) / (1.-pow(K,double(nages-1))) ,3);
  }
  wt_inc       = --mnwt(age_st+1,age_end) - mnwt(age_st,age_end-1);

  // Initialize first year
  wt_pre(styr)    = mnwt;

  // subsequent years
  for (int i=styr+1;i<=endyr;i++)
  {
    wt_pre(i,age_st) = mnwt(age_st)*mfexp(square(sigma_coh)/2.+sigma_coh*coh_eff(i));
    if (last_phase())
      wt_pre(i)(age_st+1,age_end) = ++(wt_pre(i-1)(age_st,age_end-1) + wt_inc*mfexp(square(sigma_yr)/2. + sigma_yr*yr_eff(i)));
    else
      wt_pre(i)(age_st+1,age_end) = ++(wt_pre(i-1)(age_st,age_end-1) + wt_inc*mfexp(sigma_yr*yr_eff(i)));
  }
  int iyr;

  // Fit global mean to all years...
  for (int h = 1;h<=ndat;h++) // Loop over number of data sets (assuming 1st is target)
  {
    for (int i=1;i<=nyrs_data(h);i++) // Loop over how many observations w/in a data set
    {
      iyr = yrs_data(h,i);
    	if (h>1) // First data set is correct scale, latter sets need a multiplier (d_scale)
    		wt_hat(h,i) = elem_prod(d_scale(h-1) , wt_pre(iyr) );
    	else
        wt_hat(h,i) = wt_pre(iyr);

      for (int j=age_st;j<=age_end;j++)
      {
				// cout<<nll<<endl;
        nll += square(wt_obs(h,i,j) - mnwt(j))      /(2.*square(sd_obs(h,i,j)));
        nll += square(wt_obs(h,i,j) - wt_hat(h,i,j))/(2.*square(sd_obs(h,i,j)));
      }
    }
  }
  // Random effects are on N(0,1) scale
  nll += 0.5*norm2(coh_eff);
  nll += 0.5*norm2( yr_eff);

  if (sd_phase())
  {
    wt_cur  = wt_pre(cur_yr);  
    wt_next = wt_pre(endyr-1); 
    wt_yraf = wt_pre(endyr);  
		wt_hist = wt_pre;
  }

REPORT_SECTION
  if (last_phase())
  {
  int retyr ;
  retyr = int (endyr - endyr);
  switch ( retyr ) 
  {
    case 0 : 
      report << "lof1"<<endl; report << 0 <<endl;
      report << "lof2"<<endl; report << 0 <<endl;
      report << "lof3"<<endl; report << 0 <<endl;
    break;

    case 1 : 
      report << "lof1"<<endl; report << pow(norm2(wt_obs(1,endyr-3)-wt_cur),.5) <<endl;
      report << "lof2"<<endl; report << 0 <<endl;
      report << "lof3"<<endl; report << 0 <<endl;
    break;

    case 2 : 
      report << "lof1"<<endl; report << pow(norm2(wt_obs(1,endyr-3)-wt_cur),.5) <<endl;
      report << "lof2"<<endl; report << pow(norm2(wt_obs(1,endyr-2)-wt_next),.5) <<endl;
      report << "lof3"<<endl; report << 0 <<endl;
    break;

    default : 
      report << "lof1"<<endl; report << pow(norm2(wt_obs(1,endyr-3)-wt_cur),.5) <<endl;
      report << "lof2"<<endl; report << pow(norm2(wt_obs(1,endyr-2)-wt_next),.5) <<endl;
      report << "lof3"<<endl; report << pow(norm2(wt_obs(1,endyr-1)-wt_yraf),.5) <<endl;
    break;
  }
  }
    /* */ 
  report <<"endyr" << endl;
  report << endyr << endl;
  REPORT(retro);
  REPORT(cur_yr);
  report <<"data"   << endl;
  report <<wt_obs   << endl;
  report << "yr"<<endl; 
  for (int iyr=styr;iyr<=endyr;iyr++) report <<iyr<<" "<<endl;
  report <<"wt_pre"<<endl;
  report << wt_pre<<endl;
  if (last_phase())
  {
    for (int h = 1;h<=ndat;h++)
    {
      dmatrix resid(1,nyrs_data(h),age_st,age_end);
      for (int i=1;i<=nyrs_data(h);i++)
        resid(i) = value(wt_obs(h,i)-wt_hat(h,i));
      report <<"residuals_"<<h<<endl;
      report << resid << endl;
    }
  }
  REPORT(sigma_yr);
  REPORT(yr_eff);
  REPORT(sigma_coh);
  report <<"cohort"<<endl;
  // for (int j=age_st;j>=age_end;j--)
    // report << styr-j <<" ";
  for (int i=styr+1;i<=endyr;i++)
    report << i-age_st <<" ";
  report<<endl;
  report <<"coh_eff"<<endl;
  // for (int j=age_end;j>=age_st;j--)
    // report << sigma_coh*coh_eff(styr-j) <<" ";
  for (int i=styr;i<=endyr;i++)
    report << sigma_coh*coh_eff(i) <<" ";
  report<<endl;
  REPORT(ages);
  REPORT(mnwt);
  // REPORT(winf);
  REPORT(K);
  REPORT(L1);
  REPORT(L2);
  // REPORT(t0);


GLOBALS_SECTION
  #include <admodel.h>  
  #undef REPORT
  #define REPORT(object) report << #object "\n" << setw(8)  << setprecision(6) << setfixed() << object << endl;

  #undef ReadLog
  #define ReadLog(object) Log << #object "\n" << setw(8)  << setprecision(6) << setfixed() << object << endl;

	#undef log_input
	#define log_input(object) write_input_log << "# " #object "\n" << object << endl;
  ofstream write_input_log("input.log");

  #undef COUT
  #define COUT(object) cout << #object "\n" << setw(6) << setprecision(5) << setfixed() << object << endl;


TOP_OF_MAIN_SECTION
  arrmblsize = 1000000;
  gradient_structure::set_MAX_NVAR_OFFSET(100000);
	gradient_structure::set_NUM_DEPENDENT_VARIABLES(10000);
