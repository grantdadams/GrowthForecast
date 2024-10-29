#ifdef DEBUG
  #ifndef __SUNPRO_C
    #include <cfenv>
    #include <cstdlib>
  #endif
#endif
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
#ifdef DEBUG
  #include <chrono>
#endif
#include <admodel.h>
#ifdef USE_ADMB_CONTRIBS
#include <contrib.h>

#endif
#include <df1b2fun.h>

#include <adrndeff.h>

  extern "C"  {
    void ad_boundf(int i);
  }
#include <wt.htp>

  df1b2_parameters * df1b2_parameters::df1b2_parameters_ptr=0;
  model_parameters * model_parameters::model_parameters_ptr=0;
model_data::model_data(int argc,char * argv[]) : ad_comm(argc,argv)
{
  adstring tmpstring;
  tmpstring=adprogram_name + adstring(".dat");
  if (argc > 1)
  {
    int on=0;
    if ( (on=option_match(argc,argv,"-ind"))>-1)
    {
      if (on>argc-2 || argv[on+1][0] == '-')
      {
        cerr << "Invalid input data command line option"
                " -- ignored" << endl;
      }
      else
      {
        tmpstring = adstring(argv[on+1]);
      }
    }
  }
  global_datafile = new cifstream(tmpstring);
  if (!global_datafile)
  {
    cerr << "Error: Unable to allocate global_datafile in model_data constructor.";
    ad_exit(1);
  }
  if (!(*global_datafile))
  {
    delete global_datafile;
    global_datafile=NULL;
  }
  cur_yr.allocate("cur_yr");
  styr.allocate("styr");
  endyr.allocate("endyr");
  ndat.allocate("ndat");
 nscale_parm = ndat-1;
  nyrs_data.allocate(1,ndat,"nyrs_data");
 cout <<"N yrs of data: "<<nyrs_data<<endl;
  yrs_data.allocate(1,ndat,1,nyrs_data,"yrs_data");
  age_st.allocate("age_st");
  age_end.allocate("age_end");
  wt_obs.allocate(1,ndat,1,nyrs_data,age_st,age_end,"wt_obs");
  sd_obs.allocate(1,ndat,1,nyrs_data,age_st,age_end,"sd_obs");
  ages.allocate(age_st,age_end);
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
}

void model_parameters::initializationfunction(void)
{
  L1.set_initial_value(27);
  L2.set_initial_value(46);
  log_K.set_initial_value(-.13);
  log_alpha.set_initial_value(-11.);
  log_sd_coh.set_initial_value(-1.8);
  log_sd_yr.set_initial_value(-.82);
  if (global_datafile)
  {
    delete global_datafile;
    global_datafile = NULL;
  }
}

model_parameters::model_parameters(int sz,int argc,char * argv[]) : 
 model_data(argc,argv) , function_minimizer(sz)
{
  model_parameters_ptr=this;
  initializationfunction();
  L1.allocate(5,50,1,"L1");
  L2.allocate(10,110,2,"L2");
  log_alpha.allocate(-1,"log_alpha");
  log_K.allocate(3,"log_K");
  d_scale.allocate(1,nscale_parm,age_st,age_end,phase_d_scale,"d_scale");
  wt_hat.allocate(1,ndat,1,max_nyrs_data,age_st,age_end,"wt_hat");
  #ifndef NO_AD_INITIALIZE
    wt_hat.initialize();
  #endif
  wt_pre.allocate(styr,endyr,age_st,age_end,"wt_pre");
  #ifndef NO_AD_INITIALIZE
    wt_pre.initialize();
  #endif
  alpha.allocate("alpha");
  #ifndef NO_AD_INITIALIZE
  alpha.initialize();
  #endif
  mnwt.allocate(age_st,age_end,"mnwt");
  #ifndef NO_AD_INITIALIZE
    mnwt.initialize();
  #endif
  wt_inc.allocate(age_st,age_end-1,"wt_inc");
  #ifndef NO_AD_INITIALIZE
    wt_inc.initialize();
  #endif
  sigma_coh.allocate("sigma_coh");
  sigma_yr.allocate("sigma_yr");
  K.allocate("K");
  wt_hist.allocate(styr,endyr,age_st,age_end,"wt_hist");
  wt_cur.allocate(age_st,age_end,"wt_cur");
  wt_next.allocate(age_st,age_end,"wt_next");
  wt_yraf.allocate(age_st,age_end,"wt_yraf");
  log_sd_coh.allocate(-5,10,5,"log_sd_coh");
  log_sd_yr.allocate(-5,10,6,"log_sd_yr");
  coh_eff.allocate(styr,endyr,4,"coh_eff");
  yr_eff.allocate(styr,endyr,4,"yr_eff");
  prior_function_value.allocate("prior_function_value");
  likelihood_function_value.allocate("likelihood_function_value");
  nll.allocate("nll");  /* ADOBJECTIVEFUNCTION */
}
void model_parameters::userfunction(void)
{
  nll =0.0;
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
}

void model_parameters::report(const dvector& gradients)
{
 adstring ad_tmp=initial_params::get_reportfile_name();
  ofstream report((char*)(adprogram_name + ad_tmp));
  if (!report)
  {
    cerr << "error trying to open report file"  << adprogram_name << ".rep";
    return;
  }
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
}
  long int arrmblsize=0;

int main(int argc,char * argv[])
{
#ifdef DEBUG
  auto start = std::chrono::high_resolution_clock::now();
  #ifndef __SUNPRO_C
std::feclearexcept(FE_ALL_EXCEPT);
  #endif
#endif
  ad_set_new_handler();
  ad_exit=&ad_boundf;
  arrmblsize = 1000000;
  gradient_structure::set_MAX_NVAR_OFFSET(100000);
	gradient_structure::set_NUM_DEPENDENT_VARIABLES(10000);
    gradient_structure::set_NO_DERIVATIVES();
    gradient_structure::set_YES_SAVE_VARIABLES_VALUES();
      if (!arrmblsize) arrmblsize=150000;
    df1b2variable::noallocate=1;
df1b2variable::pool = new adpool();
initial_df1b2params::varsptr = new P_INITIAL_DF1B2PARAMS[1000];
{
    df1b2_parameters mp(arrmblsize,argc,argv);
    mp.iprint = defaults::iprint;

    function_minimizer::random_effects_flag=1;
    df1b2variable::noallocate=0;
    mp.preliminary_calculations();
    mp.computations(argc,argv);
}
delete [] init_df1b2variable::list;
init_df1b2variable::list = NULL;
delete [] initial_df1b2params::varsptr;
initial_df1b2params::varsptr = NULL;
delete df1b2variable::pool;
df1b2variable::pool = NULL;
#ifdef DEBUG
  std::cout << endl << argv[0] << " elapsed time is " << std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::high_resolution_clock::now() - start).count() << " microseconds." << endl;
  #ifndef __SUNPRO_C
bool failedtest = false;
if (std::fetestexcept(FE_DIVBYZERO))
  { cerr << "Error: Detected division by zero." << endl; failedtest = true; }
if (std::fetestexcept(FE_INVALID))
  { cerr << "Error: Detected invalid argument." << endl; failedtest = true; }
if (std::fetestexcept(FE_OVERFLOW))
  { cerr << "Error: Detected overflow." << endl; failedtest = true; }
if (std::fetestexcept(FE_UNDERFLOW))
  { cerr << "Error: Detected underflow." << endl; }
if (failedtest) { std::abort(); } 
  #endif
#endif
    return 0;
}

extern "C"  {
  void ad_boundf(int i)
  {
    /* so we can stop here */
    exit(i);
  }
}

void model_parameters::preliminary_calculations(void){
  #if defined(USE_ADPVM)

  admaster_slave_variable_interface(*this);

  #endif

}

model_data::~model_data()
{}

model_parameters::~model_parameters()
{}

void model_parameters::final_calcs(void){}

void model_parameters::set_runtime(void){}

#ifdef _BORLANDC_
  extern unsigned _stklen=10000U;
#endif


#ifdef __ZTC__
  extern unsigned int _stack=10000U;
#endif

void df1b2_parameters::user_function(void)
{
  nll =0.0;
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
}
   
void df1b2_pre_parameters::setup_quadprior_calcs(void) 
{ 
  df1b2_gradlist::set_no_derivatives(); 
  quadratic_prior::in_qp_calculations=1; 
}  
  
void df1b2_pre_parameters::begin_df1b2_funnel(void) 
{ 
  (*re_objective_function_value::pobjfun)=0; 
  other_separable_stuff_begin(); 
  f1b2gradlist->reset();  
  if (!quadratic_prior::in_qp_calculations) 
  { 
    df1b2_gradlist::set_yes_derivatives();  
  } 
  funnel_init_var::allocate_all();  
}  
 
void df1b2_pre_parameters::end_df1b2_funnel(void) 
{  
  lapprox->do_separable_stuff(); 
  other_separable_stuff_end(); 
  funnel_init_var::deallocate_all(); 
} 
  
void model_parameters::begin_df1b2_funnel(void) 
{ 
  if (lapprox)  
  {  
    {  
      begin_funnel_stuff();  
    }  
  }  
}  
 
void model_parameters::end_df1b2_funnel(void) 
{  
  if (lapprox)  
  {  
    end_df1b2_funnel_stuff();  
  }  
} 

void df1b2_parameters::deallocate() 
{
  L1.deallocate();
  L2.deallocate();
  log_alpha.deallocate();
  log_K.deallocate();
  d_scale.deallocate();
  wt_hat.deallocate();
  wt_pre.deallocate();
  alpha.deallocate();
  mnwt.deallocate();
  wt_inc.deallocate();
  sigma_coh.deallocate();
  sigma_yr.deallocate();
  K.deallocate();
  wt_hist.deallocate();
  wt_cur.deallocate();
  wt_next.deallocate();
  wt_yraf.deallocate();
  log_sd_coh.deallocate();
  log_sd_yr.deallocate();
  coh_eff.deallocate();
  yr_eff.deallocate();
  prior_function_value.deallocate();
  likelihood_function_value.deallocate();
  nll.deallocate();
} 
void df1b2_parameters::allocate(void) 
{
  L1.allocate(5,50,1,"L1");
  L2.allocate(10,110,2,"L2");
  log_alpha.allocate(-1,"log_alpha");
  log_K.allocate(3,"log_K");
  d_scale.allocate(1,nscale_parm,age_st,age_end,phase_d_scale,"d_scale");
  wt_hat.allocate(1,ndat,1,max_nyrs_data,age_st,age_end,"wt_hat");
  #ifndef NO_AD_INITIALIZE
    wt_hat.initialize();
  #endif
  wt_pre.allocate(styr,endyr,age_st,age_end,"wt_pre");
  #ifndef NO_AD_INITIALIZE
    wt_pre.initialize();
  #endif
  alpha.allocate("alpha");
  #ifndef NO_AD_INITIALIZE
  alpha.initialize();
  #endif
  mnwt.allocate(age_st,age_end,"mnwt");
  #ifndef NO_AD_INITIALIZE
    mnwt.initialize();
  #endif
  wt_inc.allocate(age_st,age_end-1,"wt_inc");
  #ifndef NO_AD_INITIALIZE
    wt_inc.initialize();
  #endif
  sigma_coh.allocate("sigma_coh");
  sigma_yr.allocate("sigma_yr");
  K.allocate("K");
  wt_hist.allocate(styr,endyr,age_st,age_end,"wt_hist");
  wt_cur.allocate(age_st,age_end,"wt_cur");
  wt_next.allocate(age_st,age_end,"wt_next");
  wt_yraf.allocate(age_st,age_end,"wt_yraf");
  log_sd_coh.allocate(-5,10,5,"log_sd_coh");
  log_sd_yr.allocate(-5,10,6,"log_sd_yr");
  coh_eff.allocate(styr,endyr,4,"coh_eff");
  yr_eff.allocate(styr,endyr,4,"yr_eff");
  prior_function_value.allocate("prior_function_value");
  likelihood_function_value.allocate("likelihood_function_value");
  nll.allocate("nll");  /* ADOBJECTIVEFUNCTION */
}
