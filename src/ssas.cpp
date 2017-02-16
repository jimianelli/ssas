#ifdef DEBUG
  #ifndef __SUNPRO_C
    #include <cfenv>
    #include <cstdlib>
  #endif
#endif
  #include "admodel.h"
  const double pi = 3.141592654;
  ofstream mcmc_Fbar("mcmc_Fbar.csv");
  ofstream mcmc_hr("mcmc_hr.csv");
  ofstream mcmc_ssb("mcmc_ssb.csv");
  ofstream mcmc_B4("mcmc_B4.csv");
  ofstream mcmc_rec("mcmc_rec.csv");
  ofstream mcmc_yieldPred("mcmc_yieldPred.csv");
#include <admodel.h>
#include <contrib.h>

  extern "C"  {
    void ad_boundf(int i);
  }
#include <gdbprintlib.cpp>

#include <ssas.htp>

model_data::model_data(int argc,char * argv[]) : ad_comm(argc,argv)
{
  lambda.allocate(1,10,"lambda");
  nyrs_r.allocate("nyrs_r");
  datafile.allocate("datafile");
  R_output.allocate("R_output");
  nsel_coff.allocate("nsel_coff");
  nsel_changes.allocate("nsel_changes");
  yrs_sel_changes.allocate(1,nsel_changes,"yrs_sel_changes");
 ad_comm::change_datafile_name(datafile); 
 cout << "+----------------------+" << endl;
 cout << "| Reading data file    |" << endl;
 cout << "+----------------------+" << endl;
  nyrs.allocate("nyrs");
  nages.allocate("nages");
  nyrs_surv1.allocate("nyrs_surv1");
  nages_surv1.allocate("nages_surv1");
  nyrs_surv2.allocate("nyrs_surv2");
  nages_surv2.allocate("nages_surv2");
  obs_catch_at_age.allocate(1,nyrs,1,nages,"obs_catch_at_age");
  obs_surv1_at_age.allocate(1,nyrs_surv1,1,nages_surv1,"obs_surv1_at_age");
  obs_surv2_at_age.allocate(1,nyrs_surv2,1,nages_surv2,"obs_surv2_at_age");
  M.allocate("M");
  weight_at_age.allocate(1,nyrs,1,nages,"weight_at_age");
  mat_at_age.allocate(1,nyrs,1,nages,"mat_at_age");
  denom_res_mat.allocate(1,nyrs,1,nages,"denom_res_mat");
  denom_res_mat_surv1.allocate(1,nyrs_surv1,1,nages_surv1,"denom_res_mat_surv1");
  denom_res_mat_surv2.allocate(1,nyrs_surv2,1,nages_surv2,"denom_res_mat_surv2");
  ages.allocate(1,nages);
  ages4plus.allocate(1,nages-1);
  years.allocate(1,nyrs);
  yield_obs.allocate(1,nyrs);
  nyrs = nyrs - nyrs_r;	
  nyrs_surv1 = nyrs_surv1 - nyrs_r;
  nyrs_surv2 = nyrs_surv2 - nyrs_r;		
  cout<<nsel_changes<<endl;
  if (nsel_changes>0)
  {
    int itmp=1;
    for (int i=1;i<=nyrs; i++)
    {
      if (yrs_sel_changes(itmp)==i)
        itmp++;
      if (itmp>nsel_changes)
        itmp--;
      if (yrs_sel_changes(itmp)>=nyrs)
        nsel_changes = itmp-1;
    }
  }
  cout<<nsel_changes<<endl;
}

model_parameters::model_parameters(int sz,int argc,char * argv[]) : 
 model_data(argc,argv) , function_minimizer(sz)
{
  initializationfunction();
  log_sel_coff.allocate(1,nsel_coff,2,"log_sel_coff");
  log_sel_dev.allocate(1,nsel_changes,1,nsel_coff,-5,5,3,"log_sel_dev");
  log_mean_F.allocate(1,"log_mean_F");
  log_fy_devs.allocate(1,nyrs,-10,10,2,"log_fy_devs");
  log_mean_recruit.allocate(1,"log_mean_recruit");
  log_recruit_devs.allocate(1,nyrs,-50,50,2,"log_recruit_devs");
  log_initpop_devs.allocate(1,nages-1,2,"log_initpop_devs");
  log_q1.allocate(1,nages_surv1,-15,15,2,"log_q1");
  log_q2.allocate(1,nages_surv2,-15,15,2,"log_q2");
  log_alfa.allocate(0,1,3,"log_alfa");
  log_beta.allocate(-5,0,3,"log_beta");
  log_sel.allocate(1,nyrs,1,nages,"log_sel");
  #ifndef NO_AD_INITIALIZE
    log_sel.initialize();
  #endif
  avgsel_fsh.allocate(1,nyrs,"avgsel_fsh");
  #ifndef NO_AD_INITIALIZE
    avgsel_fsh.initialize();
  #endif
  log_fy.allocate(1,nyrs,"log_fy");
  #ifndef NO_AD_INITIALIZE
    log_fy.initialize();
  #endif
  F.allocate(1,nyrs,1,nages,"F");
  #ifndef NO_AD_INITIALIZE
    F.initialize();
  #endif
  Z.allocate(1,nyrs,1,nages,"Z");
  #ifndef NO_AD_INITIALIZE
    Z.initialize();
  #endif
  S.allocate(1,nyrs,1,nages,"S");
  #ifndef NO_AD_INITIALIZE
    S.initialize();
  #endif
  N.allocate(1,nyrs,1,nages,"N");
  #ifndef NO_AD_INITIALIZE
    N.initialize();
  #endif
  C.allocate(1,nyrs,1,nages,"C");
  #ifndef NO_AD_INITIALIZE
    C.initialize();
  #endif
  U1.allocate(1,nyrs_surv1,1,nages_surv1,"U1");
  #ifndef NO_AD_INITIALIZE
    U1.initialize();
  #endif
  U2.allocate(1,nyrs_surv2,1,nages_surv2,"U2");
  #ifndef NO_AD_INITIALIZE
    U2.initialize();
  #endif
  neg_log_like.allocate(1,10,"neg_log_like");
  #ifndef NO_AD_INITIALIZE
    neg_log_like.initialize();
  #endif
  ObjFun.allocate("ObjFun");
  prior_function_value.allocate("prior_function_value");
  likelihood_function_value.allocate("likelihood_function_value");
  recsum.allocate("recsum");
  #ifndef NO_AD_INITIALIZE
  recsum.initialize();
  #endif
  initsum.allocate("initsum");
  #ifndef NO_AD_INITIALIZE
  initsum.initialize();
  #endif
  B4.allocate(1,nyrs,"B4");
  #ifndef NO_AD_INITIALIZE
    B4.initialize();
  #endif
  hr.allocate(1,nyrs,"hr");
  #ifndef NO_AD_INITIALIZE
    hr.initialize();
  #endif
  avg_F.allocate("avg_F");
  #ifndef NO_AD_INITIALIZE
  avg_F.initialize();
  #endif
  yield_pred.allocate(1,nyrs,"yield_pred");
  #ifndef NO_AD_INITIALIZE
    yield_pred.initialize();
  #endif
  Fbar.allocate(1,nyrs,"Fbar");
  ssb.allocate(1,nyrs,"ssb");
  ssb_m3.allocate(1,nyrs-3,"ssb_m3");
  rec.allocate(1,nyrs,"rec");
  rec_m3.allocate(1,nyrs-3,"rec_m3");
  pred_rec_m3.allocate(1,nyrs-3,"pred_rec_m3");
  #ifndef NO_AD_INITIALIZE
    pred_rec_m3.initialize();
  #endif
}

void model_parameters::preliminary_calculations(void)
{

#if defined(USE_ADPVM)

  admaster_slave_variable_interface(*this);

#endif
  ages.fill_seqadd(3,1);           // vector of ages
  ages4plus.fill_seqadd(4,1);      // vector of non recruiting ages
  years.fill_seqadd(1961,1);       // fill vector of years with years
  pred_year=years[nyrs]+1;         // year of prediction = last_year +1
  yield_obs=rowsum(elem_prod(obs_catch_at_age,weight_at_age));
}

void model_parameters::initializationfunction(void)
{
  log_mean_recruit.set_initial_value(10.);
  log_mean_F.set_initial_value(-1.2);
}

void model_parameters::userfunction(void)
{
  ObjFun =0.0;
  get_mortality_and_survivial_rates();
  get_numbers_at_age();
  get_catch_at_age();
  get_survey1_at_age();
  get_survey2_at_age();
  ssb_rec_ricker();
  evaluate_the_objective_function();
}

void model_parameters::get_mortality_and_survivial_rates(void)
{
  int i, j, iselyr;
  iselyr = 1;
   //calculate the selectivity from the sel_coffs ---------------------
  for (i=1;i<=nyrs;i++)
  { 
    if (i==1)
      log_sel(1)(1,nsel_coff) = log_sel_coff;
    else
      log_sel(i)(1,nsel_coff) = log_sel(i-1)(1,nsel_coff) ;
    if( nsel_changes>0)
    {
      if(i== yrs_sel_changes(iselyr) )
      {
        for (j=1;j<=nsel_coff;j++)
          log_sel(i,j) += log_sel_dev(iselyr,j);
        if (iselyr<nsel_changes)
          iselyr++;
      }
    }
    log_sel(i)(nsel_coff+1,nages) = log_sel(i,nsel_coff);
    avgsel_fsh(i) = log(mean(mfexp(log_sel(i))));
    // Condition selectivities to have mean of zero
    log_sel(i) -= mean(log_sel(i));
    log_fy(i) = log_fy_devs(i)  + log_mean_F;
    F(i)      = mfexp(log_fy(i) + log_sel(i));  
  }
  // get the total mortality
  Z = F+M;
  // get the survival rate
  S = mfexp(-Z);
}

void model_parameters::get_numbers_at_age(void)
{
  int i, j;
  for (i=1;i<=nyrs;i++)
  {
    //Number of recruits each year
    //Modeled as deviations around a median recruitment
    N(i,1) = mfexp(log_mean_recruit+log_recruit_devs(i));
  }
  for (j=1;j<=nages-1;j++)
  {
    //Numbers at age in initial year
    //Modeled as differences
    N(1,j+1) = exp(log_mean_recruit -double(j)*M + log_initpop_devs(j));
  }
 // calculate remaining numbers-at-age
  for (i=1;i<nyrs;i++)                   //  for (i=1;i<nyrs;i++)         
  {                                      //  {                            
    for (j=1;j<nages;j++)              //    for (j=1;j<nages-1;j++) // if plus-group is calculated differently     
    {                                    //    {                          
      N(i+1,j+1)=N(i,j)*S(i,j);          //      N(i+1,j+1)=N(i,j)*S(i,j);
    }//Calculate numbers in +group (14+) //    }                          
    N(i+1,nages) += N(i,nages)*S(i,nages);
  }//-----------------------------------------------------------------------------                                                                 
}

void model_parameters::get_catch_at_age(void)
{
  int i;
  C=elem_prod(elem_div(F,Z),elem_prod(1.-S,N));
  for (i=1;i<=nyrs;i++)
  {                                               // get B4+ Biomass
    yield_pred(i) = C(i)*weight_at_age(i);
    // B4(i)=sum(elem_prod(N(i).sub(2,nages).shift(1),weight_at_age(i).sub(2,nages).shift(1)));
		B4(i) = N(i)(2,nages) * weight_at_age(i)(2,nages);
    hr(i) = yield_obs(i)/B4(i);
  }
}

void model_parameters::get_survey1_at_age(void)
{
  int i, j;
  for (j=1;j<=nages_surv1;j++)
  {
    for (i=1;i<=nyrs_surv1;i++)
    {
    //if(i<=11)
    U1(i,j)=exp(log_q1(j) + log(N(i,j)));
    //else 
    //U1(i,j)=exp(log_q1_2(j) + log(N(i,j)));
    }
  }
  //-----------------------------------------------------------------------------                                                                  
}

void model_parameters::get_survey2_at_age(void)
{
  int i, j;
  for (j=1;j<=nages_surv2;j++)
  {
    for (i=1;i<=nyrs_surv2;i++)
    {
    //if(i<=11)
    U2(i,j)=exp(log_q2(j) + log(N(i,j)));
    //else 
    //U2(i,j)=exp(log_q2_2(j) + log(N(i,j)));
    }
  }
  //-----------------------------------------------------------------------------                                                                  
}

void model_parameters::ssb_rec_ricker(void)
{
  ssb=rowsum(elem_prod(elem_prod(N,weight_at_age),mat_at_age))/1000;
  ssb_m3=ssb.sub(1,nyrs-3);
  rec=column(N,1)/1000;
  rec_m3=rec.sub(4,nyrs).shift(1);
  for (int i=1;i<=nyrs-3;i++)
  {
    pred_rec_m3(i)=exp(log_alfa)*ssb_m3(i)*exp(-exp(log_beta)*ssb_m3(i));
  }  
  //-----------------------------------------------------------------------------                                                                  
}

void model_parameters::evaluate_the_objective_function(void)
{
  avg_F=sum(F)/double(size_count(F));
  if (last_phase())
    ObjFun += .001*square(log(avg_F/.2));                    
  else                                                 
    ObjFun += 1000.*square(log(avg_F/.2));                   
  ObjFun += 20. * norm2(avgsel_fsh);
  Concentrated_Likelihoods();
  ObjFun  += sum(neg_log_like);
  if(last_phase())
   for(int i=1;i<=nyrs;i++)
      Fbar(i) = mean(F(i)(2,6));
  if(mceval_phase())
    write_mcmc();
}

void model_parameters::Concentrated_Likelihoods(void)
{
  neg_log_like.initialize();
  // Condition model to fit total catch
  neg_log_like(1) = lambda(1) * norm2(log(yield_pred(1,nyrs)+.1)-log(yield_obs(1,nyrs)+.1));
  // Fit catch age
  neg_log_like(2) = lambda(2) * 0.5*double(size_count(C)) *
                     log( sum(elem_div(square(C-obs_catch_at_age.sub(1,nyrs) ),.1+elem_prod(denom_res_mat.sub(1,nyrs),C))))          ;
  // Fit to survey data (not sure what's in denominator here...) lognormal
  neg_log_like(3) = lambda(3) * log( sum(elem_div(square(U1-obs_surv1_at_age.sub(1,nyrs_surv1)),.1+elem_prod(denom_res_mat_surv1.sub(1,nyrs_surv1),U1)))) ;
  neg_log_like(4) = lambda(4) * log( sum(elem_div(square(U2-obs_surv2_at_age.sub(1,nyrs_surv2)),.1+elem_prod(denom_res_mat_surv2.sub(1,nyrs_surv2),U2)))) ;
  // SRR likelihood? ...concentrated so variance term falls out (it is implicitly estimated)
  neg_log_like(5) = lambda(5) * (nyrs-3)*0.5*log( sum(square(rec_m3.sub(1,nyrs-3)-pred_rec_m3.sub(1,nyrs-3))));                                                                      
  // Penalty smoother for selectivities
 if (active(log_sel_dev))
  {
    // Minor penalty to keep centered on zero
    neg_log_like(6) = lambda(6) * norm2(log_sel_dev);
    // Penalty for changes over time
    for (int j=1;j<=nsel_coff;j++)  
      neg_log_like(7) += lambda(7) * norm2(trans(log_sel_dev)(j));
    // Penalty for changes over ages
    for (int i=1;i<=nsel_changes;i++)  
      neg_log_like(8) += lambda(8) * norm2(first_difference(first_difference(log_sel(yrs_sel_changes(i)))));
  }
  // impose a minor penalty deviations
  neg_log_like(9)  = lambda(9)  * norm2(log_recruit_devs);
  neg_log_like(10) = lambda(10) * norm2(log_initpop_devs);
}

void model_parameters::write_mcmc(void)
{
  // get Fbar
  for(int i=1;i<=nyrs;i++)
    Fbar(i)=mean(F(i)(2,6));
  mcmc_Fbar       << Fbar       << endl;
  mcmc_hr         << hr         << endl;
  mcmc_ssb        << ssb        << endl;
  mcmc_B4         << B4         << endl;
  mcmc_rec        << rec        << endl;
  mcmc_yieldPred  << yield_pred << endl;
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
  cout   << "+----------------------+" << endl;
  cout   << "| Done phase "<<current_phase()<<" "<<endl;
  cout   << "+----------------------+" << endl;
  report << "ssb =  " << endl;
  report << ssb <<endl;
  report << "ssb_m3 =  " << endl;
  report << ssb_m3 <<endl;
  report << "rec_m3 =  " << endl;
  report << rec_m3 <<endl;
  report << "rec =  " << endl;
  report << rec <<endl;
  report << neg_log_like <<endl;
  if (last_phase())
    write_out();
  // OUTPUT DATA OUTPUT THROUGH AN "output file"
  // rather than using the "report file" both are comparable
}

void model_parameters::write_out(void)
{
  ofstream out(R_output); //   ofstream out("sai6109_rm_pt_sr_3sel_9609.out");
  out.precision(8);
  // FORMATTING SHOULD BE AS FOLLOWS ,, from "/media/Data/ADMB/NCEAS ADMB Course files/presentations/7b Formatting output.ppt"
  // out<<"NAME (with no spaces)"<<endl;
  // out<<NROWS<<" "<<NCOLS<<endl;
  // out<<VARIABLE<<endl;
  out<<"N"<<endl;
  out<<nyrs<<" "<<nages<<endl;
  out<<N<<endl;
  out<<"ricker_log_alfa"<<endl;
  out<<1<<" "<<1<<endl;
  out<<log_alfa<<endl;
  out<<"ricker_log_beta"<<endl;
  out<<1<<" "<<1<<endl;
  out<<log_beta<<endl;
  out<<"F"<<endl;
  out<<nyrs<<" "<<nages<<endl;
  out<<F<<endl;
  out<<"Chat"<<endl;
  out<<nyrs<<" "<<nages<<endl;
  out<<C<<endl;
  out<<"Cobs"<<endl;
  out<<nyrs<<" "<<nages<<endl;
  out<<obs_catch_at_age.sub(1,nyrs)<<endl;
  out<<"U1obs"<<endl;
  out<<nyrs_surv1<<" "<<nages_surv1<<endl;
  out<<obs_surv1_at_age.sub(1,nyrs_surv1)<<endl;
  out<<"U1hat"<<endl;
  out<<nyrs_surv1<<" "<<nages_surv1<<endl;
  out<<U1<<endl;
  out<<"U2obs"<<endl;
  out<<nyrs_surv2<<" "<<nages_surv2<<endl;
  out<<obs_surv2_at_age.sub(1,nyrs_surv2)<<endl;
  out<<"U2hat"<<endl;
  out<<nyrs_surv2<<" "<<nages_surv2<<endl;
  out<<U2<<endl;
  out<<"W"<<endl;
  out<<nyrs<<" "<<nages<<endl;
  out<<weight_at_age.sub(1,nyrs)<<endl;
  out<<"Mat"<<endl;
  out<<nyrs<<" "<<nages<<endl;
  out<<mat_at_age.sub(1,nyrs)<<endl;
  out<<"log_q1"<<endl;
  out<<1<<" "<<nages_surv1<<endl;
  out<<log_q1<<endl;
  out<<"log_q2"<<endl;
  out<<1<<" "<<nages_surv2<<endl;
  out<<log_q2<<endl;
  out<<"denom_res_mat"<<endl;
  out<<1<<" "<<nages<<endl;
  out<<denom_res_mat[1]<<endl;
  out<<"denom_res_mat_surv1"<<endl;
  out<<1<<" "<<nages_surv1<<endl;
  out<<denom_res_mat_surv1[1]<<endl;
  out<<"loglikeComp"<<endl;
  out<<1<<" "<<10   <<endl;
  out<<neg_log_like <<endl;
  out<<"Catch_pred" <<endl;
  out<<1<<" "<<nyrs <<endl;
  out<<yield_pred(1,nyrs)  <<endl;
  out<<"Catch_obs"  <<endl;
  out<<1<<" "<<nyrs <<endl;
  out<<yield_obs(1,nyrs)    <<endl;
  out<<"rec"        <<endl;
  out<<1<<" "<<nyrs <<endl;
  out<<rec          <<endl;
  out<<"SSB"        <<endl;
  out<<1<<" "<<nyrs <<endl;
  out<<ssb          <<endl;
  out<<"sel"        <<endl;
  out<<nages<<" "<<nyrs <<endl;
  out<<exp(log_sel) <<endl;
  out<<"Fbar"        <<endl;
  out<<1 <<" "<<nyrs <<endl;
  out<< Fbar        <<endl;
  out<<"nage_s1"   <<endl;
  out<<1 <<" "<<1  <<endl;
  out<< nages_surv1 <<endl;
  out<<"nage_s2"   <<endl;
  out<<1 <<" "<<1  <<endl;
  out<< nages_surv2 <<endl;
  out<<"nyrs_s1"   <<endl;
  out<<1 <<" "<<1  <<endl;
  out<< nyrs_surv1 <<endl;
  out<<"nyrs_s2"   <<endl;
  out<<1 <<" "<<1  <<endl;
  out<< nyrs_surv2 <<endl;
  out.close();
}

void model_parameters::set_runtime(void)
{
  dvector temp("{1.e-1,1.e-08}");
  convergence_criteria.allocate(temp.indexmin(),temp.indexmax());
  convergence_criteria=temp;
  dvector temp1("{100,100,2500}");
  maximum_function_evaluations.allocate(temp1.indexmin(),temp1.indexmax());
  maximum_function_evaluations=temp1;
}

void model_parameters::final_calcs()
{
  ofstream out(R_output,ios::app); 
  out<<"SSB.sd"        <<endl;
  out<<1 <<" "<<nyrs <<endl;
  out<< ssb.sd        <<endl;
  out<<"rec.sd"        <<endl;
  out<<1 <<" "<<nyrs <<endl;
  out<< rec.sd        <<endl;
  out<<"Fbar.sd"        <<endl;
  out<<1 <<" "<<nyrs <<endl;
  out<< Fbar.sd        <<endl;
  out.close();
}

model_data::~model_data()
{}

model_parameters::~model_parameters()
{}

#ifdef _BORLANDC_
  extern unsigned _stklen=10000U;
#endif


#ifdef __ZTC__
  extern unsigned int _stack=10000U;
#endif

  long int arrmblsize=0;

int main(int argc,char * argv[])
{
    ad_set_new_handler();
  ad_exit=&ad_boundf;
  gradient_structure::set_NUM_DEPENDENT_VARIABLES(500);
    gradient_structure::set_NO_DERIVATIVES();
#ifdef DEBUG
  #ifndef __SUNPRO_C
std::feclearexcept(FE_ALL_EXCEPT);
  #endif
#endif
    gradient_structure::set_YES_SAVE_VARIABLES_VALUES();
    if (!arrmblsize) arrmblsize=15000000;
    model_parameters mp(arrmblsize,argc,argv);
    mp.iprint=10;
    mp.preliminary_calculations();
    mp.computations(argc,argv);
#ifdef DEBUG
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
