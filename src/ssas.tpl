DATA_SECTION
  init_vector lambda(1,10);     // number of retrospective years  
  init_int nyrs_r;				 			// number of retrospective years	
	init_adstring datafile;
	init_adstring R_output;
  init_int nsel_coff;          // Number of selectivity coeficients to estimate (i.e., last value set as plus)
  init_int nsel_changes;
  init_ivector yrs_sel_changes(1,nsel_changes);
	!! ad_comm::change_datafile_name(datafile); 
	!! cout << "+----------------------+" << endl;
	!! cout << "| Reading data file    |" << endl;
	!! cout << "+----------------------+" << endl;

  init_int nyrs;                                 			// the number of years in catch data		c@a
  init_int nages;                                			// the number of age classess in catch data	c@a
  init_int nyrs_surv1;                            			// the number of years in survey 1		u@a
  init_int nages_surv1;                           			// the number of age classess in survey 1	u@a
  init_int nyrs_surv2;                            			// the number of years in survey 2		u@a
  init_int nages_surv2;                           			// the number of age classess in survey 2	u@a
  init_matrix obs_catch_at_age(1,nyrs,1,nages);  			// observed catch-at-age data			c@a
  init_matrix obs_surv1_at_age(1,nyrs_surv1,1,nages_surv1);	        // observed survey1-at-age data			u@a
  init_matrix obs_surv2_at_age(1,nyrs_surv2,1,nages_surv2);	        // observed survey2-at-age data			u@a
  init_number M;                                 			// estimate of natural mortality rate
  //init_vector relwt(2,nages);                    			// need to have relative weight-at-age to calculate B2+
  init_matrix weight_at_age(1,nyrs,1,nages);     			// observed weigth-at-age data
  init_matrix mat_at_age(1,nyrs,1,nages);        			// maturity-at-age data
  init_matrix denom_res_mat(1,nyrs,1,nages);     			// age-based denominator for catch residuals (Cobs/Chat)^2/denom_res
  init_matrix denom_res_mat_surv1(1,nyrs_surv1,1,nages_surv1);	        // age-based denominator for survey1 residuals (Uobs/Uhat)^2/denom_res
  init_matrix denom_res_mat_surv2(1,nyrs_surv2,1,nages_surv2);	        // age-based denominator for survey2 residuals (Uobs/Uhat)^2/denom_res
  vector ages(1,nages);                          			// ages of data
  vector ages4plus(1,nages-1);                   			// non-recruiting ages in the population
  vector years(1,nyrs);                          			// years of data
  int pred_year;                                 			// prediction year
  vector yield_obs(1,nyrs);
 LOCAL_CALCS
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
 END_CALCS
//  vector pred_rec(1,nyrs)
//  vector pred_rec_m3(1,nyrs-3)
  
PARAMETER_SECTION
  init_vector             log_sel_coff(1,nsel_coff,2)	
  init_bounded_matrix     log_sel_dev(1,nsel_changes,1,nsel_coff,-5,5,3) 

  init_number log_mean_F(1)
  init_bounded_dev_vector log_fy_devs(1,nyrs,-10,10,2)  // -.45,5.45original log_fy_coff(1,nyrs,-2.,2.,3)         
  init_number log_mean_recruit(1)			//median recruitment (log scale) 
  init_bounded_dev_vector log_recruit_devs(1,nyrs,-50,50,2)	//deviations in annual recruitment (log scale)
  init_vector log_initpop_devs(1,nages-1,2)	                //deviations in abundance at age in first year (log scale)

  init_bounded_vector log_q1(1,nages_surv1,-15,15,2)		    //survey1-catchability parameter (log scale) 1st period 
  //init_bounded_vector log_q1_2(1,nages_surv1,-15,15,2)		//survey1-catchability parameter (log scale) 2nd period
  init_bounded_vector log_q2(1,nages_surv2,-15,15,2)		    //survey2-catchability parameter (log scale) 1st period 
  //init_bounded_vector log_q2_2(1,nages_surv2,-15,15,2)		//survey2-catchability parameter (log scale) 2nd period
  init_bounded_number log_alfa(0,1,3)				//alfa parameter ssb-rec Ricker (log scale) 
  init_bounded_number log_beta(-5,0,3)				//beta parameter ssb-rec Ricker (log scale) 
  matrix log_sel(1,nyrs,1,nages)
  vector avgsel_fsh(1,nyrs)
  vector log_fy(1,nyrs)
  matrix F(1,nyrs,1,nages)                     // instantaneous fishing mortality
  matrix Z(1,nyrs,1,nages)                     // instantaneous total mortality
  matrix S(1,nyrs,1,nages)                     // survival rate
  matrix N(1,nyrs,1,nages)                     // predicted numbers-at-age
  matrix C(1,nyrs,1,nages)                     // predicted catch-at-age
  matrix U1(1,nyrs_surv1,1,nages_surv1)	       // predicted survey-at-age1
  matrix U2(1,nyrs_surv2,1,nages_surv2)	       // predicted survey-at-age2
	vector neg_log_like(1,10)
  objective_function_value ObjFun
  number recsum
  number initsum
  vector B4(1,nyrs)                            // Vector with B4 biomass for each year
  vector hr(1,nyrs)                            // Vector with Harvest Ratio Yield/B4+ biomass for each year
  number avg_F
  vector yield_pred(1,nyrs)
  sdreport_vector Fbar(1,nyrs)                          // Vector with Fbar(F48) for each year
  sdreport_vector ssb(1,nyrs)
  sdreport_vector ssb_m3(1,nyrs-3)
  sdreport_vector rec(1,nyrs)
  sdreport_vector rec_m3(1,nyrs-3)
  vector pred_rec_m3(1,nyrs-3)
   
PRELIMINARY_CALCS_SECTION
  ages.fill_seqadd(3,1);           // vector of ages
  ages4plus.fill_seqadd(4,1);      // vector of non recruiting ages
  years.fill_seqadd(1961,1);       // fill vector of years with years
  pred_year=years[nyrs]+1;         // year of prediction = last_year +1
  yield_obs=rowsum(elem_prod(obs_catch_at_age,weight_at_age));

INITIALIZATION_SECTION
  log_mean_recruit 10.
  log_mean_F   -1.2
//  log_alfa .4
//  log_beta -4.0
PROCEDURE_SECTION
  get_mortality_and_survivial_rates();
  get_numbers_at_age();
  get_catch_at_age();
  get_survey1_at_age();
  get_survey2_at_age();
  ssb_rec_ricker();
  evaluate_the_objective_function();

FUNCTION get_mortality_and_survivial_rates
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

FUNCTION get_numbers_at_age
  int i, j;

//Loops set initial population values
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

FUNCTION get_catch_at_age
  int i;
  C=elem_prod(elem_div(F,Z),elem_prod(1.-S,N));
  for (i=1;i<=nyrs;i++)
  {                                               // get B4+ Biomass
    yield_pred(i) = C(i)*weight_at_age(i);
    // B4(i)=sum(elem_prod(N(i).sub(2,nages).shift(1),weight_at_age(i).sub(2,nages).shift(1)));
		B4(i) = N(i)(2,nages) * weight_at_age(i)(2,nages);
    hr(i) = yield_obs(i)/B4(i);
  }

FUNCTION get_survey1_at_age
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

FUNCTION get_survey2_at_age
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

FUNCTION ssb_rec_ricker
  ssb=rowsum(elem_prod(elem_prod(N,weight_at_age),mat_at_age))/1000;
  ssb_m3=ssb.sub(1,nyrs-3);
  rec=column(N,1)/1000;
  rec_m3=rec.sub(4,nyrs).shift(1);

  for (int i=1;i<=nyrs-3;i++)
  {
    pred_rec_m3(i)=exp(log_alfa)*ssb_m3(i)*exp(-exp(log_beta)*ssb_m3(i));
  }  
  //-----------------------------------------------------------------------------                                                                  

FUNCTION evaluate_the_objective_function
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

FUNCTION Concentrated_Likelihoods
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

FUNCTION write_mcmc
  // get Fbar
  for(int i=1;i<=nyrs;i++)
    Fbar(i)=mean(F(i)(2,6));

  mcmc_Fbar       << Fbar       << endl;
  mcmc_hr         << hr         << endl;
  mcmc_ssb        << ssb        << endl;
  mcmc_B4         << B4         << endl;
  mcmc_rec        << rec        << endl;
  mcmc_yieldPred  << yield_pred << endl;

REPORT_SECTION
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
FUNCTION write_out
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

TOP_OF_MAIN_SECTION
  gradient_structure::set_NUM_DEPENDENT_VARIABLES(500);
RUNTIME_SECTION
  convergence_criteria 1.e-1,1.e-08
  maximum_function_evaluations 100,100,2500

GLOBALS_SECTION
  #include "admodel.h"
  const double pi = 3.141592654;
  ofstream mcmc_Fbar("mcmc_Fbar.csv");
  ofstream mcmc_hr("mcmc_hr.csv");
  ofstream mcmc_ssb("mcmc_ssb.csv");
  ofstream mcmc_B4("mcmc_B4.csv");
  ofstream mcmc_rec("mcmc_rec.csv");
  ofstream mcmc_yieldPred("mcmc_yieldPred.csv");

FINAL_SECTION
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

