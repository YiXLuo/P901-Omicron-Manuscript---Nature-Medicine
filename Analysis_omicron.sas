*******************************************************
Purpose: 
    Moderna variant study Omicron analysis: 
    describe the study population, baseline characterizes 
    and conducting the analysis for study outcomes.
    Refer to SAP table shells for more details.
*******************************************************;

* Cohort data;
    libname YT "~/omicron3" access=readonly;

* set output passway;
    libname YL "~/Omicron";

*** Call macros;
    %include "~/Omicron/Analysis_omicron_macro.sas";




**********************************
** Characteristics of SARS-CoV-2 test-positive cases and test-negative controls (2-dose analysis), by variant;

    * Create subgroup flag for variant subgroup (matching sets);
    proc sql; create table subgroupflag as
    select a.*, b.variant as variant_subgroup
      from YT.omicron_analysis_revise2 a left join YT.omicron_analysis_revise2(where=(COVID19 ='P')) b 
        on a.pair_id=b.pair_id and a.cohort=b.cohort
     where a.flag_matched=1;
    quit;

option mprint;
    %TNDbaselinetable(subgroupflag,2-dose,Delta,outcome);
    %TNDbaselinetable(subgroupflag,2-dose,Omicron,outcome);

    proc sql; create table table1(drop=rowname temp) as
    select Delta.*, ' ' as separator,
           Omicron.pos_Omicron, Omicron.neg_Omicron, Omicron.p_Omicron, Omicron.ABS_Omicron 
      from Delta left join Omicron on Delta.Rowname=Omicron.Rowname
     order by Delta.order;
    quit;

    data table1_2; set table1;
    pd=.; ad=.; po=.;ao=.;
    pd=p_delta; ad=abs_delta; po=p_omicron; ao=abs_omicron;
    if pd^=. then p_delta=put(pd, 4.2);
    if p_delta in ('<0.0001','0.00') then p_delta='<0.01';
    if po^=. then p_omicron=put(po, 4.2);
    if p_omicron in ('<0.0001','0.00') then p_omicron='<0.01';
    if ad^=. then abs_delta=put(ad, 4.2);
    if ao^=. then abs_omicron=put(ao, 4.2);
    p_abs_delta=tranwrd(catx(' / ', p_delta, abs_delta),'N/A / N/A','N/A');
    p_abs_omicron=tranwrd(catx(' / ', p_omicron, abs_omicron),'N/A / N/A','N/A');
    run;

    proc print data=table1_2; run;





**********************************
** Characteristics of SARS-CoV-2 test-positive cases and test-negative controls (3-dose analysis), by variant;

    %TNDbaselinetable(subgroupflag,3-dose,Delta,outcome);
    %TNDbaselinetable(subgroupflag,3-dose,Omicron,outcome);

    proc sql; create table table2(drop=rowname temp) as
    select Delta.*, ' ' as separator,
           Omicron.pos_Omicron, Omicron.neg_Omicron, Omicron.p_Omicron, Omicron.ABS_Omicron 
      from Delta left join Omicron on Delta.Rowname=Omicron.Rowname
     order by Delta.order;
    quit;

    * HF wanted to combine p-value and ASD together to save some space;
    data table2_2; set table2;
    pd=.; ad=.; po=.;ao=.;
    pd=p_delta; ad=abs_delta; po=p_omicron; ao=abs_omicron;
    if pd^=. then p_delta=put(pd, 4.2);
    if p_delta in ('<0.0001','0.00') then p_delta='<0.01';
    if po^=. then p_omicron=put(po, 4.2);
    if p_omicron in ('<0.0001','0.00') then p_omicron='<0.01';
    if ad^=. then abs_delta=put(ad, 4.2);
    if ao^=. then abs_omicron=put(ao, 4.2);
    p_abs_delta=tranwrd(catx(' / ', p_delta, abs_delta),'N/A / N/A','N/A');
    p_abs_omicron=tranwrd(catx(' / ', p_omicron, abs_omicron),'N/A / N/A','N/A');
    run;

    proc print data=table2_2; run;



**********************************
** Characteristics of SARS-CoV-2 test-positive cases and test-negative controls (1-dose analysis), by variant;
    %TNDbaselinetable(subgroupflag,1-dose,Delta,outcome);
    %TNDbaselinetable(subgroupflag,1-dose,Omicron,outcome);

    proc sql; create table Var_for_1dose(drop=rowname temp) as
    select Delta.*, ' ' as separator,
           Omicron.pos_Omicron, Omicron.neg_Omicron, Omicron.p_Omicron, Omicron.ABS_Omicron 
      from Delta left join Omicron on Delta.Rowname=Omicron.Rowname
     order by Delta.order;
    quit;

    data Var_for_1dose_2; set Var_for_1dose;
    pd=.; ad=.; po=.;ao=.;
    pd=p_delta; ad=abs_delta; po=p_omicron; ao=abs_omicron;
    if pd^=. then p_delta=put(pd, 4.2);
    if p_delta in ('<0.0001','0.00') then p_delta='<0.01';
    if po^=. then p_omicron=put(po, 4.2);
    if p_omicron in ('<0.0001','0.00') then p_omicron='<0.01';
    if ad^=. then abs_delta=put(ad, 4.2);
    if ao^=. then abs_omicron=put(ao, 4.2);
    p_abs_delta=tranwrd(catx(' / ', p_delta, abs_delta),'N/A / N/A','N/A');
    p_abs_omicron=tranwrd(catx(' / ', p_omicron, abs_omicron),'N/A / N/A','N/A');
    run;





**********************************
** Odds ratio of test-positive SARS-CoV-2 cases among mRNA-1273 vaccinated vs. unvaccinated individuals;

** Regroup variables;
    data Regrouped; set subgroupflag;
    * BMI;
    	if BMI>=30 then obese='Yes';
		if .<BMI<30 then obese='No';
		if BMI=. then obese='Unknown';	 
    run;

* Potential confounders will be determined by absolute standardized difference (ASD) >0.1 and p-value<0.1; 
options mprint;

* Control list;
%let varlist = hist_covid_test preventive_care VA_grp obese charlson_wt frailty Specimen_Type_n IC hist_covid_diag;
%let interaction = ;
* 1-dose;
    %TND(Regrouped, Delta,   1-dose, covid_dx, Delta_dx_1);
    %TND(Regrouped, Omicron, 1-dose, covid_dx, Omicron_dx_1);
* 2-dose;
    %TND(Regrouped, Delta,   2-dose, covid_dx, Delta_dx_2);
    %TND(Regrouped, Omicron, 2-dose, covid_dx, Omicron_dx_2);
* 3-dose;
    %TND(Regrouped, Delta,   3-dose, covid_dx, Delta_dx_3);
    %TND(Regrouped, Omicron, 3-dose, covid_dx, Omicron_dx_3);

** Hospitalization;
* Remove COVID infection sets, the infection is not considered case or control here;
    proc sql; create table hospitalsets as
    select distinct a.* from Regrouped a, Regrouped(where=(covid_hosp)) b
     where a.pair_id=b.pair_id and a.cohort=b.cohort;
    quit;*;


* 1-dose;
    %let varlist = hist_covid_test preventive_care obese charlson_wt hist_covid_diag;
        %TND(hospitalsets, Delta,   1-dose, covid_hosp, Delta_ip_1); 
        %TND(hospitalsets, Omicron, 1-dose, covid_hosp, Omicron_ip_1);
        * No enough case, get counts only;

* 2-dose;
    %let varlist = hist_covid_test preventive_care obese charlson_wt hist_covid_diag;
        %TND(hospitalsets, Delta,   2-dose, covid_hosp, Delta_ip_2);
    %let varlist = hist_covid_test preventive_care charlson_wt IC hist_covid_diag;
        %TND(hospitalsets, Omicron, 2-dose, covid_hosp, Omicron_ip_2);

* 3-dose;
    *Delta;
    %let varlist = hist_covid_test preventive_care obese charlson_wt hist_covid_diag;
        %TND(hospitalsets, Delta,   3-dose, covid_hosp, Delta_ip_3, condition_model=0);

    *Omicron;
    %let varlist = hist_covid_test preventive_care charlson_wt IC ;
        %TND(hospitalsets, Omicron, 3-dose, covid_hosp, Omicron_ip_3, condition_model=0);


    * Put results together;
    data blank; variant=' '; run;
    data table3; 
    length variant $25.;
    set results_Delta_dx_1 results_Omicron_dx_1 results_Delta_dx_2 results_Omicron_dx_2 results_Delta_dx_3 results_Omicron_dx_3 
        blank 
        results_Delta_IP_1 results_Omicron_IP_1 results_Delta_IP_2 results_Omicron_IP_2 results_Delta_IP_3 results_Omicron_IP_3;
    drop estimate lower_ci upper_ci;
    run;

    proc print data=table3; run;





**********************************
*** Odds ratio of test-positive SARS-CoV-2 cases by time since vaccination among mRNA-1273 (2 doses) vaccinated vs. unvaccinated individuals;

    * Check time from vaccination distributions in cases;
    data subgroupflag2; length arm $17.; set subgroupflag;
    if arm='Vaccinated' and 14 <= days_interval <= 90 then arm='Vac1 14-90 days';
    if arm='Vaccinated' and 91 <= days_interval <=180 then arm='Vac2 91-180 days';
    if arm='Vaccinated' and 181<= days_interval <=270 then arm='Vac3 181-270 days';
    if arm='Vaccinated' and 271<= days_interval <=365 then arm='Vac4 271-365 days';
    	if BMI>=30 then obese='Yes';
		if .<BMI<30 then obese='No';
		if BMI=. then obese='Unknown';	 
    run;
    proc freq; where cohort='2-dose'; table arm; run;

* Control list;
%let varlist = hist_covid_test preventive_care VA_grp obese charlson_wt frailty Specimen_Type_n IC hist_covid_diag;
    %TNDsubgroup(2-dose,Delta, 17);
    %TNDsubgroup (2-dose,Omicron,17);
    data blank; variant=' '; run;
    data table4; 
    length variant $25.;
    set results_Delta blank results_Omicron;
    drop estimate lower_CI upper_CI;
    run;

    proc print data=table4; run;

    data YL.table4_figuredata_02102022;
    length variant $25. Time $17.;
    set results_Delta results_Omicron;
    Time_order=substr(arm,4,1);
    Time=substr(arm,6, length(arm)-5);
    if lower_CI='0.0' and upper_CI='100.0' then delete;
    run;

	ods graphics / imagefmt=pdf; 
	ods pdf file="~/Omicron/figure2.pdf" compress=0;
    Options nodate nonumber;
    proc sgplot data=YL.table4_figuredata_02102022;
    series x=time y=estimate / datalabel=estimate group=variant;
    scatter x=time y=estimate / yerrorlower=lower_CI yerrorupper=upper_CI group=variant markerattrs=(symbol=circlefilled);
    label estimate="Vaccine effectiveness (VE), %";
    label time="Time since vaccination";
    legend position=(bottom left inside);
    xaxis type=discrete offsetmin=0.1 offsetmax=0.1;
    yaxis values=(0 to 100 by 20);
	format estimate 4.1; 
    run;
	ods pdf close;



*** Time varying for 3-dose;
    data subgroupflag2; length arm $17.; set subgroupflag;
    if arm='Vaccinated' and 14 <= days_interval <= 60 then arm='Vac1 14-60 days';
    if arm='Vaccinated' and 61 <= days_interval <=365 then arm='Vac2 61-365 days';
    	if BMI>=30 then obese='Yes';
		if .<BMI<30 then obese='No';
		if BMI=. then obese='Unknown';	 
    run;

* Control list;
%let varlist = hist_covid_test preventive_care VA_grp obese charlson_wt frailty Specimen_Type_n IC hist_covid_diag;
    %TNDsubgroup (3-dose,Delta, 17);
    %TNDsubgroup (3-dose,Omicron,17);
    data blank; variant=' '; run;
    data table5; 
    length variant $25.;
    set results_Delta blank results_Omicron;
    drop estimate lower_CI upper_CI;
    run;

    proc print data=table5; run;




*** Time varying for 3-dose, without IC;
    data subgroupflag2; length arm $17.; set subgroupflag;
    where IC=0;
    if arm='Vaccinated' and 14 <= days_interval <= 60 then arm='Vac1 14-60 days';
    if arm='Vaccinated' and 61 <= days_interval <=365 then arm='Vac2 61-365 days';
    	if BMI>=30 then obese='Yes';
		if .<BMI<30 then obese='No';
		if BMI=. then obese='Unknown';	 
    run;

* Control list;
%let varlist = hist_covid_test preventive_care VA_grp obese charlson_wt frailty Specimen_Type_n IC hist_covid_diag;
    %TNDsubgroup (3-dose,Delta, 17);
    %TNDsubgroup (3-dose,Omicron,17);
    data blank; variant=' '; run;
    data table6; 
    length variant $25.;
    set results_Delta blank results_Omicron;
    drop estimate lower_CI upper_CI;
    run;

    proc print data=table6; run;

    data YL.table6_figuredata_02102022;
    length variant $25. Time $17.;
    set results_Delta results_Omicron;
    Time_order=substr(arm,4,1);
    Time=substr(arm,6, length(arm)-5);
    if lower_CI='0.0' and upper_CI='100.0' then delete;
    if time='61-365 days' then time='>60 days';
    run;

	ods graphics / imagefmt=pdf; 
	ods pdf file="~/Omicron/figure3.pdf" compress=0;
    Options nodate nonumber;
    proc sgplot data=YL.table6_figuredata_02102022;
    series x=time y=estimate / datalabel=estimate group=variant;
    scatter x=time y=estimate / yerrorlower=lower_CI yerrorupper=upper_CI group=variant markerattrs=(symbol=circlefilled);
    label estimate="Vaccine effectiveness (VE), %";
    label time="Time since vaccination";
    legend position=(bottom left inside);
    xaxis type=discrete offsetmin=0.2 offsetmax=0.2;
    yaxis values=(0 to 100 by 20);
	format estimate 4.1; 
    run;
	ods pdf close;




*********************;
*** Subgroup analysis;

    * Age;
    %let varlist = hist_covid_test preventive_care VA_grp obese charlson_wt frailty Specimen_Type_n IC hist_covid_diag;
    %let interaction =;
        *<65;
        data subgroupflag2; length arm $17.; set Regrouped;
        where age_grp in ('18-44 yo','45-64 yo');
        run;
        %TND(subgroupflag2, Delta,   3-dose, covid_dx, Delta_age_1);
        %TND(subgroupflag2, Omicron, 3-dose, covid_dx, Omicron_age_1);
        *>=65;
        data subgroupflag2; length arm $17.; set Regrouped;
        where age_grp not in ('18-44 yo','45-64 yo');
        run;
        %TND(subgroupflag2, Delta,   3-dose, covid_dx, Delta_age_2);
        %TND(subgroupflag2, Omicron, 3-dose, covid_dx, Omicron_age_2);
        *interaction;
        data subgroupflag3; length arm $17.; set Regrouped;
        age_strata=(age_grp in ('18-44 yo','45-64 yo'));
        run;
        %let interaction = age_strata*arm age_strata*hist_covid_test age_strata*preventive_care age_strata*VA_grp age_strata*obese age_strata*charlson_wt 
                           age_strata*frailty age_strata*Specimen_Type_n age_strata*IC age_strata*hist_covid_diag;
        %TND(subgroupflag3, Delta,   3-dose, covid_dx, Delta_age_i);
        %TND(subgroupflag3, Omicron, 3-dose, covid_dx, Omicron_age_i);


    * Sex;
    %let interaction =;
        *Female;
        data subgroupflag2; length arm $17.; set Regrouped;
        where sex in ('Female');
        run;
        %TND(subgroupflag2, Delta,   3-dose, covid_dx, Delta_sex_1);
        %TND(subgroupflag2, Omicron, 3-dose, covid_dx, Omicron_sex_1);
        *Male;
        data subgroupflag2; length arm $17.; set Regrouped;
        where sex not in ('Female');
        run;
        %TND(subgroupflag2, Delta,   3-dose, covid_dx, Delta_sex_2);
        %TND(subgroupflag2, Omicron, 3-dose, covid_dx, Omicron_sex_2);
        *interaction;
        %let interaction = sex*arm sex*hist_covid_test sex*preventive_care sex*VA_grp sex*obese sex*charlson_wt 
                           sex*frailty sex*Specimen_Type_n sex*IC sex*hist_covid_diag;
        %TND(Regrouped, Delta,   3-dose, covid_dx, Delta_sex_i);
        %TND(Regrouped, Omicron, 3-dose, covid_dx, Omicron_sex_i);


    * Race/ethnicity;
    %let interaction =;
        * Hisp;
        data subgroupflag2; length arm $17.; set Regrouped;
        where race_eth in ('Hispanic');
        run;
        %TND(subgroupflag2, Delta,   3-dose, covid_dx, Delta_hisp);
        %TND(subgroupflag2, Omicron,   3-dose, covid_dx, Omicron_hisp);
        * Non-hisp;
        data subgroupflag2; length arm $17.; set Regrouped;
        where race_eth not in ('Hispanic');
        run;
        %TND(subgroupflag2, Delta,   3-dose, covid_dx, Delta_nonhisp);
        %TND(subgroupflag2, Omicron,   3-dose, covid_dx, Omicron_nonhisp);
        *interaction;
        data subgroupflag3; length arm $17.; set Regrouped;
        race_strata=(race_eth in ('Hispanic'));
        run;
        %let interaction = race_strata*arm race_strata*hist_covid_test race_strata*preventive_care race_strata*VA_grp race_strata*obese race_strata*charlson_wt 
                           race_strata*frailty race_strata*Specimen_Type_n race_strata*IC race_strata*hist_covid_diag;
        %TND(subgroupflag3, Delta,   3-dose, covid_dx, Delta_race_i);
        %TND(subgroupflag3, Omicron, 3-dose, covid_dx, Omicron_race_i);


    * Immunocompromise status;
        * IC;
        %let interaction =;
        data subgroupflag2; set Regrouped;
        where IC;
        run;
        %let varlist = hist_covid_test preventive_care obese charlson_wt frailty Specimen_Type_n IC hist_covid_diag;
        %TNDsubgroup (3-dose,Delta, 10);
        %let varlist = hist_covid_test preventive_care VA_grp obese charlson_wt frailty Specimen_Type_n IC hist_covid_diag;
        %TNDsubgroup (3-dose,Omicron,10);
        *Non IC;
        %let varlist = hist_covid_test preventive_care VA_grp obese charlson_wt frailty Specimen_Type_n IC hist_covid_diag;
        data subgroupflag2; set Regrouped;
        where ^IC;
        run;
        %TNDsubgroup (3-dose,Delta, 10);
        %TNDsubgroup (3-dose,Omicron,10);
        *interaction;
        data subgroupflag2; set Regrouped;
        run;
        %let varlist = hist_covid_test preventive_care obese charlson_wt frailty Specimen_Type_n IC hist_covid_diag;
        %let interaction = IC*arm IC*hist_covid_test IC*preventive_care IC*obese IC*charlson_wt 
                           IC*frailty IC*Specimen_Type_n IC*IC IC*hist_covid_diag
                           IC*age_grp IC*sex IC*race_eth;
        %TNDsubgroup (3-dose,Delta, 10);
        %let varlist = hist_covid_test preventive_care VA_grp obese charlson_wt frailty Specimen_Type_n IC hist_covid_diag;
        %let interaction = IC*arm IC*hist_covid_test IC*preventive_care IC*VA_grp IC*obese IC*charlson_wt 
                           IC*frailty IC*Specimen_Type_n IC*IC IC*hist_covid_diag
                           IC*age_grp IC*sex IC*race_eth;
        %TNDsubgroup (3-dose,Omicron,10);



**********
* Analysis removing those with history of COVID, 
   evaluate VE against omicron infection, stratified by age (<65 years and >=65 years);

    %let interaction =;

    * No history of COVID and age<65;
    data subgroupflag2; set Regrouped;
    where hist_covid_diag=0 and age_grp in ('18-44 yo','45-64 yo');
    run;
    %TNDsubgroup (1-dose,Omicron,10);
    %TNDsubgroup (2-dose,Omicron,10);
    %TNDsubgroup (3-dose,Omicron,10);
            

    * No history of COVID and age>=65;
    data subgroupflag2; set Regrouped;
    where hist_covid_diag=0 and age_grp in ('65-74 yo', '75 and older');
    run;
    %TNDsubgroup (1-dose,Omicron,10);
    %TNDsubgroup (2-dose,Omicron,10);
    %TNDsubgroup (3-dose,Omicron,10);





