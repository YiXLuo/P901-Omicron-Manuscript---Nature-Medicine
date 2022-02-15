*******************************************************
Name: Analysis_omicron_supporting macro

Purpose: 
    Support Moderna variant study Omicron analysis: 
    Generate baseline table per study table shell
    Conduct test negative analysis and variable selection based on SAP
    Conduct subgroup analysis with matching broken

*******************************************************;




***************************************************;
* Macro for calculating standardized differences   ;
***************************************************;


/******************************************************************************/
/* Program : stddiff.sas
/* Purpose : SAS macro to calculate the Standardized Difference
/* Authors: Dongsheng Yang and Jarrod E. Dalton, department of quantitative health sciences and oucomes research
/*             Cleveland Clinic, Cleveland, OH, USA
/* Reference Paper: A Unified approach to measuring the effect size between two groups using SAS
/* Usage : %stddiff(inds = Studydata, groupvar = dex,
/* numvars = age bmi/r glucose,
/* charvars = female surgtype,
/* stdfmt = 8.5,
/* outds = std_result);
/*******************************************************************************/
/* NOTE: All binary variables must be coded as 0 and 1 in the dataset
/* PARAMETERS:
/* inds:       input dataset
/* groupvar:   a binary variable, must be coded as 0 and 1
/* numvars:    a list of continuous variables.
/*             "/r" denotes to use the rank-based mean and SD to calculate Stddiff
/* charvars:   a list of categorical variables. If a variable is a binary categorical variable,
/*             it must be coded as 0 and 1 since we use the level = 0 as the reference level.
/* wtvar:      a weight variable.
/* stdfmt = 8.5 the format of Standardized Difference
/* outds output result dataset
/*********************************************************************************/

/*options  symbolgen mlogic mprint; */

%macro     stddiff( inds = ,groupvar =,numvars = ,charvars = ,wtvar = ,stdfmt = 8.4,outds = stddiff_result ); 

/* create a table to store stddiff */
proc sql; 
   create table &outds.  
       (VarName char(32), 
           Stddiff char (10)
       ); 
quit; 

/* delete records if the group variable is missing */

data base_data; 
     set &inds.; 
     where &GroupVar. ne .; 
run; 

/* remove leading or tailing blanks */
%let groupvar = %sysfunc(strip(&GroupVar.)); 

                    /****************************************/
                    /* part 1: compare continuous variables */

%if %length(&numvars.) > 0 %then %do; 

/* remove multiple blanks and get the total number of continuous variables */
    %let numvar = %sysfunc(compbl(&numvars.)); 
    %let numvar = %sysfunc(strip(&numvar.)); 
    %let n_convar = %sysfunc(countc(&numvar.,' ')); 
    %let n_convar = %eval(&n_convar. + 1); 

/* summarize variables one-by-one */
    %do ii = 1 %to &n_convar.; 
        %let convar = %sysfunc(scan(&numvar.,&ii.,' ')); 

    /* if requires rank-based mean and std for skewed variables */
        %if %index(&convar., /r) > 0 %then %do; 
            %let convar = %sysfunc(scan(&convar.,1,'/')); 
            %let convar = %sysfunc(strip(&convar.)); 

            data temp_1; 
                 set base_data (keep = &groupvar. &convar. &wtvar.); 
            run; 

    /* rank a variable */
            proc rank data=temp_1 out=temp_2; 
                  var &convar.; 
                ranks rank_&convar.; 
            run; 

    /* get ranked-mean and sd */

            proc means data = temp_2;
                class &groupvar.;
                var rank_&convar.;
                weight &wtvar.;
                output out = temp_3 mean = _mean_  std = _std_;
            run;

            data  temp_3;
                set temp_3;
                where _type_ = 1;
            run;

            proc sort data = temp_3;
                by &groupvar.;
            run;
             %end; 
   
    /* for normal-distributed variable */

        %else %do; 
            %let convar = %sysfunc(strip(&convar.)); 
            data temp_1; 
                 set base_data (keep = &groupvar. &convar. &wtvar.); 
            run; 
            data temp_2; 
                 set temp_1; 
            run; 

        /* get mean and sd */

            proc means data = temp_2;
                class &groupvar.;
                var &convar.;
                weight &wtvar.;
                output out = temp_3 mean = _mean_  std = _std_;
            run;

            data  temp_3;
                set temp_3;
                where _type_ = 1;
            run;

            proc sort data = temp_3;
                by &groupvar.;
            run;

            %end; 

/* calculate stddiff */   
       proc sql; 
            create table temp_4 as  
                select (a._mean_ - b._mean_)/ 
                sqrt((a._std_**2 + b._std_**2)/2) as d 
                from temp_3(where = (&groupvar = 1)) as a, 
                       temp_3(where = (&groupvar = 0)) as b; 
       quit; 
         
       data temp_5; 
               set temp_4; 
            stddiff = compress(put(d,&stdfmt.)); 
            keep stddiff; 
        run; 

    /* insert into std table */
       proc sql noprint; 
            select stddiff into: std_value from temp_5; 
            insert into &outds.  values("&convar.", "&std_value."); 
       quit; 

    /* delete temporary data sets */

       proc datasets lib = work nodetails nolist; 
        delete  temp_1 - temp_5; 
       quit; 
          %end;  
%end; 

        /**********************************************/
        /* part 2: compare categorical variables      */

%if %length(&charvars.) > 0 %then %do; 
    %let n_charvar = %sysfunc(countw(&charvars.)); 

/* get column percents for each levels of the variable by the group */
    %do jj = 1 %to &n_charvar.; 
           %let char_var = %scan(&charvars., &jj.); 
           %let char_var = %sysfunc(strip(&char_var.)); 
          data temp_1; 
               set base_data (keep = &groupvar. &char_var. &wtvar.); 
          run; 
    
          proc sql; 
               create table temp_2 as 
               select distinct &char_var. as &char_var.
               from temp_1
            where &char_var. is not missing; 
          quit; 

          proc sql noprint; 
               select count(*) into :_mylevel_ from temp_2; 
          quit; 

           %let _mylevel_ = %sysfunc(strip(&_mylevel_.)); 

          data temp_3; 
               set temp_2; 
                   do &groupvar. = 0,1 ; 
                output; 
                   end; 
          run;

        ods output CrossTabFreqs = temp_4; 
          proc freq data = temp_1; 
               table &char_var. * &groupvar.; 
            %if %length(&wtvar.) > 0 %then %do;
                weight &wtvar.;
                %end;
          run; 

          proc sql; 
               create table  temp_5 as 
               select a.*, b.ColPercent 
               from temp_3 as a 
               left join temp_4 as b 
               on     a.&groupvar. = b.&groupvar. and  
                 a.&char_var. = b.&char_var.; 
          quit; 

          data temp_6; 
               set temp_5; 
               if ColPercent = . then ColPercent = 0; 
          run; 

          proc sort data = temp_6 out = catfreq; 
               by &groupvar. &char_var.; 
          run; 
  
          proc datasets lib = work nodetails nolist; 
               delete  temp_1 - temp_6; 
          quit; 

/* if a categorical variable only has one level: 0 or 1 */
/* stddiff = 0 */
        %if &_mylevel_. = 1 %then %do; 
              proc sql noprint; 
                   insert into &outds.  values("&char_var.", "0"); 
              quit; 
           %end; 

/* if a categorical variable  has two level: 0 and 1 */
/* it is a binary variable, using two sample proportation formula */
        %else %if &_mylevel_. = 2 %then %do; 

              data temp_7; 
                   set catfreq; 
                  where &char_var. = 1; 
                   ColPercent = ColPercent/100; 
              run; 

              proc sql; 
                   create table temp_8 as  
                   select (a.ColPercent - b.ColPercent)/(sqrt((a.ColPercent*(1- 
                        a.ColPercent) +  
                         b.ColPercent*(1-b.ColPercent))/2)) as d 
                   from temp_7(where = (&groupvar = 1)) as a, 
                      temp_7(where = (&groupvar = 0)) as b; 
              quit; 
     
              data temp_9; 
                  set temp_8; 
                  stddiff = compress(put(d,&stdfmt.)); 
                keep stddiff; 
            run; 

              proc sql noprint; 
                   select  stddiff into: std_value from temp_9; 
                       insert into &outds.  values("&char_var.", "&std_value."); 
              quit; 

              proc datasets lib = work nodetails nolist; 
                   delete  temp_7 temp_8 temp_9; 
              quit; 
        %end; 
/* if a categorical variable  has more than two level such as a, b and c */
        %else %if &_mylevel_. > 2 %then %do; 
               %let _k_ = %eval(&_mylevel_. - 1); 
               %let _k_ = %sysfunc(strip(&_k_.)); 
              data temp_7; 
                   set catfreq; 
                  by &groupvar.; 
                   if last.&groupvar. then delete; 
                   ColPercent = ColPercent/100; 
              run; 

              proc sql noprint; 
                   select ColPercent into :tlist separated by ' '  
                from temp_7 where &groupvar. = 1; 

                   select ColPercent into :clist separated by ' '  
                from temp_7 where &groupvar. = 0; 
              quit; 

/* vector T, C and T-C */
              data t_1; 
                   array t{*}  t1- t&_k_.   (&tlist.); 
                   array c{*}  c1- c&_k_.   (&clist.); 
                   array tc{*} tc1 - tc&_k_. ; 
                   do i = 1 to dim(t); 
                    tc{i} = t{i} - c{i}; 
                   end; 
               drop i; 
              run; 

/* each column has one element of a S covariance matrix (k x k) */

            %let _dm = ; 
            %let _dm = %eval(&_k_.*&_k_.); 
              data covdata; 
                   array t{*}  t1- t&_k_.  (&tlist.); 
                   array c{*}  c1- c&_k_.   (&clist.); 
                   array cv{&_k_.,&_k_.} x1 -x&_dm.; 
                   do i = 1 to &_k_.; 
                    do j = 1 to &_k_.; 
                         if i = j then do; 
                              cv{i,j} = 0.5*(t{i}*(1-t{i}) + c{i}*(1-c{i})); 
                              end; 
                         else do; 
                              cv{i,j} = -0.5 * (t[i] * t[j] + c[i] * c[j]); 
                              end; 
                        if cv{&_k_.,&_k_.] ne . then output; 
                    end; 
                  end; 
              run; 

              proc transpose data = covdata(keep = x1 -x&_dm.) out = covdata_1; 
              run; 

              data covdata_2; 
                   set covdata_1; 
                   retain id gp 1; 
                   if mod(_n_ - 1,&_k_.) = 0 then gp = gp + 1; 
              run; 

              proc sort data = covdata_2 ; 
                   by gp id; 
              run;   

            data covdata_3; 
                   set covdata_2; 
                   by gp id; 
                   retain lp; 
                   if first.gp then lp = 0; 
                   lp = lp+1; 
              run; 

/* transpose to a S variance-covariance matrix format */
           
              data covdata_4; 
                   set covdata_3; 
                   retain y1-y&_k_.; 
                   array cy{1:&_k_.} y1-y&_k_.; 
                   by gp id; 
                   if first.gp then do; 
                    do k = 1 to &_k_.; 
                         cy{k} = .; 
                    end; 
                   end; 
                   cy{lp} = col1; 
                   if last.gp then output; 
                   keep y:; 
              run; 

/* get inverse of S matrix */
          data A_1; 
           set covdata_4; 
           array _I{*} I1-I&_k_.; 
           do j=1 to &_k_.; 
            if j=_n_ then _I[j]=1;  
            else _I[j]=0; 
           end; 
           drop j; 
          run; 

/* solve the inverse of the matrix */

  %macro inv; 
        %do j=1 %to &_k_.; 
            proc orthoreg data=A_1 outest=A_inv_&j.(keep=y1-y&_k_.) 
                 noprint singular=1E-16; 
                 model I&j=y1-y&_k_. /noint; 
            run; 
            quit; 
        %end; 

           data A_inverse; 
            set %do j=1 %to &_k_.; 
             A_inv_&j 
         %end;; 
           run; 
  %mend; 
  %inv; 

          proc transpose data=A_inverse out=A_inverse_t; 
          run; 

   /* calculate the mahalanobis distance */
          data t_2; 
               set A_inverse_t; 
               array t{*}  t1- t&_k_.  (&tlist.); 
               array c{*}  c1- c&_k_.  (&clist.); 
               i = _n_; 
               trt = t{i}; 
               ctl = c{i}; 
               tc = t{i} - c{i}; 
          run; 
 
        data t_3; 
               set t_2; 
               array aa{&_k_.} col1 - col&_k_.; 
               array bb{&_k_.} bb1- bb&_k_.; 
               do i = 1 to &_k_.; 
                bb{i} = aa{i}*tc; 
               end; 
          run; 

          proc summary data = t_3 ; 
               var bb1-bb&_k_.; 
               output out = t_4 sum =; 
          run; 

          data t_5; 
               merge t_1 t_4; 
               array d1{*} tc1- tc&_k_. ; 
               array d2{*} bb1-bb&_k_.; 
               array d3{*} y1-y&_k_.; 
               do i = 1 to &_k_.; 
                   d3{i} = d1{i}*d2{i}; 
               end; 
               d = sqrt(sum(of y1-y&_k_.)); 
               stddiff = compress(put(d,&stdfmt.));      
               keep stddiff; 
          run; 

          proc sql noprint; 
               select  stddiff into: std_value from t_5; 
               insert into &outds.  values("&char_var.", "&std_value."); 
          quit; 
   
          proc datasets lib = work nodetails nolist; 
               delete  covdata covdata_1 covdata_2 covdata_3 covdata_4 
              A_1 A_inverse A_inverse_t t_1 t_2 t_3 t_4 t_5
             A_inv_:; 
          quit; 
       %end; 
    %end; 
%end; 

proc datasets lib = work nodetails nolist; 
  delete Catfreq  Base_data temp_7; 
quit; 

proc print data = &outds.; 
    title 'Calculated Standardized Difference';
run; 

title;

%mend;





*******************************************;
* Macro for descriptive table and p-values ;
*******************************************;

%MACRO stat_TABLE(DSN=_last_, ID=, BY=, VAR=, TYPE=, OUTDOC=, OUTDAT=, LIST=N, PRINT=Y, NUMBER=N, 
             CSTATS=n mean sd median quartiles range, DSTATS=n percent, pfoot=n,
             SURV=, SCENSOR=0, TOTAL=Y, PVALUES=Y, PTYPE=, CI=, POP=, INCMISS=N, INCMISS1=, 
             COMMENTS=, DVAR=, DLINE=, TTITLE1=, TTITLE2=, TTITLE3=, TTITLE4=, 
             DECPCT=1, DECSTD=2, DECMEAN=1, DECMEDIAN=1, DECMIN=1, DECMAX=1, DECQ1=1, DECQ3=1, 
             FOOT=, FOOT2=, FOOT3=, FOOT4=, FOOT5=, DATE=Y, PAGE=portrait, RULES=groups, FRAME=hsides, 
             TITLESZ=10, BODYSZ=10, FOOTSZ=10, TITLEBLD=Y, TITLEFNT=Times New Roman, HEADFNT=Times New Roman,
             BODYFNT=Times New Roman, CIWD=, PVALWD=100, LEVELWD=, DATAWD=, SPACE=, DEBUG=N, PDEC=4, pcttype=col, PCTPRSNT=1, MEANPRSNT=1
             );

%let validvn=%sysfunc(getoption(validvarname));
%let sdate=%sysfunc(getoption(nodate));
%let snotes=%sysfunc(getoption(nonotes));
%let snumb=%sysfunc(getoption(nonumber));
%let number=%upcase(&number);

%let debug=%upcase(&debug);
%if (&debug=Y) %then %do;  options mprint mtrace macrogen notes linesize=132 ps=58; %end;
               %else %do;  options nonotes nomprint nomacrogen nomtrace nosymbolgen nomlogic linesize=132 ps=58; %end;

options nodate %if &number=Y %then %do; number %end; %else %do; nonumber %end; validvarname=v7 ;


/* creates defaults and conversions */
%if (&space=) %then %let space=2;
%else %if (&space=1)  %then %let space=0;
%else %if (&space=2)  %then %let space=4;

%let dsn=%upcase(&dsn);
%let by=%upcase(&by);
%let pvalues=%upcase(&pvalues);
%let total=%upcase(&total);
%let incmiss=%upcase(&incmiss);
%let list=%upcase(&list);
%let print=%upcase(&print);
%let titlebld=%upcase(&titlebld);
%let date=%upcase(&date);
%let pfoot=%upcase(&pfoot);
%let pcttype=%upcase(&pcttype);
  %if "&pcttype"="PCT_ROW" | "&pcttype"="ROW" | "&pcttype"="ROW_PCT" %then %let pct_type=PCT_COL;
  %if "&pcttype"="PCT_COL" | "&pcttype"="COL" | "&pcttype"="COLUMN" | "&pcttype"="COL_PCT" %then %let pct_type=PCT_ROW;
%if &pvalues=N %then %let pfoot=N;
%if &pfoot=Y %then %do;
  %let pfnum=0; %let pt1=0; %let pt2=0; %let pt3=0; %let pt4=0; %let pt5=0; %let pt6=0; %let pt7=0; %let pt8=0;
  %let pnm1=Chi-Square; %let pnm2=Fisher Exact; %let pnm3=Kruskal Wallis; %let pnm4=Exact Kruskal Wallis;
  %let pnm5=Wilcoxon; %let pnm6=Exact Wilcoxon; %let pnm7=ANOVA F-Test; %let pnm8=Log-Rank; %let pft= ;
%end;

%if (&by=) %then %do; 
  %let pvalues=N;
  %let total=N;
  %let by=_MYBY_;
  %let noby=1;
  %let ci=;
%end;
%else %do; %let noby=0; %end;

      %if (&decpct=4) %then %do; %let decpct=0.0001; %let decpctf=8.4; %end;
%else %if (&decpct=3) %then %do; %let decpct=0.001;  %let decpctf=7.3; %end;
%else %if (&decpct=2) %then %do; %let decpct=0.01;   %let decpctf=6.2; %end;
%else %if (&decpct=1) %then %do; %let decpct=0.1;    %let decpctf=5.1; %end;
%else %if (&decpct=0) %then %do; %let decpct=1;      %let decpctf=3.0; %end;
%else %do; %let decpct=0.1;    %let decpctf=5.1; %end;

      %if (&decstd=4)  %then %do; %let decstd=0.0001; %let decstdf=13.4; %end; 
%else %if (&decstd=3)  %then %do; %let decstd=0.001;  %let decstdf=12.3; %end; 
%else %if (&decstd=2)  %then %do; %let decstd=0.01;   %let decstdf=11.2; %end; 
%else %if (&decstd=1)  %then %do; %let decstd=0.1;    %let decstdf=10.1; %end; 
%else %if (&decstd=0)  %then %do; %let decstd=1;      %let decstdf=8.0 ; %end;
%else %do; %let decstd=0.01;   %let decstdf=11.2; %end; 

      %if (&decmean=4)  %then %do; %let decmean=0.0001; %let decmeanf=13.4; %end; 
%else %if (&decmean=3)  %then %do; %let decmean=0.001;  %let decmeanf=12.3; %end; 
%else %if (&decmean=2)  %then %do; %let decmean=0.01;   %let decmeanf=11.2; %end; 
%else %if (&decmean=1)  %then %do; %let decmean=0.1;    %let decmeanf=10.1; %end; 
%else %if (&decmean=0)  %then %do; %let decmean=1;      %let decmeanf=8.0 ; %end;
%else %do; %let decmean=0.1;    %let decmeanf=10.1; %end; 

      %if (&decmedian=4) %then %do; %let decmedian=0.0001;  %let decmedianf=13.4; %end;
%else %if (&decmedian=3) %then %do; %let decmedian=0.001;   %let decmedianf=12.3; %end;
%else %if (&decmedian=2) %then %do; %let decmedian=0.01;    %let decmedianf=11.2; %end;
%else %if (&decmedian=1) %then %do; %let decmedian=0.1;     %let decmedianf=10.1; %end;
%else %if (&decmedian=0) %then %do; %let decmedian=1;       %let decmedianf=8.0 ; %end;
%else %do; %let decmedian=0.1;     %let decmedianf=10.1; %end;

      %if (&decmin=4) %then %let decminf=13.4; 
%else %if (&decmin=3) %then %let decminf=12.3; 
%else %if (&decmin=2) %then %let decminf=11.2; 
%else %if (&decmin=1) %then %let decminf=10.1; 
%else %if (&decmin=0) %then %let decminf=8.0 ; 
%else %let decminf=10.1; 

      %if (&decmax=4) %then %let decmaxf=13.4; 
%else %if (&decmax=3) %then %let decmaxf=12.3; 
%else %if (&decmax=2) %then %let decmaxf=11.2; 
%else %if (&decmax=1) %then %let decmaxf=10.1; 
%else %if (&decmax=0) %then %let decmaxf=8.0 ; 
%else %let decmaxf=10.1; 

      %if (&decq1=4) %then %let decq1f=13.4; 
%else %if (&decq1=3) %then %let decq1f=12.3; 
%else %if (&decq1=2) %then %let decq1f=11.2; 
%else %if (&decq1=1) %then %let decq1f=10.1; 
%else %if (&decq1=0) %then %let decq1f=8.0 ; 
%else %let decq1f=10.1; 

      %if (&decq3=4) %then %let decq3f=13.4; 
%else %if (&decq3=3) %then %let decq3f=12.3; 
%else %if (&decq3=2) %then %let decq3f=11.2; 
%else %if (&decq3=1) %then %let decq3f=10.1; 
%else %if (&decq3=0) %then %let decq3f=8.0 ; 
%else %let decq3f=10.1; 



      %if (&pdec=4) %then %let ppct=0.0001;
%else %if (&pdec=3) %then %let ppct=0.001;
%else %if (&pdec=2) %then %let ppct=0.01;
%else %if (&pdec=1) %then %let ppct=0.1;
%else %let ppct=0.0001;

%let errors=0;
%let bytype=char;
%let byf=$char5.;
%if %scan(&surv,2)^= %then %let snum=1;
%else %let snum=0;
%if %scan(&scensor,2)^= %then %let scen=1;
%else %let scen=0;

%let cn=0; %let cmn=0; %let csd=0; %let cmd=0; %let cq=0; %let cr=0; %let cnmiss=0;
%do i=1 %to 7;
  %let myscan=%upcase(%scan(&cstats,&i));
  %if (&myscan^=) %then %do;
    %if &myscan=N %then %do; %let cn=1; %end; 
    %if &myscan=MEAN %then %let cmn=1; 
    %if &myscan=SD %then %let csd=1; 
    %if &myscan=MEDIAN %then %let cmd=1;
    %if &myscan=QUARTILES %then %let cq=1;
    %if &myscan=RANGE %then %let cr=1;
    %if &myscan=NMISS %then %let cnmiss=1;
  %end;
  %else %let i=7;
%end;

%if &cmn=1 & &csd=1 %then %let cms=1; 
%else %if &cmn=1 %then %let cms=2; 
%else %if &csd=1 %then %let cms=3; 
%else %let cms=0;

%let dn=0; %let dp=0;
%do i=1 %to 2;
  %if (%scan(&dstats,&i)^=) %then %do;
  %if %upcase(%scan(&dstats,&i))=N %then %let dn=1;
  %if %upcase(%scan(&dstats,&i))=PERCENT %then %let dp=1;
  %end;
%end;

%if &dn=1 & &dp=1 %then %let dnp=1; 
%else %if &dn=1 %then %let dnp=2; 
%else %if &dp=1 %then %let dnp=3; 
%else %do; %put ERROR: incorrect statistics chosen in DSTATS, defaults will be used; %let dnp=1; %end;


proc format;
  picture fpvalue low-high='9.9999';

/* creates an analysis master file */
%if ^(%sysfunc(exist(&dsn))) %then %do; %let errors=1;  %let errorwhy=Dataset &dsn does not exist; %end;
%else %do;

data _master (keep=&id &var &surv &by _by_ );
  set &dsn;

  /* creates a character variable to replace the by variable */
  length _by_ $ 40;
  %if &noby=1 %then %do; &by='Total'; %end; 
  _by_ = trim(&by); 

  /* selects the correct analysis population */
  %if ("&pop"^="")   %then %do; if &pop;  %end;
run;

%if (&type=) %then %do; %let errors=1;  %let errorwhy=Variable TYPE was not defined; %end; %let opn=%sysfunc(open(_master));  
%if &opn %then %do;
  %if (%sysfunc(attrn(&opn,NOBS))=0) %then %do; %let errors=1;  %let errorwhy=Error creating master dataset from dataset &dsn; %end;
  %let rc=%sysfunc(close(&opn));
%end;
    
%end; 
  
%if &errors^=1 %then %do;
  
/* creates the template for the listings */
proc template;
  define style lsttable;
  style table /
     frame=&frame
     cellpadding=4
     cellspacing=2
     rules=&rules
     asis=on
     borderwidth=2;
  style data /
     font_face="&bodyfnt"
     asis=on
     font_size=&bodysz.pt;
  style Body /
     font_face="&bodyfnt"
     asis=on
     font_size=&bodysz.pt;
  style SystemTitle /
     font_face="&bodyfnt"
     font_weight=bold
     asis=on
     font_size=&bodysz.pt;
  style Header /
     font_face="&bodyfnt"
     asis=on
     font_size=&bodysz.pt;
  style BodyDate /
     font_face="&bodyfnt"
     asis=on
     font_size=&bodysz.pt;
  style Byline /
     font_face="&bodyfnt"
     asis=on
     font_size=&bodysz.pt;
  style SystemFooter /
     font_face="&bodyfnt"
     asis=on
     font_size=&bodysz.pt;
  style SysTitleAndFooterContainer /
     font_face="&bodyfnt"
     asis=on
     font_size=&bodysz.pt;
  style Obs /
     font_face="&bodyfnt"
     asis=on
     font_size=&bodysz.pt;
  style IndexItem /
     font_face="&bodyfnt"
     asis=on
     font_size=&bodysz.pt;
  style Rowheader /
     font_face="&bodyfnt"
     asis=on
     font_size=&bodysz.pt;
  end;
  run;

options linesize=120 pagesize=54;

%if %length(%scan(&outdoc,-1,%str(.)))=3 %then %let doctype=%scan(&outdoc,-1,%str(.)); %else %let doctype=0;

%if &list=Y %then %do;
  %if "&doctype"="0" %then %do; %let docnm=&outdoc._lst.doc; %end;
  %else %do; %let docnm=%scan(&outdoc,1,.)_lst.&doctype; %end;

  ods listing close;
  ods rtf file="&docnm"  style=lsttable;

  data _tmp_;
    set _master;
  
  proc sort data=_tmp_;
    by &by &id;

  %let docname=%scan(&outdoc,-1,/);

  proc print data=_tmp_ label uniform split='#';
    %if &noby=0 %then %do; by &by; %end;
    var &id &var;
    title "Listing of Data for &docname created on &sysdate9";
  run;
  
  ods trace off;

  ods rtf close;
  ods listing;
  %put Created File: &docnm;
%end;

options linesize=72 pagesize=60;

/* creates tables of variable characteristics */
/* as macro variables                         */
proc sql;
  create table _c_ as
    select *
     from sashelp.vcolumn
     where libname="WORK"
      and memname="_MASTER";

  create table _c2_ as
    select *
     from sashelp.vtable
     where libname="WORK"
      and memname="_MASTER";
  quit;

/* determines the number of observations in the analysis file */
     /* nobs = number of observations in the analysis file */
data _null_;
  set _c2_;
  call symput('nobs',trim(left(put(nobs,5.))));

/* gets variable type, format, and label for the by variable */
     /* bytype = char or num = by variable type */
     /* byf = format for the by variable */
     /* bylbl = label for the by variable */
     
 data _null_;
  set _c_;
  if (label=' ') then label=name;
  if (upcase(name)="&by") then do;
        call symput('bytype',trim(left(type)));
        call symput('byf',trim(left(format)));
        call symput('bylbl',trim(left(label)));
        call symput('bylength',trim(left(put(length,3.))));
  end;
  run;
   
/* distribution of by variable (including missing)     */
/* and creates macro variables of the different levels */
  proc freq data=_master noprint;
    table &by / out=_d01 missing ;

  proc sort data=_d01; by &by; run;

      /* by1, by2, etc are identifiers of the by levels */
      /* byn1, ... are the totals in each by group */
      /* nby is the number of levels of the by variable */

  data _null_;
    set _d01;
    %if (&byf^=) %then %do;
      call symput("fby" || trim(left(_n_)),trim(put(&by,&byf.)));
    %end;
    %else %do; call symput("fby" || trim(left(_n_)),trim(&by)); %end;
    call symput("by" || trim(left(_n_)),trim(&by));
    call symput("byn" || trim(left(_n_)),trim(left(count)));
    call symput("nby",trim(left(_n_)));
  run;
  data _null_;
    set _d01;
    %if &bytype=num %then %do; where &by^=.; %end;
    %if &bytype=char %then %do; where &by^=""; %end;
    call symput("nby_miss",trim(left(_n_)));
  run;

/** does the analyses for each variable **/
%let num=1;                                    /* num = the number of the variable being processed */
%let v1=%upcase(%scan(&var,&num));             /* v1 = the variable name for the current discrete variable */
%let ind=%scan(&type,&num);                    /* ind=1 for cont, 2 for discrete, 3 for both */

%if ("&ptype"^="") %then %let pval=%scan(&ptype,&num);                  /* ptype=type of pvalue desired */
%if ("&ptype"="") %then %if &ind=1 | &ind=7 | &ind=8 %then %let pval=3;
  %else %if &ind=3 %then %let pval=5;
  %else %if &ind=4 %then %let pval=8; 
  %else %let pval=1;
  
%if &pfoot=Y %then %do;
%if &pval>0 %then %do;
%if &&pt&pval=0 %then %do;
  %let pt&pval=1;
  %let pfnum=%eval(&pfnum+1);
  %let pnum&pval=&pfnum;
  %let pft=~{super &pfnum} &&pnm&pval;
%end; %end;
%end; 

data _mst; run;
   
%do %while (&v1^=);

  /* gets the variable type, format, and label */
  %let v1type=;
  %let v1f=;
  %let v1lbl=;
  %let hasit=0;
  
  data _null_;
    set _c_;
    if (label=' ') then label=name;
    if (upcase(name)="&v1") then do;
          call symput('v1type',trim(left(type)));        /* v1type = current variable type   */
          call symput('v1f',trim(left(format)));         /* v1f    = current variable format */
          call symput('v1lbl',trim(left(%str(label))));  /* v1lbl  = current variable label  */
          call symput('hasit',1); 
    end;
    run;
  
   %if &hasit=0 %then %do; %let errorvar=1; %put WARNING: Variable &v1 not found in dataset &dsn; %end;
   %else %let errorvar=0; 

   %do i=1 %to 100;
     %if (%scan(&incmiss1,&i)^=) %then %do;
       %if %scan(&incmiss1,&i)=&num
       %then %do; %let mis=1; %let i=100; %end;
       %else %do; %let mis=0; %end;
     %end;
     %else %do; %let mis=0; %let i=100; %end;
   %end;
   
   %if &nby>=2 %then %do;
   %do i=1 %to 100;
     %if (%scan(&ci,&i)^=) %then %do;
       %if %scan(&ci,&i)=&num
       %then %do;
         %if &type=1 %then %do; %let cip=1; %end;
         %else %do; %let cip=2; %end;
         %let i=100;
       %end;
       %else %do; %let cip=0; %end;
     %end;
     %else %do; %let cip=0; %let i=100; %end;
   %end;
   %end; %else %do; %let ci=; %let cip=0; %end;

  %if &errorvar=0 %then %do;

  /** analyses for discrete variables */
  %if &ind=2 or &ind=3 or &ind=5 or &ind=6 or "&ind"="p1" or "&ind"="p2" or "&ind"="p3"  %then %do;

    /* overall distribution of analysis variable */
    data _master;
      set _master;
      format &by;

    proc freq data=_master noprint;
      table &v1 / out=_d02
                  %if (&incmiss=Y) | (&mis=1)
                  %then %do;  missing  %end; ;
    run;

    
/* creates a matrix of all possible combinations of the analysis variable and the by variable */
/* including missing values                                                                   */
           /* l1-l8 are the levels of the by variable */
           /* n1-n8 are the total counts at each by level */
           /* nall is the total number of observations */

    data _level (keep=l1-l&nby. n1-n&nby. nall _merge);
      set _d01 end=eof;
      %if (&bytype=char) %then %do; length l1-l&nby. $ &bylength.; %end;
      %if (&byf^=)       %then %do; format &by l1-l&nby. &byf; %end;
     
        retain l1-l&nby. n1-n&nby. ;
        _merge=1;
        %do i=1 %to &nby;
          if _n_=&i then do;
            l&i. = &by;
            n&i. = count;
          end;
      %end;

      if eof then do;
         nall=sum(of n1-n&nby.);
         output;
      end;

    /* puts the overall by level values on each record */
    data _dall;
      set _d02;
      if (_n_=1) then set _level;

    data _dall (keep=&by &v1);
      set _dall end=eof;
      %if (&bytype=char) %then %do; length &by. $ &bylength.; %end;
      %if (&byf^=) %then %do; format &by &byf; %end;
      %if (&bytype=num) %then %do;
         %do i = 1 %to &nby.;
           if l&i.^=. then do; &by=l&i.; output; end;
         %end;
      %end;
      %if (&bytype=char) %then %do; %do i=1 %to &nby.;
        if l&i.^=' ' then do; &by=l&i.; output; end;
      %end; %end;
      if eof then call symput ('ny',trim(left(_n_)));
    run;

    %let nx=&nby;

    /* sets p-values initially to missing */
    data _d2;
      p_pchi=.; p_exact2=.; output;

    /* subgroups distributions excluding missing */
    proc freq data=_master noprint;
      table &by * &v1 / out=_d1 outpct
                        %if (&incmiss=Y) | (&mis=1)
                        %then %do;  missing %end;

      /* the following criteria are for whether or not
         to do an exact test */
                      %if &nby_miss>1 & "&pvalues"="Y" %then %do;  
                        %if (%eval(&ny<=8) & %eval(&nx<=5) & %eval(&nobs<100)) |
                            (%eval(&ny<=5) & %eval(&nx<=3)) & %eval(&nobs<=50)
                        %then %do; exact %end;
                        %else %do; chisq %end; 
                      %end;  
                        ;
       %if &nby_miss>1 & "&pvalues"="Y"  %then %do;                   
          output out=_d2 n nmiss pchi 
              %if (%eval(&ny<=8) & %eval(&nx<=5) & %eval(&nobs<100)) |
                  (%eval(&ny<=5) & %eval(&nx<=3)) & %eval(&nobs<=50)
              %then %do; exact %end;
              %else %do; chisq %end; ;
       %end;
    run;

    data _d1;
      set _d1;
      %if (&byf^=) %then %do; format &by &byf; %end;
    run;      
      
    data _d2;  *create a by variable to merge on later on;
      set _d2; _merge=1; run;

    /* determines if there are enough levels to run the Kruskal-Wallis tests */  
   data _check;
      set _d1;
      if (count=0) then delete;
      %if (&bytype=num) %then %do;
         if (&by=.) then delete;
      %end;
      %else %do;
         if (&by=' ') then delete;
      %end;

    proc sort data=_check;
      by &by;

    data _check;
      set _check;
      by &by;
      if (first.&by);

    %let varby=1;

    data _check;
      set _check;
      if (_n_>1) then call symput('varby','2');
      run;

    /* Kruskal-Wallis p-value excluding missing */
    %if (&v1type=num) & (&varby^=1)  & "&pvalues"="Y" %then %do;
        proc npar1way data=_master wilcoxon noprint;
          var &v1;
          class &by;

          * the following criteria are for whether;
          * or not to do an exact test;
          %if (%eval(&ny<=8) & %eval(&nx<=5) & %eval(&nobs<50)) |
              (%eval(&ny<=5) & %eval(&nx<=3)) | %eval(&nobs<=20)
          %then %do; exact wilcoxon; %end;
          output out=_f1 wilcoxon anova;
    %end;
    %else %do;
         data _f1;
           p_kw=.;
           xp_kw=.;
           p_f=.;
           p2_wil=.;
           xp2_wil=.;
           output;
    %end;
  
    data _f1; set _f1; _merge=1; run;
   
    proc sort data=_dall;
      by &v1 &by;

    proc sort data=_d1;
      by &v1 &by;

    data _final;
      merge _dall (in=a) _d1 (in=b);
      by &v1 &by;

    data _tmp;
      set _master;
      if (_by_^=' ');

    proc freq data=_tmp noprint;
      table &v1 / out=_d02  %if (&incmiss=Y) |
                                (&mis=&num)
                            %then %do;  missing  %end; ;

    data _d02;
      set _d02;
      %if (&bytype=num) %then %do; &by=99999999; %end;
                        %else %do; &by="zzzzzzz"; %end;
      total='y';

    data _final;
      set _final _d02;

    proc sort data=_final;
      by &v1 &by;

    data _d2;
      merge _d2 _level  _f1;
      by _merge;
/*      %if &pfoot=Y %then %do; length pvalue $ 17; %end; %else %do; length pvalue $ 7; %end; */
      format p1 p2 p3 p4 p5 p6 p7 fpvalue.;
      label p1='Chi-Square p-value'
            p2="Fisher's Exact p-value"
            p3='Kruskal-Wallis p-value'
            p4="Exact Kruskal-Wallis p-value"
            p5='Wilcoxon p-value'
            p6='Exact Wilcoxon p-value'
            p7='ANOVA F-test p-value';
      p1=round(p_pchi,&ppct.); p2=round(xp2_fish,&ppct.);
      p3=round(p_kw,&ppct.); p4=round(xp_kw,&ppct.);
      p5=round(p2_wil,&ppct.); p6=round(xp2_wil,&ppct.);
      p7=round(p_f,&ppct.);
      %if &pval>0 /*(1<=&pval<=7)*/ %then %do;
        if (p&pval.<&ppct. & p&pval.^=.) then do;
           %if &pfoot=Y %then %do; pvalue="<&ppct.~{super &&pnum&pval}"; %end;
           %else %do; pvalue="<&ppct."; %end;
        end;
        else do;
           %if &pfoot=Y %then %do; pvalue=trim(left(put(round(p&pval.,&ppct.),6.&pdec.))) || "~{super &&pnum&pval}"; %end;
           %else %do; pvalue=put(round(p&pval.,&ppct.),6.&pdec); %end;
        end;
      %end;
      %else %do; pvalue=' '; %end;
      keep l1-l&nby. n1-n&nby. nall pvalue p1-p7;
    run;

    data _final;
      set _d2 _final;
      if (count=.) then count=0;
      %if (&bytype=num) & (&v1type=num) %then %do;
          if (percent=.) & (&v1^=.) & (&by^=.) then percent=0;
          if (pct_col=.) & (&v1^=.) & (&by^=.) then pct_col=0;
          if (pct_row=.) & (&v1^=.) & (&by^=.) then pct_row=0;
      %end;
      %if (&bytype=char) & (&v1type=num) %then %do;
          if (percent=.) & (&v1^=.) & (&by^=' ') then percent=0;
          if (pct_col=.) & (&v1^=.) & (&by^=' ') then pct_col=0;
      if (pct_row=.) & (&v1^=.) & (&by^=' ') then pct_row=0;
      %end;
      %if (&bytype=num) & (&v1type=char) %then %do;
          if (percent=.) & (&v1^=' ') & (&by^=.) then percent=0;
          if (pct_col=.) & (&v1^=' ') & (&by^=.) then pct_col=0;
          if (pct_row=.) & (&v1^=' ') & (&by^=.) then pct_row=0;
      %end;
      %if (&bytype=char) & (&v1type=char) %then %do;
          if (percent=.) & (&v1^=' ') & (&by^=' ') then percent=0;
          if (pct_col=.) & (&v1^=' ') & (&by^=' ') then pct_col=0;
          if (pct_row=.) & (&v1^=' ') & (&by^=' ') then pct_row=0;
     %end;

      length value vlbl vtype var $ 50;
      nvar=.; ny=.;
      nvar=&num;
      %if (&v1f=)  & (("&v1type"="num") | ("&v1type"="char"))
        %then %do; value=trim(&v1); %end;
      %if (&v1f^=) & (("&v1type"="num") | ("&v1type"="char"))
        %then %do; value=trim(put(&v1,&v1f)); %end;
      vlbl="&v1lbl";
      vtype="&v1type";
      ny=&ny;
      var="&v1";
      %if (&byf^=) %then %do; format &by &byf; %end;

    /* creates a file that can be easily printed with proc print */
    data _tmp (keep=nvar _line level ci c1-c&nby. p1-p7
                    ctotal pvalue);
      length level $ 50;
      length ci c1-c&nby. ctotal $ 40;
 /*     %if &pfoot=Y %then %do; length pvalue $ 17; %end; %else %do; length pvalue $7; %end; */
      set _final end=last;
      %if (&by1=) %then %do; label c1='Missing'; %end;
                  %else %do; label c1="&fby1";    %end;
      label level='#'
            ci='     '
            ctotal="Total";
      %if (&nby>1) %then %do; 
        %do i = 2 %to &nby;
          label c&i.="&&fby&i.";
        %end;
      %end;

      /* creates the first line as the label and p-values */
      if (_n_=1) then do;
         level=vlbl; ci='      '; cctotal=' ';
         %if (&nby>1) %then %do; %do i = 2 %to &nby.;  c&i.=' '; %end; %end;
         _line=1;
         pci1=.; pci2=.; nci1=0; nci2=0;
         output;
      end;
      else do;
         pvalue=' '; p1=.; p2=.; p3=.; p4=.; p5=.; p6=.; p7=.;
         if (&v1=' ') then level='    Missing';
                      else level="    " || trim(left(value));
         %if ((&by=' ') | (&v1=' ')) & ("&incmiss"^="Y")
         %then %do i = 1 %to &nby.;
            if (&by="&&by&i.") then c&i.=trim(left(count));
         %end;
         %else %do i = 1 %to &nby;
           %if &ind=5 | &dnp=2 %then %do;
             if (&by="&&by&i.") then c&i.=trim(left(count));
           %end;
           %else %if &ind=6 | &dnp=3 %then %do;
             if (&by="&&by&i.") then c&i.=trim(left(put(round(&pct_type,&decpct),&decpctf))) || "%";
           %end;
           %else %if &PCTPRSNT=1 %then %do;
             if (&by="&&by&i.") then c&i.=trim(left(count)) || ", " || trim(left(put(round(&pct_type,&decpct),&decpctf))) || "%";
           %end;
           %else %if &PCTPRSNT=2 %then %do;
             if (&by="&&by&i.") then c&i.=trim(left(count)) || " (" || trim(left(put(round(&pct_type,&decpct),&decpctf))) || "%)";
          %end;
           %if (&nby>1) %then %do; 
             if (&by="&by1") then do; pci1=&pct_type; nci1=nci1+count; end;
             if (&by="&by2") then do; pci2=&pct_type; nci2=nci2+count; end;
           %end;
         %end;
      end;
      if (total='y') then do;
        %do i = 1 %to &nby;
          %if &PCTPRSNT=1 %then %do;
          if (c&i.=' ') then c&i.= "0" || ", " || trim(left(put(round(0,&decpct),&decpctf))) || "%";
          /* different formats if the values are missing */
          if (&v1=' ') & ("&incmiss"^='Y') then ctotal=trim(left(count)); else do;
            %if &ind=5 | &dnp=2 %then %do; ctotal=trim(left(count)); %end;
            %else %if &ind=6 | &dnp=3 %then %do; ctotal=trim(left(put(round(percent,&decpct),&decpctf))) || "%"; %end;
            %else %if "&pct_type"="PCT_COL" %then %do; ctotal=trim(left(count)); %end;
            %else %do; ctotal=trim(left(count)) || ", " || trim(left(put(round(percent,&decpct),&decpctf))) || "%"; %end;
          end;
          %end;
          %if &PCTPRSNT=2 %then %do;
          if (c&i.=' ') and &PCTPRSNT=2 then c&i.= "0" || " (" || trim(left(put(round(0,&decpct),&decpctf))) || "%)";
          /* different formats if the values are missing */
          if (&v1=' ') & ("&incmiss"^='Y') then ctotal=trim(left(count)); else do; 
            %if &ind=5 | &dnp=2 %then %do; ctotal=trim(left(count)); %end;
            %else %if &ind=6 | &dnp=3 %then %do; ctotal=trim(left(put(round(percent,&decpct),&decpctf))) || "%"; %end;
            %else %if "&pct_type"="PCT_COL" %then %do; ctotal=trim(left(count)); %end;
            %else %do; ctotal=trim(left(count)) || " (" || trim(left(put(round(percent,&decpct),&decpctf))) || "%)"; %end;
          end;
          %end;
        %end;
        _line=_line+1;
        output;
      end;

      /* inserts a 95% confidence interval if requested */
      if (last) & ("&cip"="2") then do;
         mdiff=pci1-pci2; sd=sqrt(pci1*(100-pci1)/nci1 + pci2*(100-pci2)/nci2);
         lowerci=round(mdiff-1.96*sd,0.1);
         upperci=round(mdiff+1.96*sd,0.1);
         mdiff=round(mdiff,0.1);
         _line=_line+1;
         level=' ';
         %do i = 1 %to &nby; c&i.=' '; %end;
         ctotal=' ';
         _line=_line+1;
         output;
         level='Difference (95% CI)';
         ci=trim(left(mdiff)) || ' (' || trim(left(lowerci)) || ', ' || trim(left(upperci)) || ')';
         output;
      end;

      /* inserts a blank line after each variable */
      if (last) then do;
         level=' ';
         %do i = 1 %to &nby.; c&i.=' '; %end;
         ctotal=' '; pvalue=' '; ci=' ';
         _line=_line+1;
         output;
      end;
      retain ci c1-c&nby. ctotal _line pci1 pci2 nci1 nci2 p1-p4;
    run;

    /* removes lines if only parts of variables are printed */
    data _tmp;
      set _tmp end=last;
       if ("&ind"="p1") | ("&ind"="n1") then do;
          if (_line>2) & not last then delete;
       end;
       if ("&ind"="p2") | ("&ind"="n2") then do;
          if (_line=2) then delete;
       end;
       if ("&ind"="p3") then do;
          if (_line=2) | (_line=3) then delete;
       end;

    /* appends this variable to the main analysis file */
    %if (&num=1) %then %do;
        data _mst;
          set _tmp;
    %end;
    %else %do;
        data _mst;
          set _mst _tmp;
    %end;
  %end;  ** end of discrete/ordinal analysis **;

  %if (&ind=1) | (&ind=7) | (&ind=8) | (&ind=9) %then %do;

    /** analyses for continuous variables **/
    /* overall distributions without missing values */
    proc univariate data=_master noprint;
      var &v1;
      where (_by_ ^= ' ');
      output out=_e0 mean=mean median=median std=std min=min max=max n=n nmiss=nmiss q1=q1 q3=q3;

    data _e0;
      set _e0;
      total='y';

    /* distributions by group (including missing as a group) */
    proc sort data=_master;
      by &by;

    proc univariate data=_master noprint;
      var &v1;
      by &by;
      output out=_e mean=mean median=median std=std min=min max=max n=n nmiss=nmiss var=var q1=q1 q3=q3;

    data _all;
      set _e _e0;
      std=round(std,&decstd);
      mean=round(mean,&decmean);
      median=round(median,&decmedian);
      format mean std median 10.1;
      
    /* p-values without missing values */
    data _f;
/*       %if &pfoot=Y %then %do; length pvalue $ 17; %end; %else %do; length pvalue $ 7; %end; */
      _merge=1; &by=.;
      pvalue=' '; p1=.; p2=.; p3=.; p4=.; p5=.; p6=.; p7=.;
    run;

    %if &nby_miss>1 & "&pvalues"="Y" %then %do;
      proc npar1way data=_master wilcoxon anova noprint;
        var &v1;
        class &by;
        where _by_ ^= ' ';
        %if %eval(&nobs<=20) %then %do; exact wilcoxon; %end;
        output out=_f wilcoxon anova;

      data _f;
        set _f;
  /*      %if &pfoot=Y %then %do; length pvalue $ 17; %end; %else %do; length pvalue $ 7; %end; */
        _merge=1;
        p1=.; p2=.;
        p3=round(p_kw,&ppct.);
        p4=round(xp_kw,&ppct.);
        p5=round(p2_wil,&ppct.);
        p6=round(xp2_wil,&ppct.);
        p7=round(p_f,&ppct.);
        %if (&pval>0) %then %do;
          if (p&pval.<&ppct. & p&pval.^=.) then do;
            %if &pfoot=Y %then %do; pvalue="<&ppct.~{super &&pnum&pval}"; %end; %else %do; pvalue="<&ppct."; %end; end;
          else do; %if &pfoot=Y %then %do; pvalue=trim(left(put(round(p&pval.,&ppct.),6.&pdec.))) || "~{super &&pnum&pval}"; %end;
                   %else %do; pvalue=put(round(p&pval.,&ppct.),6.&pdec.); %end; end;
        %end;
        %else pvalue=' ';
      run;
    %end;

    /* creates a file that can be easily printed with proc print */

     data _t;
      length level $ 50;
      length ci c1-c&nby. ctotal $ 40;
      _merge=1;
      level="&v1lbl"; _line=1;
      ci='      '; ctotal=' ';
      %if &nby>1 %then %do; %do i = 2 %to &nby.;  c&i.=' '; %end;  %end;
      output;

    data _p (keep=_line level ci c1-c&nby. ctotal pvalue p1-p7);
      merge _t _f; by _merge;

    %let ln=2;
  
    data _n (keep=_line level c1-c&nby. ctotal);
      set _all end=last;
      length level $ 50;
      length c1-c&nby. ctotal $ 40;
      retain c1-c&nby. ctotal;
      level='    N';  _line=&ln;
      %if (&cn=1 & &ind=1) | &ind=7 | &ind=8 %then %do; %let ln=%eval(&ln+1); %end;
      %if &cip=1 %then %do;
        if (&by="&by1") then do; call symput("cin1",trim(left(N))); end;
        if (&by="&by2") then do; call symput("cin2",trim(left(N))); end;
      %end;
      %do i = 1 %to &nby.;
        if (&by="&&by&i.") & total^='y' then c&i.=trim(left(N));
      %end;
      if (last) then ctotal=trim(left(N));
      if (last) then output;
      
    data _nmiss (keep=_line level c1-c&nby. ctotal);
      set _all end=last;
      length level $ 50;
      length c1-c&nby. ctotal $ 40;
      retain c1-c&nby. ctotal;
      level='    Missing';  _line=&ln;
      %if (&cnmiss=1 & &ind=1) %then %do; %let ln=%eval(&ln+1); %end;
      %do i = 1 %to &nby.;
        if (&by="&&by&i.")  & total^='y' then c&i.=trim(left(Nmiss));
      %end;
      if (last) then ctotal=trim(left(N));
      if (last) then output;

    data _mean (keep=_line level c1-c&nby. ctotal);
      set _all end=last;
      length level $ 50;
      length c1-c&nby. ctotal $ 40;
      retain c1-c&nby. ctotal;
      if (&cms=2 & &ind=1) then level='    Mean';
      else if (&cms=3 & &ind=1) then level='    SD';
      else level='    Mean (SD)'; _line=&ln; 
      %if (&cms>0 & &ind=1) | &ind=8 %then %do; %let ln=%eval(&ln+1); %end;
      %if &cip=1 %then %do;
        if (&by="&by1") then do; 
          call symput("cimean1",trim(left(MEAN))); 
          call symput("cisd1",trim(left(STD))); 
        end;
        if(&by="&by2") then do;  
          call symput("cimean2",trim(left(MEAN))); 
          call symput("cisd2",trim(left(STD))); 
        end;
      %end;
      %do i = 1 %to &nby.;
        if (&by="&&by&i.")  & total^='y' then do;
           %if &cms=2 & &ind=1  %then %do; c&i.=trim(left(put(MEAN,&decmeanf))); %end; 
           %else %if &cms=3 & &ind=1 %then %do; c&i.=trim(left(put(STD,&decstdf))); %end; 
           %else %if &MEANPRSNT=1 %then %do;
             c&i.=trim(left(put(MEAN,&decmeanf))) || ', ' || trim(left(put(STD,&decstdf)));
           %end;
           %else %if &MEANPRSNT=2 %then %do;
             c&i.=trim(left(put(MEAN,&decmeanf))) || ' (' || trim(left(put(STD,&decstdf))) || ')';
           %end;
        end;
      %end;
      if (last) then do;
           %if &cms=2 & &ind=1 %then %do; ctotal=trim(left(put(MEAN,&decmeanf))); %end; 
           %else %if &cms=3 & &ind=1 %then %do; ctotal=trim(left(put(STD,&decstdf))); %end;  
           %else %if &MEANPRSNT=1 %then %do;
             ctotal=trim(left(put(MEAN,&decmeanf))) || ', ' || trim(left(put(STD,&decstdf)));
           %end;
           %else %if &MEANPRSNT=2 %then %do;
             ctotal=trim(left(put(MEAN,&decmeanf))) || ' (' || trim(left(put(STD,&decstdf))) || ')';
           %end;
      end;
      if (last) then output;

    data _med (keep=_line level c1-c&nby. ctotal);
      set _all end=last;
      length level $ 50;
      length c1-c&nby. ctotal $ 40;
      retain c1-c&nby. ctotal;
      level='    Median';  _line=&ln;  
      %if (&cmd=1 & &ind=1) | &ind=7 %then %do; %let ln=%eval(&ln+1); %end;
      %do i = 1 %to &nby.;
        if (&by="&&by&i.") & total^='y'  then c&i.=trim(left(put(MEDIAN,&decmedianf)));
      %end;
      if (last) then ctotal=trim(left(put(MEDIAN,&decmedianf)));
      if (last) then output;

    data _quart (keep=_line level c1-c&nby. ctotal);
      set _all end=last;
      length level $ 50;
      length c1-c&nby. ctotal $ 40;
      retain c1-c&nby. ctotal;
      level='    Q1, Q3'; _line=&ln;
      %if (&cq=1 & &ind=1) | &ind=7 %then %do; %let ln=%eval(&ln+1); %end;
      %do i = 1 %to &nby.;
        if (&by="&&by&i.")  & total^='y' then c&i.=trim(left(put(Q1,&decq1f))) || ', ' || trim(left(put(Q3,&decq3f)));
      %end;
      if (last) then ctotal=trim(left(put(Q1,&decq1f))) || ', ' || trim(left(put(Q3,&decq3f)));
      if (last) then output;

    data _range (keep=_line level c1-c&nby. ctotal);
      set _all end=last;
      length level $ 50;
      length c1-c&nby. ctotal $ 40;
      retain c1-c&nby. ctotal;
      level='    Range'; _line=&ln; 
      %if (&cr=1 & &ind=1) | &ind=8 %then %do; %let ln=%eval(&ln+1); %end;
      %do i = 1 %to &nby.;
        if (&by="&&by&i.")  & total^='y' then c&i.= trim(left(put(MIN,&decminf))) || ', ' || trim(left(put(MAX,&decmaxf)));
      %end;
      if (last) then ctotal=trim(left(put(MIN,&decminf))) || ', ' || trim(left(put(MAX,&decmaxf)));
      if (last) then output;

    data _blank (keep=_line level ci c1-c&nby. ctotal pvalue p1-p7);
      length level $ 50;
      length ci c1-c&nby. ctotal $ 40;
      %do i = 1 %to &nby.;  c&i.=' '; %end;
      level=' '; ci='      '; ctotal=' ';
      pvalue=' '; p1=.; p2=.;  p3=.; p4=.;  p5=.; p6=.; p7=.;
      _line=&ln;

    data _final;
      set                                                  _p 
        %if (&cn=1 & &ind=1) | &ind=7 | &ind=8 %then %do;  _n     %end;
        %if (&cnmiss=1) %then %do;                         _nmiss %end;
        %if (&cmn=1 & &ind=1) | &ind=8 %then %do;          _mean  %end;
        %if (&cmd=1 & &ind=1) | &ind=7 %then %do;          _med   %end;
        %if (&cq=1 & &ind=1) | &ind=7 %then %do;           _quart %end;
        %if (&cr=1 & &ind=1) | &ind=8 %then %do;           _range %end; 
                                                           _blank ;
      nvar=&num;
    run;

   %if &cip=1 %then %do;
     data _ci (keep=_line nvar level c1-c&nby. ci ctotal p1-p7 );
       length level $ 50;
       length c1 ci c2 c3 c4 c5 c6 c7 c8 ctotal $ 40;
       nvar=&num;
       _line=7;
       mdiff=&cimean2-&cimean1;
       se=sqrt(&cisd1*&cisd1/&cin1 + &cisd2*&cisd2/&cin2);
       lowerci=mdiff-1.96*se;
       upperci=mdiff+1.96*se;
       mdiff=round(mdiff,0.1);
       lowerci=round(lowerci,0.1);
       upperci=round(upperci,0.1);
       level="Difference (95% CI)"; 
       ci=trim(left(mdiff)) || ' (' || trim(left(lowerci)) || ', ' || trim(left(upperci)) || ')'; 
       output;
       %do i = 1 %to &nby.;  c&i.=' '; %end;
       level=' '; ci='      '; ctotal=' ';
       pvalue=' '; p1=.; p2=.;  p3=.; p4=.;  p5=.; p6=.; p7=.;
       _line=8; output;

     data _final;
       set _final _ci;
     run;
   %end;
      
    %if (&num=1) %then %do;
        data _mst;
        set _final;
    %end;
    %else %do;
       data _mst;
       set _mst _final;
    %end;
  %end;  ** end of continuous analysis **;

    /** analyses for survival variables **/
  %if (&ind=4) %then %do;
     %if &snum=0 %then %let event=&surv; 
     %else %do; %if (%scan(&surv,&snum)^=) %then %do; %let event=%scan(&surv,&snum); %let snum=%eval(&snum+1); %end; %end;
     %if &scen=0 %then %let cen_vl=&scensor;
     %else %do; %if (%scan(&scensor, %scen)^=) %then %do; %let cen_vl=%scan(&scensor,&scen); %let scen=%eval(&scen+1);
     %end; %end;

     %let time=&v1;
     %let errorflg = 0;

     %if &time=  %then %do;
       %put  ERROR - Variable <time> not defined;
       %let  errorflg = 1;
       %end;

     %if &event=  %then %do;
       %put  ERROR - Variable <event> not defined;
       %let  errorflg = 1;
       %end;

     data _tmp_;  set _master;
       keep &by &time &event;
       where &by is  not missing;
       if &time=. or &time < 0 then do;
          error "ERROR - &time= " &time ' - not used.';
          &time = .;
          &event = .;
          end;
       if &event > &cen_vl+1 or &event < &cen_vl then do;
          error "ERROR - &event= " &event ' - not used.';
          &time = .;
          &event = .;
          end;

     proc sort; by &by &time;
     proc means noprint data=_tmp_;
       var &time;
       by &by;
       output out=_counts_ n=nrisk max=maxtime nmiss=tl_miss;

     data _sumry_ (keep=&by total cum_ev cum_cen median tl_miss);
       set _tmp_ nobs=nobs; by &by &time;
       retain pt nevent _kt_ ncensor nrisk cum_ev cum_cen
              total median firstmed;
       label cum_ev = "Cumulative events including (t)"
             cum_cen = "Cumulative censors including (t)"
             median = "Median Survival"
             tl_miss = "Total Missing";

       _ft_=first.&time;
       _lt_=last.&time;

              *do if the first observation per by group;
       if first.&by=0 then go to notfirst;

       set _counts_;
       &time=0; nevent=0;
       _kt_=0; ncensor=0; cum_ev=0; cum_cen=0; pt=1;
       years=0; total=nrisk; median=.; firstmed=.;

               *do for each observation in the dataset;
       NOTFIRST:

       if _ft_ then do;  *do for the first obs. per time;
         nevent=0;
         _kt_=0;
         ncensor=0;
       end;

                         *do for each observation;
       if &time ^= . then do;
         if &event = &cen_vl+1
         then do;
           nevent=nevent+1;
         end;
         else do;
           ncensor=ncensor+1;
         end;
         _kt_=_kt_+1;
       end;

       if _lt_ then do;  *do for the last observation per time;
         if _kt_=0 then go to next3;
         if nrisk>0 then pt=pt*(1-nevent/nrisk); else pt=.;
         cum_ev = cum_ev + nevent;
         cum_cen = cum_cen + ncensor;
         nrisk=nrisk-_kt_;
         if _kt_ = 0 then go to next3;
         if ABS(pt-0.5)<=0.00001 then do;
           if firstmed = . then firstmed = &time;
         end;
         if median=. and round(pt,0.00001) < 0.5  then do;
           if firstmed ^=.
           then median = (&time + firstmed)/2.0;
           else median = &time;
         end;
                             *output summary data;
         next3:

         if last.&by=1 then output _sumry_;
       end;
     run;

     ******* outputs totals for median, etc. **********;

     proc sort data=_tmp_; by &time;
     proc means noprint data=_tmp_;
       var &time;
       output out=_counts_ n=nrisk max=maxtime nmiss=tl_miss;

     data _tot_  (keep=tottot totcumev totcen totmed tl_miss);
       set _tmp_ nobs=nobs end=end; by &time;
       retain pt nevent _kt_ ncensor nrisk totcumev totcen
              tottot totmed firstmed;
       label totcumev = "Cumulative events including (t)"
             totcen = "Cumulative censors including (t)"
             totmed = "Median Survival"
             tl_miss = "Total Missing";

       _ft_=first.&time; _lt_=last.&time;

       if _n_^=1 then go to notfirst;
       set _counts_;
       &time=0;
       nevent=0; _kt_=0; ncensor=0; totcumev=0; totcen=0;
       pt=1; years=0; tottot=nrisk; totmed=.; firstmed=.;

       NOTFIRST:
       if _ft_ then do;  *do for the first obs. per time;
         nevent=0;
         _kt_=0;
         ncensor=0;
       end;
       if &time ^= . then do;
         if &event = &cen_vl+1
         then do;
           nevent=nevent+1;
         end;
         else do;
           ncensor=ncensor+1;
         end;
         _kt_=_kt_+1;
       end;
       if _lt_ then do;  *do for the last observation per time;
         if _kt_=0 then go to next3;
         if nrisk>0 then pt=pt*(1-nevent/nrisk); else pt=.;
         totcumev = totcumev + nevent;
         totcen = totcen + ncensor;
         nrisk=nrisk-_kt_;
         if _kt_ = 0 then go to next3;
         if ABS(pt-0.5)<=0.00001 then do;
           if firstmed = . then firstmed = &time;
         end;
         if totmed=. and round(pt,0.00001) < 0.5  then do;
           if firstmed ^=.
           then totmed = (&time + firstmed)/2.0;
           else totmed = &time;
         end;
         next3:
         if end then output _tot_;
       end;
     run;

     proc sort data=_tmp_; by &by; run;

     data _sumry_;
       set _sumry_ end=last; by &by;
       line_num=put(_n_,3.);

     data _tmp_;
       merge _sumry_ (in=in1) _tmp_(in=in2);  by &by;
       keep &by &time &event line_num;
       if in1 and in2;

     %survlrk(data=_tmp_,time=&time,death=&event,
                censor=&cen_vl,strata=line_num,out=_x1);

     data _sumry_;
       merge _x1 (in=in1) _sumry_(in=in2);  by line_num;
       format observed expected o_e 8.1 rr 5.3
              chisq 8.2 pvalue 8.4;
       drop line_num;
       if in1 and in2;
       if chisq = 0 then pvalue=.N;
     run;

    proc transpose data=_sumry_ out=_tmp_;
      var total cum_ev median;
    run;

    proc transpose data=_tot_ out=_tmp1_ prefix=val;
      var tottot totcumev totmed;
    run;

    data _t;
     length level $ 50;
     length ci c1-c&nby. ctotal $ 40;
     _merge=1;
     level="&v1lbl"; _line=1;
     ci='      '; ctotal=' ';
     %do i = 1 %to &nby.;  c&i.=' '; %end;
     output;


    data _p (keep=pval _merge);
      set _sumry_;
      _merge=1;
      if pvalue^=. then do;
        if pvalue<&ppct. then do;
            %if &pfoot=Y %then %do; pval="<&ppct.~{super &pnum8}"; %end; %else %do; pval="<&ppct."; %end;
        end; else do; %if &pfoot=Y %then %do; pval=put(round(pvalue,&ppct.),6.&pdec.) || "~{super &pnum8}"; %end;
                   %else %do; pval=put(round(pvalue,&ppct.),6.&pdec.); %end; 
        end;
        output;
      end;
   
    data _tot (keep=_line level ci c1-c&nby. ctotal pvalue _merge);
      merge _t _p; by _merge;
      pvalue=pval;

    data _tmp_; set _tmp_; _merge=1; run;
    data _tmp1_; set _tmp1_; _merge=1; run;
    
    data _all (keep=_line level c1-c&nby. ctotal _merge);
      merge _tmp_ _tmp1_; by _merge;
      length level $ 50;
      length c1-c&nby. ctotal $ 40;
      retain c1-c&nby. ctotal;
      if _n_=1 then do;
        level='    Number of Patients'; _line=2;
        tot=0;
        %do i = 1 %to &nby.;
          %let ii=%eval(&i+1);
          %if (&by1=) %then %do; c&ii.=trim(left(col&i)); %end;
          %else %do; c&i.=trim(left(col&i)); %end;
        %end;
        ctotal=trim(left(val1));
        output;
      end;
      if _n_=2 then do;
        level='    Number of Events'; _line=3;
        tot=0;
        %do i= 1 %to &nby.;
          %let ii=%eval(&i+1);
          %if (&by1=) %then %do; c&ii.=trim(left(col&i)); %end;
          %else %do; c&i.=trim(left(col&i)); %end;
        %end;
        ctotal=trim(left(val1));
        output;
      end;
      if _n_=3 then do;
        level='    Median Survival Time'; _line=4;
        %do i= 1 %to &nby.;
          %let ii=%eval(&i+1);
          %if (&by1=) %then %do; c&ii.=trim(left(col&i)); %end;
          %else %do; c&i.=trim(left(col&i)); %end;
        %end;
        ctotal=trim(left(val1));
        output;
      end;

    data _blank (keep=_line level ci c1-c&nby.
                      ctotal pvalue p1-p7 _merge);
      length level $ 50;
      length ci c1-c&nby. ctotal $ 40;
      _merge=1;
      %do i = 1 %to &nby.;  c&i.=' '; %end;
      level=' '; ci='      '; ctotal=' ';
      pvalue=' '; p1=.; p2=.;  p3=.; p4=.;  p5=.; p6=.; p7=.;
      _line=6;

    data _final;
      set _tot _all _blank;
      by _merge;
      nvar=&num;

    %if (&num=1) %then %do;
        data _mst;
        set _final;
    %end;
    %else %do;
       data _mst;
       set _mst _final;
    %end;

  %end; ** end of survival analysis **;


  %if (&debug=Y) %then %do;
    %put _all_;

    proc print data=_final;
      title '_final';

    proc print data=_mst;
      title '_mst';
  %end;

  %end; *end of errorvar - skips analysis if the variable does not exist;

  /* reads in the next variable */
  %let num=%eval(&num+1);
  %let v1=%upcase(%scan(&var,&num));
  %let indold=&ind;
  %let ind=%scan(&type,&num);
  %if (&ind=) %then %do; %let ind=&indold; %end;
  %if "&ptype"^="" %then %let pval=%scan(&ptype,&num);                  /* ptype=type of pvalue desired */
  %if "&ptype"="" %then %if &ind=1 | &ind=7 | &ind=8 %then %let pval=3;
  %else %if &ind=3 %then %let pval=5;
  %else %let pval=1;
  

  %if (&v1^=) & ( &pval>0 & &pfoot=Y) %then %do;
  %if &&pt&pval=0 %then %do;
    %let pt&pval=1;
    %let pfnum=%eval(&pfnum+1);
    %let pnum&pval=&pfnum;
    %let pft=&pft ~{super &pfnum}&&pnm&pval;
  %end; %end;
     
  
%end;

%let num=%eval(&num-1);


/* removes the final blank line */
data _mst;
  set _mst end=last;
  if (last) then delete;


/* adds in any comments for the table */
/* and deletes specified lines        */

data _mst;
  set _mst %if (&comments^=) %then %do; &comments %end; ;
  %do i = 1 %to 50;
    %let vr=%scan(&dvar,&i);
    %let dl=%scan(&dline,&i);
    %if (&vr^=) %then %do;
      if (&vr=nvar) & (&dl=_line) then delete;
    %end;
    %else %do; %let i = 50; %end;
  %end;

proc sort data=_mst;
  by nvar _line;

/* creates the table template */
proc template;
  define style newtable;
  style cellcontents /
     nobreakspace=on
     font_face="&bodyfnt."
     font_weight=medium
     font_style=roman
     font_size=&bodysz. pt
     just=center
     vjust=center
     asis=on
     font_size=1;
  style lhead /
     nobreakspace=on
     font_face="&headfnt."
     font_weight=bold
     font_size=&bodysz. pt
     font_style=roman
     just=center
     vjust=center
     asis=on
     font_size=1;
  style table /
     frame=&frame
     asis=on
     cellpadding=&space.
     cellspacing=&space.
     just=center
     rules=&rules
     borderwidth=2;
  style Body /
     font_face="&headfnt."
     asis=on
     font_size=&bodysz. pt;
  style BodyDate /
     font_face="&headfnt."
     asis=on
     font_size=&bodysz. pt;
  style SysTitleAndFooterContainer /
     font_face="&headfnt."
     asis=on
     font_size=&bodysz. pt;
  style SystemFooter /
     font_face="&headfnt."
     asis=on
     font_size=&bodysz. pt;
  style data /
     font_face="&headfnt."
     font_size=&bodysz. pt;
  style SystemTitle /
     font_face="&headfnt."
     font_size=&bodysz. pt;
  style ByLine /
     font_face="&headfnt."
     asis=on
     font_size=&bodysz. pt;
  style Header /
     font_face="&headfnt."
     asis=on
     font_size=&bodysz. pt;
  end;
  run;


%if ("&outdoc"^="") %then %do;

ods listing close;

options orientation=&page.;

%if "&doctype"="0" %then %do; ods rtf file="&outdoc..doc" style=newtable; %end;
%else %do; ods rtf file="&outdoc." style=newtable; %end;


%if &pfoot=Y %then %do; ods escapechar='~'; %end;

title ' ';

%let titles=&ttitle1;
%if "&ttitle2"^="" %then %do; %let titles=&titles.#&ttitle2.; %end;
%if "&ttitle3"^="" %then %do; %let titles=&titles.#&ttitle3.; %end;
%if "&ttitle4"^="" %then %do; %let titles=&titles.#&ttitle4.; %end;
%if &date=Y %then %do; %let fdate=(report generated on " sysdate9 "); %end; %else %do; %let fdate=; %end;


proc template;
  define table summ;
  mvar sysdate9;
  column level c1 c2-c&nby. ci ctotal pvalue;
  header table_header_1;
  footer table_footer_1 table_footer_2;
  define table_header_1;
     text "&titles.# ";
     style=header{font_size=&titlesz. pt font_face="&titlefnt." %if &titlebld=Y %then %do; font_weight=bold %end; };
     split='#';
  end;
  define table_footer_1;
     %if (&date=Y) | ("&foot"^="") | &pfoot=Y %then %do; 
       %if &pfoot=Y %then %do; text " &fdate#&pft#&foot"; %end;
       %else %do; text " &fdate#&foot"; %end;
     %end; %else %do; text ""; %end;
     split='#';
     just=left;
     style=header{font_size=&footsz. pt font_face="&headfnt."};
  end;
  
  define table_footer_2;
     %if ("&foot2"^="") | ("&foot3"^="") | ("&foot4"^="") | ("&foot5"^="") %then %do; text " &foot2#&foot3#&foot4#&foot5"; %end; 
       %else %do; text ""; %end;
     split='#';
     just=left;
     style=header{font_size=&footsz. pt font_face="&headfnt."};
  end;

  define header header;
     split='#';
  end;
  define column level;
         generic=on;
         vjust=top;
         just=left;
         header=" ";
         cellstyle substr(_val_,1,1)^=' ' as
           cellcontents{font_weight=bold   font_size=&bodysz. pt font_face="&headfnt." 
           %if (&levelwd^=) %then %do; cellwidth=&levelwd. %end; },
         substr(_val_,1,1)=' ' as
           cellcontents{font_weight=medium font_size=&bodysz. pt font_face="&headfnt."
           %if (&levelwd^=) %then %do; cellwidth=&levelwd. %end; };
  end;

  %do i = 1 %to &nby.;
    define column c&i.;
         generic=on;
         style=cellcontents{font_size=&bodysz. pt font_face="&bodyfnt."
                            %if (&datawd^=) %then %do; cellwidth=&datawd. %end; };
         vjust=top;
         just=center;
         %if (&&by&i.=) %then %do;
           header="Missing                                     (N=&&byn&i.)";
         %end;
         %else %do;
           header="&&fby&i.                                        (N=&&byn&i.)";
         %end;
         end;
  %end;

  define column ci;
         generic=on;
         style=cellcontents{font_size=&bodysz. pt font_face="&bodyfnt."
                            %if (&ciwd^=) %then %do; cellwidth=&ciwd. %end; };
         vjust=top;
         just=center;
         header=" ";
         end;

  define column ctotal;
         generic=on;
         style=cellcontents{font_size=&bodysz. pt font_face="&bodyfnt."
                            %if (&datawd^=) %then %do; cellwidth=&datawd. %end; };
         vjust=top;
         just=center;
         header="Total                                       (N=&nobs.)";
         end;
  define column pvalue;
         generic=on;
         style=cellcontents{font_size=&bodysz. pt font_face="&bodyfnt."
                            %if (&pvalwd^=) %then %do; cellwidth=&pvalwd. %end; };
         vjust=top;
         just=center;
         header="p value";
         end;
  end;
  run;

  data _null_;
    set _mst;

    file print ods=(template='summ'
         columns=(level=level(generic=on)
                 %do i = 1 %to &nby.;
                    c&i.=c&i.(generic=on)
                 %end;
                 %if (&ci^=) %then %do;  ci=ci(generic=on) %end;
                 %if (&total=Y)   %then %do;  ctotal=ctotal(generic=on)  %end;
                 %if (&pvalues=Y) %then %do;  pvalue=pvalue(generic=on)  %end;
                 ));
    put _ods_;
    run;

  ods rtf close;
  ods listing;
  
  %if "&doctype"="0" %then %do; %put Created File: &outdoc..doc; %end;
  %else %do; %put Created File: &outdoc.; %end;

%end;

/* prints the final table */
options linesize=128 pagesize=56;
%if &print=Y %then %do;
  ods trace on;
  proc print data=_mst label noobs split='*';
    var %if (&debug=Y) %then %do;  nvar _line  %end;
        level c1-c&nby 
                 %if (&ci^=) %then %do;  ci %end;
                 %if (&total=Y)   %then %do;  ctotal  %end;
                 %if (&pvalues=Y) %then %do;  pvalue  %end;;
    title "&ttitle1";
    %if "&ttitle2"^="" %then %do; title2 "&ttitle2"; %end;
  run;
  ods trace off;

%end;

%if (&outdat^=) %then %do;
  data &outdat;
    length tit1-tit4 ft1-ft5 $100. ;
    set _mst;
    tit1=''; tit2=''; tit3=''; tit4=''; ft1=''; ft2=''; ft3=''; ft4=''; ft5='';
    if _n_=1 then do;
      oldline=_line;
      _line=0; 
      nby=&nby;
      %do i = 1 %to &nby.; byn&i.=&&byn&i; 
      %if &bytype=num %then %do; %if (&&by&i.^=) %then %do; by&i=&&by&i; %end; %end;
      %else %if &bytype=char %then %do; by&i="&&by&i"; %end; fby&i="&&fby&i"; %end;
      nobs=&nobs.; pfoot="&pfoot"; %if &pfoot=Y %then %do; pft="&pft"; %end;
      tit1="&ttitle1"; tit2="&ttitle2"; tit3="&ttitle3"; tit4="&ttitle4";
      ft1="&foot"; ft2="&foot2"; ft3="&foot3"; ft4="&foot4"; ft5="&foot5";
      %if (&date=Y) %then %do; fdate="Y"; %end; %else %do;  fdate="N"; %end;
      %if (&ci=) %then %do; haveci=0; %end;  %else %do; haveci=1; %end;  
      %if (&total=N)   %then %do; havetot=0; %end; %else %do; havetot=1; %end; 
      %if (&pvalues=N) %then %do; havep=0; %end; %else %do;havep=1; %end; 
      output;
      _line=oldline;
    end;
    output;
  run;
%end;

/* cleans up the temporary data sets */
proc datasets nolist;
    delete _c2_ _check _c_ _d01 _d02 _d1 _d2 _dall _f1 _n _f _e _e0
          _final _level _master _mst _tmp _tmp_ _sumry_ _x1
          _t _all _tot _p _mean _med _range _blank _tot_
          _counts_ _lr1 _lr2 _lrmstr _tmp1_ _nmiss _quart;
    run;
    quit;

run;


%end;
%else %do; %put ERROR: &errorwhy; %end;

options validvarname=&validvn &sdate &snotes &snumb;

%mend stat_table;






*******************************;
* Macro for TND baseline table ;
*******************************;

%macro TNDbaselinetable(indata, cohort, subgroup, outcome);
    proc univariate data=&indata;
    where cohort="&cohort" and variant_subgroup="&subgroup";
    var Frailty_index;
    ods output Quantiles=Qs;
    run;
    proc sql; 
    select estimate into: Q1 from Qs where quantile='25% Q1';
    select estimate into: Q2 from Qs where quantile='50% Median';
    select estimate into: Q3 from Qs where quantile='75% Q3';
    quit;
    data table; set &indata;
    if cohort="&cohort" and variant_subgroup="&subgroup";
    * Create frailty category based on top qualtile;
        if Frailty_index >=&Q3 then Frailty=4;
        else if Frailty_index >=&Q2 then Frailty=3;
        else if Frailty_index >=&Q1 then Frailty=2;
        else Frailty=1;
    * Adjust variable format and create new variables for stddiff macro;
        outcome=(COVID19='P');
        sex_n=(sex='Male');
        Specimen_Type_n=(Specimen_Type='Saliva');
        index_date_c=put(index_date, monyy7.);
    run;
* Get subcategories;
    * IC stuatus;
    data table_IC; set table; where IC; run;
    proc sql; select (count(*)>0) into: IC from table_IC; quit;
    * Autoimmune conditions;
    data table_AI; set table; where Autoimmune; run;
    proc sql; select (count(*)>0) into: AI from table_AI; quit;
    * Pregnancy trimester;
    data table_PT; set table; where Pregnancy; run;
    proc sql; select (count(*)>0) into: PT from table_PT; quit;
* Get means, frequencies and p-values;
    %stat_TABLE(DSN=table, ID=mrn, by=&outcome,
                var=age age_grp sex_n Race_eth BMI_grp smoking Charlson_wt, 
                type=1 2 2 2 2 2 2, pvalues=Y, 
                decmean=2, decstd=2, decmedian=0, decmin=0, decmax=0, decq1=0, decq3=0, PCTPRSNT=2, MEANPRSNT=2, outdat=Part1);
    %stat_TABLE(DSN=table, ID=mrn, by=&outcome,
                var=Frailty,
                type=2, pvalues=Y, 
                decmean=2, decstd=2, decmedian=2, decmin=2, decmax=2, decq1=2, decq3=2, PCTPRSNT=2, MEANPRSNT=2, outdat=Part12);
    %stat_TABLE(DSN=table, ID=mrn, by=&outcome,
                var=Kidney Heart Lung Liver Diabetes,
                type=2 2 2 2 2, pvalues=Y, PCTPRSNT=2, outdat=Part13);
    %stat_TABLE(DSN=table, ID=mrn, by=&outcome,
                var=IC,
                type=2, pvalues=Y, PCTPRSNT=2, outdat=Part14);
    %if &IC=1 %then %do;
    %stat_TABLE(DSN=table_IC, ID=mrn, by=&outcome,
                var=IC_HIV IC_DX IC_TX IC_MEDS,
                type=2 2 2 2, pvalues=N, PCTPRSNT=2, outdat=Part142);
    %end;
    %stat_TABLE(DSN=table, ID=mrn, by=&outcome,
                var=Autoimmune, 
                type=2, pvalues=Y, PCTPRSNT=2, outdat=Part2);
    %if &AI=1 %then %do;
    %stat_TABLE(DSN=table_AI, ID=mrn, by=&outcome,
                var=RA IBD PPA MS SLE, 
                type=2 2 2 2 2, pvalues=N, PCTPRSNT=2, outdat=Part22);
    %end;
    %stat_TABLE(DSN=table, ID=mrn, by=&outcome,
                var=Pregnancy, 
                type=2, pvalues=Y, PCTPRSNT=2, outdat=Part3);
    %if &PT=1 %then %do;
    %stat_TABLE(DSN=table_PT, ID=mrn, by=&outcome,
                var=trimester, 
                type=2 2 2, pvalues=N, PCTPRSNT=2, outdat=Part32);
    %end;
    %stat_TABLE(DSN=table, ID=mrn, by=&outcome,
                var=hist_covid_diag hist_covid_test VA_grp ED_grp IP_grp 
                    preventive_care Medicaid h_med_income KPSC_employee 
                    medical_center_and_area index_date_c Specimen_Type_n, 
                type=2 2 2 2 2  2 2 2 2  2 2 2, pvalues=Y, 
                PCTPRSNT=2, outdat=Part4);
** Combine and reorder the data to be consistent with the table contents;
* Structure each parts;
    data Part1; set Part1;
    * Create a total line at begining;
    if _N_=3 then do; order=1; level='Total'; c1='N='||c1; c2='N='||c2; ctotal='N='||ctotal;  end;
    else if level='    N' then delete;
    run;
    data Part13; set Part13(rename=(level=level_o c1=c1_o c2=c2_o ctotal=ctotal_o pvalue=pvalue_o));
    length level $100.;
    retain level pvalue;
    if _line=0 then do; level='Comorbidities'; output; end;
    if _line=1 then do level=level_o; pvalue=pvalue_o; end;
    if _line=3 then do c1=c1_o; c2=c2_o; ctotal=ctotal_o; output; end;
    run;
    %if &IC=1 %then %do;
    data Part142; set Part142(rename=(level=level_o));
    length level $100.;
    a=lag(lag(level_o)); b=lag(level_o);
    if a='' then level=b; /*for situation where row 0 does not exist*/
    else level=a; 
    if level_o='    1' then output;
    run;
    %end;
    %if &AI=1 %then %do;
    data Part22; set Part22(rename=(level=level_o));
    length level $100.;
    a=lag(lag(level_o)); b=lag(level_o);
    if a='' then level=b; /*for situation where row 0 does not exist*/
    else level=a; 
    if level_o='    1' then output;
    run;
    %end;
    %if &PT=1 %then %do;
    data Part32; set Part32;
    if level in ('    1st trimester','    2nd trimester','    3rd trimester') ;
    run;
    %end;
    * Combine and reorder;
    data table_counts(keep=Var negative positive total pvalue order); 
    length level $100.;
    set Part1 Part12 Part13 
        Part14 %if &IC=1 %then %do; Part142 %end;
        Part2  %if &AI=1 %then %do; Part22 %end;
        Part3  %if &PT=1 %then %do; Part32 %end;
        Part4;
    order=_N_;
    * Move total line to top line;
    if level='Total' then order=0;
    * Remove useless title lines;
    if _line=0 and level^='Comorbidities' then delete; 
    * Remove other lines that report total counts;
    if level='    N' then delete;
    * Remove empty lines;
    if level=' ' then delete;
    * Reorder race-ethnicity categories to match table shell;
    if level='    Non-Hispanic White' then order=order-3;
    if level='    Non-Hispanic Black' then order=order-1;
    if level='    Hispanic' then order=order+2;
    if level='    Non-Hispanic Asian' then order=order+2;
    * Reorder BMI categories to match table shell;
    if level='    <18.5' then order=order-5.5;
    * Reorder Smoking categories to match table shell;
    if level='    Unknown' and _line=3 then order=order+1.5;
    * Recorder household income to match table sheel;
    if level='    <$40,000' then order=order-3.5;
    rename level=var c1=negative c2=positive ctotal=total;
    if level in ('age_grp', 'sex_n', 'race_eth') then pvalue='N/A';
    run;
* Get absolute standardized difference;
    %stddiff(inds =table, groupvar = &outcome, 
             numvars  = age ,
             charvars = age_grp sex_n Race_eth BMI_grp smoking Charlson_wt Frailty
                        Kidney Heart Lung Liver Diabetes IC Autoimmune Pregnancy hist_covid_diag hist_covid_test 
                        VA_grp ED_grp IP_grp preventive_care Medicaid h_med_income KPSC_employee medical_center_and_area index_date_c Specimen_Type_n,
             stdfmt = 8.4, outds = table_std ); 

* Merge to get all information for Table;
    proc sql; create table temp as 
    select a.order, a.var, a.positive as pos_&subgroup., a.negative as neg_&subgroup., a.pvalue as p_&subgroup., 
           case when a.pvalue='N/A' then 'N/A'
                else compress(b.Stddiff,'-') 
           end as ABS_&subgroup.
      from table_counts a left join table_std b on upcase(a.var)=upcase(b.VarName)
     order by a.order;
    quit;
    data &subgroup; set temp;
    length Rowname $100. temp $100.;
    retain temp;
    if p_&subgroup.^='' then temp=var;
    if p_&subgroup.^='' or _N_=1 then Rowname=var;
    else Rowname=catt(temp,var);
    run;
    proc datasets lib=work;
    delete part: table table_AI table_IC table_PT table_std temp;
    run; 
%mend TNDbaselinetable;






*******************************;
* Macro for TND analysis ;
*******************************;

%macro TND(indata, variant, cohort, outcome, output, condition_model=1);
    proc univariate data=&indata;
    where cohort="&cohort" and variant_subgroup="&variant";
    var Frailty_index;
    ods output Quantiles=Qs;
    run;
    proc sql; 
    select estimate into: Q1 from Qs where quantile='25% Q1';
    select estimate into: Q2 from Qs where quantile='50% Median';
    select estimate into: Q3 from Qs where quantile='75% Q3';
    quit;
    data analysis; set &indata;
    where cohort="&cohort" and variant_subgroup="&variant";
    * Create frailty category based on top qualtile;
        if Frailty_index >=&Q3 then Frailty=4;
        else if Frailty_index >=&Q2 then Frailty=3;
        else if Frailty_index >=&Q1 then Frailty=2;
        else Frailty=1;
    * Define outcome;
        covid_dx=(COVID19='P');
    * Define variable used for ASD macro, so it can be called for model adjustment based on ASD value;
        sex_n=(sex='Male');
        Specimen_Type_n=(Specimen_Type='Saliva');
        index_date_c=put(index_date, monyy7.);
    run;
    * Get frequency counts;
    proc freq data=analysis;
    table &outcome.*ARM/nocol nopercent or;
    ods output CrossTabFreqs=counts;
    run;
    data vac unvac; set counts;
    where _type_='11';
    counts=cat(frequency, ' (', compress(put(rowpercent, 5.1),' '), '%)');
    if ARM="Vaccinated" then output vac;
    if ARM='Unvaccinated' then output unvac;
    run;
    * Get unadjusted OR;
    proc logistic data=analysis desc;
    class arm(ref='Unvaccinated')/param=reference;
    model &outcome.=arm / maxiter=100;
    %if &condition_model=1 %then %do; strata pair_id; %end;
    ods output OddsRatios=OddRatio;
    run;
    data OR; set OddRatio;
    if OddsRatioEst<=1 then VEest=round((1-OddsRatioEst)*100,0.1);
    if OddsRatioEst>1 then VEest=round((1/OddsRatioEst-1)*100,0.1);
    if LowerCL<=1 then VEUpperCL=round((1-LowerCL)*100,0.1);
    if LowerCL>1 then VEUpperCL=round((1/LowerCL-1)*100,0.1);
    if UpperCL<=1 then VELowerCL=round((1-UpperCL)*100,0.1);
    if UpperCL>1 then VELowerCL=round((1/UpperCL-1)*100,0.1);
    run;
    * Get adjusted OR;
    title "&outcome. - &cohort. &variant.";
    proc logistic data=analysis desc;
    class arm(ref='Unvaccinated') age_grp sex Race_eth medical_center_and_area(ref='05') &varlist/param=reference;
    model &outcome.=arm &varlist  &interaction %if &condition_model=0 %then %do; age_grp sex Race_eth %end; / maxiter=100;
    %if &condition_model=1 %then %do; strata pair_id; %end;
    ods output ParameterEstimates=PE OddsRatios=AdjOddRatio;
    run;
    title;
    data AdjOR; set AdjOddRatio;
    if OddsRatioEst<=1 then VEest=round((1-OddsRatioEst)*100,0.1);
    if OddsRatioEst>1 then VEest=round((1/OddsRatioEst-1)*100,0.1);
    if LowerCL<=1 then VEUpperCL=round((1-LowerCL)*100,0.1);
    if LowerCL>1 then VEUpperCL=round((1/LowerCL-1)*100,0.1);
    if UpperCL<=1 then VELowerCL=round((1-UpperCL)*100,0.1);
    if UpperCL>1 then VELowerCL=round((1/UpperCL-1)*100,0.1);
    run;
    ** Put all results together;
    * Get counts for analysis with no case available;
    proc sql; create table counts_&output as
    select "&output" as Variant,
           a.counts as Pvaccinated, b.counts as Punvaccinated, '' as separater1, 
           c.counts as Nvaccinated, d.counts as Nunvaccinated
      from vac(where=(&outcome=1)) a,   
           unvac(where=(&outcome=1)) b, 
           vac(where=(&outcome=0)) c,   
           unvac(where=(&outcome=0)) d;
    quit;    * Get all info;
    proc sql; create table results_&output as
    select "&output" as Variant,
           a.counts as Pvaccinated, b.counts as Punvaccinated, '' as separater1, 
           c.counts as Nvaccinated, d.counts as Nunvaccinated, '' as separater2,
           cat(put(e.OddsRatioEst, 5.3), ' (', put(e.LowerCL, 5.3), ', ', put(e.UpperCL,5.3), ')') as OR,
           cat(put(f.OddsRatioEst, 5.3), ' (', put(f.LowerCL, 5.3), ', ', put(f.UpperCL,5.3), ')') as adjOR, 
           '' as separater3,
           cat(strip(put(e.VEest,5.1)), '% (', strip(put(e.VELowerCL,5.1)), '%, ', strip(put(e.VEUpperCL,5.1)), '%)') as VE,
           cat(strip(put(f.VEest,5.1)), '% (', strip(put(f.VELowerCL,5.1)), '%, ', strip(put(f.VEUpperCL,5.1)), '%)') as adjVE,
           '' as separater4,
           f.VEest as estimate, 
           f.VELowerCL as lower_CI, 
           f.VEUpperCL as upper_CI
      from vac(where=(&outcome=1)) a,   
           unvac(where=(&outcome=1)) b, 
           vac(where=(&outcome=0)) c,   
           unvac(where=(&outcome=0)) d, 
           OR(where=(upcase(substr(effect,1,3))='ARM')) e,
           AdjOR(where=(upcase(substr(effect,1,3))='ARM')) f,
           PE(where=(upcase(variable)='ARM')) g;
    quit;
    proc datasets lib=work;
    delete vac unvac  AdjOR PE;
    run; 
%mend TND;





*******************************;
* Macro for TND subgroup analysis ;
*******************************;

%macro TNDsubgroup(cohort,variant,length);
    proc univariate data=subgroupflag2;
    where cohort="&cohort" and variant_subgroup="&variant";
    var Frailty_index;
    ods output Quantiles=Qs;
    run;
    proc sql; 
    select estimate into: Q1 from Qs where quantile='25% Q1';
    select estimate into: Q2 from Qs where quantile='50% Median';
    select estimate into: Q3 from Qs where quantile='75% Q3';
    quit;
    data analysis; set subgroupflag2;
    where cohort="&cohort" and variant_subgroup="&variant";
    * Create frailty category based on top qualtile;
        if Frailty_index >=&Q3 then Frailty=4;
        else if Frailty_index >=&Q2 then Frailty=3;
        else if Frailty_index >=&Q1 then Frailty=2;
        else Frailty=1;
    * Define outcome;
        covid_dx=(COVID19='P');
    * Define variable used for ASD macro, so it can be called for model adjustment based on ASD value;
        sex_n=(sex='Male');
        Specimen_Type_n=(Specimen_Type='Saliva');
        index_date_c=put(index_date, monyy7.);
    run;
    * Get frequency counts;
    proc sort data=analysis; by covid_dx;
    proc freq data=analysis;
    by covid_dx;
    table arm;
    ods output OneWayFreqs=counts;
    run;
    data vac unvac; set counts;
    counts=cat(frequency, ' (', compress(put(percent, 5.1),' '), '%)');
    if ARM^="Unvaccinated" then output vac;
    if ARM='Unvaccinated' then output unvac;
    run;
    * Get unadjusted OR;
    proc logistic data=analysis desc;
    class arm(ref='Unvaccinated')/param=reference;
    model covid_dx=arm / maxiter=100;
/*    strata pair_id;*/
    ods output OddsRatios=OddRatio;
    run;
    data OR; set OddRatio;
    if OddsRatioEst<=1 then VEest=round((1-OddsRatioEst)*100,0.1);
    if OddsRatioEst>1 then VEest=round((1/OddsRatioEst-1)*100,0.1);
    if LowerCL<=1 then VEUpperCL=round((1-LowerCL)*100,0.1);
    if LowerCL>1 then VEUpperCL=round((1/LowerCL-1)*100,0.1);
    if UpperCL<=1 then VELowerCL=round((1-UpperCL)*100,0.1);
    if UpperCL>1 then VELowerCL=round((1/UpperCL-1)*100,0.1);
    run;
    * Get adjusted OR;
    title "&cohort. &variant.";
    proc logistic data=analysis desc;
    class arm(ref='Unvaccinated') age_grp sex Race_eth medical_center_and_area(ref='05') &varlist /param=reference;
    model covid_dx=arm &varlist &interaction age_grp sex Race_eth/ maxiter=100;
/*    strata pair_id; */
    ods output ParameterEstimates=PE OddsRatios=AdjOddRatio;
    run;
    title;
    data AdjOR; set AdjOddRatio;
    if OddsRatioEst<=1 then VEest=round((1-OddsRatioEst)*100,0.1);
    if OddsRatioEst>1 then VEest=round((1/OddsRatioEst-1)*100,0.1);
    if LowerCL<=1 then VEUpperCL=round((1-LowerCL)*100,0.1);
    if LowerCL>1 then VEUpperCL=round((1/LowerCL-1)*100,0.1);
    if UpperCL<=1 then VELowerCL=round((1-UpperCL)*100,0.1);
    if UpperCL>1 then VELowerCL=round((1/UpperCL-1)*100,0.1);
    run;    
    * Put all results together;
    proc sql; create table unvac2 as
    select b.arm, a.covid_dx, a.counts
      from unvac a, vac b;
    quit;    
    proc sql; create table results_&variant as
    select distinct "&variant" as variant, a.arm,
           a.counts as Pvaccinated, b.counts as Punvaccinated, '' as separater1, 
           c.counts as Nvaccinated, d.counts as Nunvaccinated, '' as separater2,
           cat(put(e.OddsRatioEst, 5.3), ' (', put(e.LowerCL, 5.3), ', ', put(e.UpperCL,5.3), ')') as OR,
           cat(put(f.OddsRatioEst, 5.3), ' (', put(f.LowerCL, 5.3), ', ', put(f.UpperCL,5.3), ')') as adjOR, 
           '' as separater3,
           cat(strip(put(e.VEest,5.1)), '% (', strip(put(e.VELowerCL,5.1)), '%, ', strip(put(e.VEUpperCL,5.1)), '%)') as VE,
           cat(strip(put(f.VEest,5.1)), '% (', strip(put(f.VELowerCL,5.1)), '%, ', strip(put(f.VEUpperCL,5.1)), '%)') as adjVE,
           '' as separater4,
           f.VEest as estimate, 
           f.VELowerCL as lower_CI, 
           f.VEUpperCL as upper_CI
      from vac(where=(covid_dx=1)) a left join unvac2(where=(covid_dx=1)) b on a.arm=b.arm
                                     left join vac(where=(covid_dx=0)) c on a.arm=c.arm
                                     left join unvac2(where=(covid_dx=0)) d on a.arm=d.arm
                                     left join OR(where=(upcase(substr(effect,1,3))='ARM')) e on a.arm=strip(substr(e.effect,5,&length))
                                     left join AdjOR(where=(upcase(substr(effect,1,3))='ARM')) f on a.arm=strip(substr(f.effect,17,&length))
                                     left join PE(where=(variable='arm')) g on a.arm=g.classval0;
    quit;
%mend TNDsubgroup;



