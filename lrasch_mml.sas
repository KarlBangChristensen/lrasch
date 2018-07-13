/**************************************************************************

MML2d: a SAS macro that can be used to fit longitudinal polytomous Rasch 
models using marginal maximum likelihood

***************************************************************************

Created by 

Maja Olsbjerg 
Karl Bang Christensen
Department of Biostatistics
University of Copenhagen.

this version 1.0 published April 23, 2013

documentation is avalable from

http://biostat.ku.dk/~kach/macro/lrasch_mml_documentation.pdf

**************************************************************************


DATA: the data set containing the items (scored 0,1, .. ,'max'). 
	  We assume that the items at the two time points are the same. 
	  The number of response categories 'max' can differ between items (within time point).  
	  

ITEM_NAMES: data set that contains information about the items. This data set should contain the variables 

		   item_names1: item names of time 1 items
		   item_names2: item names of time 2 items
		   item_text (optional): item text for plots 
		   max1: maximum score for time 1 items.
		   max2: maximum score for time 2 items. 

ANCHOR: a list (separated by blanks) of items (given by their time 1 'item_name') for which the item location is time invariant
 
SPLIT: a list (separated by blanks) of items (given by their time 1 'item_name') to split according to the corresponding time 1 items.

OUT: the name (default MML) given to output files

		data set 'out_thresholds' contains the item parameters (using PCM parametrization).
		data set 'out_poppar' contains population parameter estimates 
		data set 'out_logl' contains the maximum value of the log likelihood

**************************************************************************/

%macro lrasch_mml(DATA=data, ITEM_NAMES=item_names, ANCHOR=ALL, SPLIT=NONE, ICC=NO, OUT=LRASCH);

	options nomprint nonotes errors=0;
	ods listing close;
	goptions reset=all;
	title ' ';
	/**************************************************/
	/* small help macro to split for local dependence */  
	/**************************************************/
	%macro item_split_ld(data, item_names);
		/* data: data set with item responses (one variable for each item) */
		/* item_names: dataset with names of items to be splitted */
		proc sql noprint;
			select count(distinct(item_names1)), count(distinct(item_names2))
			into :__nitems1, :__nitems2
			from &item_names.;
		quit;
		%let __nitems1=&__nitems1.; %let __nitems2=&__nitems2.;
		
		proc sql noprint;
			select distinct(item_names1), max1
			into :__item1_1-:__item1_&__nitems1., :__max1_1-:__max1_&__nitems1.
			from &item_names.;
		quit;

		proc sql noprint;
			select distinct(item_names2), max2
			into :__item2_1-:__item2_&__nitems2., :__max2_1-:__max2_&__nitems2.
			from &item_names.;
		quit;

		%do _i=1 %to &__nitems1.; %let __max1_&_i=&&__max1_&_i; %end; 
		%do _i=1 %to &__nitems2.; %let __max2_&_i=&&__max2_&_i; %end;
		data &data._split;
			set &data.;
			%do _i=1 %to &__nitems1.; %do _h=0 %to &&__max1_&_i;
				if &&__item1_&_i=&_h. then &&__item2_&_i.._&_h.=&&__item2_&_i.;
				else &&__item2_&_i.._&_h.=.;
			%end; %end;	 
		run; 
	%mend item_split_ld;
	/***************************************************/
	/* end of help macro to split for local dependence */
	/***************************************************/
	%put;
	%put ------------------------;
	%put Longitudinal Rasch model;
	%put ------------------------;
	%put;

	/*************************************************/
	/* save names of anchor items as macro variables */
	/*************************************************/
	%if %upcase(&anchor.)=ALL %then %do;
		%put all items anchored;
		%put ;
	%end;
	%if %upcase(&anchor.)=NONE %then %do;
		%put no items anchored;
		%put ;
		%let _nanchor=0;
	%end;
	%if %upcase(&anchor.)^=ALL %then %do;
		%if %upcase(&anchor.)^=NONE %then %do;
			%let _anch1=%scan(&anchor.,1,' ');
			%let _nanchor=1;
			* if more than one item;
			%if %scan(&anchor.,1,' ')^=%scan(&anchor.,-1,' ') %then %do;
				%let _i=1;
				%do %until (%scan(&anchor.,%eval(&_i.),' ')=%scan(&anchor.,-1,' '));
					%let _anch%eval(&_i.+1)=%scan(&anchor.,%eval(&_i.+1),' ');
					%let _i=%eval(&_i.+1);
				%end;
				%let _nanchor=%eval(&_i.);
			%end;
			%put &_nanchor item(s) anchored;
			%do _a=1 %to &_nanchor; %put - &&_anch&_a; %end;
		%end;
	%end;
	/*********************************************/
	/* save split items as macro variables       */
	/*********************************************/
	%if %upcase(&split.)=NONE %then %do;
		%put no items are split;
		%put ;
	%end;
	%if %upcase(&split.)^=NONE %then %do;
		%let _spl1=%scan(&split.,1,' ');
		%let _nsplit=1;
		/* if more than one item */
		%if %scan(&split.,1,' ')^=%scan(&split.,-1,' ') %then %do;
			%let _i=1;
			%do %until (%scan(&split.,%eval(&_i.),' ')=%scan(&split.,-1,' '));
				%let _spl%eval(&_i.+1)=%scan(&split.,%eval(&_i.+1),' ');
				%let _i=%eval(&_i.+1);
			%end;
			%let _nsplit=%eval(&_i.);
		%end;
			%put &_nsplit item(s) are split;
			%do _s=1 %to &_nsplit; %put - &&_spl&_s; %end;
			%put ;
	%end;
	/* register the ordering of the input data */
	data _data; set &data.; order=_n_; run;
	/* Add column to 'ITEM_NAMES' dataset indicating whether item should be split */
	data &item_names._2;
		set &item_names.;
		format _split _anchor $1.;
		%if %upcase(&split.)=NONE %then %do; _split=' '; %end;
		%else %do; %do _s=1 %to &_nsplit.; if lowcase(item_names1)=%lowcase("&&_spl&_s") then _split='Y'; %end; %end;
		%if %upcase(&anchor.)=ALL %then %do; _anchor='Y'; %end;
		%else %do;
			%do _t=1 %to &_nanchor.; if lowcase(item_names1)=%lowcase("&&_anch&_t") then _anchor='Y'; %end;
		%end;
	run;
	/********************************************************************************************************/
	/* numbering of items, find maximum length of values - max(max) = highest number of response categories */
	/********************************************************************************************************/
	proc sql noprint;
		select count(distinct(item_names1))
		into :_nitems1
		from &item_names._2
		where item_names1^=' ';
	quit; 
	proc sql noprint;
		select count(distinct(item_names2))
		into :_nitems2
		from &item_names._2
		where item_names2^=' ';
	quit;
	%let _nitems1=&_nitems1.;
	%let _nitems2=&_nitems2.;
	proc sql noprint;
		select distinct(item_names1), item_names2, max1, max2, _split, _anchor
		into :_item1_1-:_item1_&_nitems1., 
			 :_item2_1-:_item2_&_nitems2., 
			 :_max1_1-:_max1_&_nitems1., 
			 :_max2_1-:_max2_&_nitems2., 
			 :_split1-:_split&_nitems2., 
			 :_anchor1-:_anchor&_nitems2.
		from &item_names._2;
	quit;
	%if %upcase(&split.)^=NONE %then %do;
		/**************************************************************************************************/
		/* names of time1- and time2-items of the 'split' items (needed for 'item_split_ld' macro)       */
		/**************************************************************************************************/
		proc sql noprint;
			select distinct(item_names1), item_names2
			into :_indep1-:_indep&_nsplit., :_dep1-:_dep&_nsplit.
			from &item_names._2
			where _split='Y';
		quit;
		%item_split_ld(data=_data, item_names=&item_names._2(where=(_split='Y')));
		/****************************/
		/* new 'item_names' dataset */
		/****************************/
		data item_names_spl1(drop=item_names2);
			set &item_names._2;
			format split_cat 2.;
			%do _i=1 %to &_nitems2.;
				if _split='Y' then do;
					%do _h=0 %to &&_max1_&_i;
						item_names2_temp=strip(item_names2)||"_"||"cat&_h.";
						split_cat="&_h.";
						output;
					%end;
				end;
				else do;
					item_names2_temp=item_names2;
					output;
				end;
			%end;
		run;

		proc sql;
			create table &item_names._split as select distinct(item_names2_temp) as 
			item_names2,
			item_names1,
			max1,
			max2,
			split_cat
			from item_names_spl1; 
		quit;
	%end;
	%else %do;
		data &item_names._split; set &item_names.; run;
		data _data_split; set _data; run;
	%end;

	/**************************************************************************************/
	/* data set with one line for each item response (more than one line for each person) */
	/**************************************************************************************/
	data _new; 
		format value 3. item $20.;
		set _data_split; 
		%do _i=1 %to &_nitems1.; 
			item="&&_item1_&_i"; 
			value=&&_item1_&_i; 
			person=_N_; 
			output;
		%end;
		%do _i=1 %to &_nitems2.; 
			%if %upcase(&&_split&_i)^=Y %then %do;
				item="&&_item2_&_i"; 
				value=&&_item2_&_i; 
				person=_N_; 
				output;
			%end;
			%if %upcase(&&_split&_i)=Y %then %do;
				%do _h=0 %to &&_max1_&_i;
					item="&&_item2_&_i.._&_h"; 
					value=&&_item2_&_i.._&_h; 
					person=_N_; 
					output;
				%end;
			%end;
		%end;
	run;
	/*******************************************/
	/* start of item parameter estimation part */
	/*******************************************/
	%put estimating item parameters using 2-dimensional MML;
	%put;
	/*******************/
	/* starting values */
	/*******************/
		proc nlmixed data=_new;
		ods output nlmixed.parameterestimates=_pe0;
		parms %do _i=1 %to &_nitems1.; %do _h=1 %to &&_max1_&_i; eta1_&_i._&_h=0, %end; %end;
			%do _i=1 %to &_nitems2.; 
				%if &&_split&_i^=Y %then %do;
					%if &&_anchor&_i^=Y %then %do; %do _h=1 %to &&_max2_&_i; eta2_&_i._&_h=0, %end; %end;
				%end;
				%if &&_split&_i=Y %then %do;
					%if &&_anchor&_i^=Y %then %do; 
						%do _h1=0 %to &&_max1_&_i; %do _h2=1 %to &&_max1_&_i; 
							eta2_&_i._&_h1._&_h2=0, 
						%end; %end;
					%end;
					%if &&_anchor&_i=Y %then %do; 
						%do _h1=0 %to &&_max1_&_i; %do _h2=1 %to &&_max1_&_i; 
							eta2_&_i._&_h1._&_h2=0, 
						%end; %end;
					%end;
				%end;
			%end;
			; 
			_theta1=0; _theta2=mu;
		/* time 1 likelihood */
		%do _i=1 %to &_nitems1.; 
			_denom=1 %do _k=1 %to &&_max1_&_i; +exp(&_k*_theta1+eta1_&_i._&_k) %end;;
			if item="&&_item1_&_i" and value=0 then ll=-log(_denom);
			%do _h=1 %to &&_max1_&_i;
				if item="&&_item1_&_i" and value=&_h then ll=(&_h*_theta1+eta1_&_i._&_h.)-log(_denom);
			%end; 
		%end;
		/* time 2 likelihood */
		%do _i=1 %to &_nitems2.;
			%if &&_split&_i=Y %then %do;
				%if &&_anchor&_i=Y %then %do;
					%do _h1=0 %to &&_max1_&_i;	
						_denom=1 %do _h2=1 %to &&_max2_&_i; +exp(&_h2*_theta2+eta2_&_i._&_h1._&_h2) %end;;
						if item="&&_item2_&_i.._&_h1" and value=0 then ll=-log(_denom);
						%do _h2=1 %to &&_max2_&_i;
							if item="&&_item2_&_i.._&_h1" and value=&_h2 then ll=(&_h2*_theta2+eta2_&_i._&_h1._&_h2.)-log(_denom);
						%end; 
					%end;
				%end;
				%if &&_anchor&_i^=Y %then %do;
					%do _h1=0 %to &&_max1_&_i;	
						_denom=1 %do _h2=1 %to &&_max2_&_i; +exp(&_h2*_theta2+eta2_&_i._&_h1._&_h2) %end;;
						if item="&&_item2_&_i.._&_h1" and value=0 then ll=-log(_denom);
						%do _h2=1 %to &&_max2_&_i;
							if item="&&_item2_&_i.._&_h1" and value=&_h2 then ll=(&_h2*_theta2+eta2_&_i._&_h1._&_h2.)-log(_denom);
						%end; 
					%end;
				%end;
			%end;
			%if &&_split&_i^=Y %then %do;
				%if &&_anchor&_i=Y %then %do;
					_denom=1 %do _h2=1 %to &&_max2_&_i; +exp(&_h2*_theta2+eta1_&_i._&_h2) %end;;
					if item="&&_item2_&_i.." and value=0 then ll=-log(_denom);
					%do _h2=1 %to &&_max2_&_i;
						if item="&&_item2_&_i.." and value=&_h2 then ll=(&_h2*_theta2+eta1_&_i._&_h2.)-log(_denom);
					%end; 
				%end;
				%if &&_anchor&_i^=Y %then %do;
					_denom=1 %do _h2=1 %to &&_max2_&_i; +exp(&_h2*_theta2+eta2_&_i._&_h2) %end;;
					if item="&&_item2_&_i.." and value=0 then ll=-log(_denom);
					%do _h2=1 %to &&_max2_&_i;
						if item="&&_item2_&_i.." and value=&_h2 then ll=(&_h2*_theta2+eta2_&_i._&_h2.)-log(_denom);
					%end; 
				%end;
			%end;
		%end;
		model value~general(ll);
	run;
	data __pe0;
		Parameter='rho'; Estimate=0.5;
	run;
	data _pe0; set _pe0 __pe0; run;
	/***********************/
	/* use starting values */
	/***********************/
	
	proc nlmixed data=_new;
		ods output nlmixed.convergencestatus=_cs;
		ods output nlmixed.parameterestimates=_pe;
		ods output nlmixed.additionalestimates=_ae;
		ods output nlmixed.fitstatistics=_logl;
		parms / data=_pe0;
		/* time 1 likelihood */
		%do _i=1 %to &_nitems1.; 
			_denom=1 %do _k=1 %to &&_max1_&_i; +exp(&_k*_theta1+eta1_&_i._&_k) %end;;
			if item="&&_item1_&_i" and value=0 then ll=-log(_denom);
			%do _h=1 %to &&_max1_&_i;
				if item="&&_item1_&_i" and value=&_h then ll=(&_h*_theta1+eta1_&_i._&_h.)-log(_denom);
			%end; 
		%end;
		/* time 2 likelihood */
		%do _i=1 %to &_nitems2.;
			%if &&_split&_i=Y %then %do;
				%if &&_anchor&_i=Y %then %do;
					%do _h1=0 %to &&_max1_&_i;	
						_denom=1 %do _h2=1 %to &&_max2_&_i; +exp(&_h2*_theta2+eta2_&_i._&_h1._&_h2) %end;;
						if item="&&_item2_&_i.._&_h1" and value=0 then ll=-log(_denom);
						%do _h2=1 %to &&_max2_&_i;
							if item="&&_item2_&_i.._&_h1" and value=&_h2 then ll=(&_h2*_theta2+eta2_&_i._&_h1._&_h2.)-log(_denom);
						%end; 
					%end;
				%end;
				%if &&_anchor&_i^=Y %then %do;
					%do _h1=0 %to &&_max1_&_i;	
						_denom=1 %do _h2=1 %to &&_max2_&_i; +exp(&_h2*_theta2+eta2_&_i._&_h1._&_h2) %end;;
						if item="&&_item2_&_i.._&_h1" and value=0 then ll=-log(_denom);
						%do _h2=1 %to &&_max2_&_i;
							if item="&&_item2_&_i.._&_h1" and value=&_h2 then ll=(&_h2*_theta2+eta2_&_i._&_h1._&_h2.)-log(_denom);
						%end; 
					%end;
				%end;
			%end;
			%if &&_split&_i^=Y %then %do;
				%if &&_anchor&_i=Y %then %do;
					_denom=1 %do _h2=1 %to &&_max2_&_i; +exp(&_h2*_theta2+eta1_&_i._&_h2) %end;;
					if item="&&_item2_&_i.." and value=0 then ll=-log(_denom);
					%do _h2=1 %to &&_max2_&_i;
						if item="&&_item2_&_i.." and value=&_h2 then ll=(&_h2*_theta2+eta1_&_i._&_h2.)-log(_denom);
					%end; 
				%end;
				%if &&_anchor&_i^=Y %then %do;
					_denom=1 %do _h2=1 %to &&_max2_&_i; +exp(&_h2*_theta2+eta2_&_i._&_h2) %end;;
					if item="&&_item2_&_i.." and value=0 then ll=-log(_denom);
					%do _h2=1 %to &&_max2_&_i;
						if item="&&_item2_&_i.." and value=&_h2 then ll=(&_h2*_theta2+eta2_&_i._&_h2.)-log(_denom);
					%end; 
				%end;
			%end;
		%end;
		model value~general(ll);
		random _theta1 _theta2 ~ normal([0,mu],[sigma1*sigma1,rho*sigma1*sigma2,sigma2*sigma2]) subject=person;	
		/*****************************/
		/* estimate statement time 1 */
		/*****************************/
		%do _i=1 %to &_nitems1.; 
			%do _h=1 %to &&_max1_&_i; 
				estimate "&&_item1_&_i.|&_h." -eta1_&_i._&_h. %if &_h.>1 %then %do; +eta1_&_i._%eval(&_h.-1) %end;; 
			%end; 
		%end;		
		/***********************************************/
		/* estimate statement time 2 (for split items) */
		/***********************************************/
		%do _i=1 %to &_nitems1.; 
			%if &&_split&_i=Y %then %do;
				%if &&_anchor&_i=Y %then %do;
					%do _h1=0 %to &&_max1_&_i;	
						%do _h2=1 %to &&_max1_&_i;
							estimate "&&_item2_&_i.|&_h2. (&&_item1_&_i = &_h1.)" -eta2_&_i._&_h1._&_h2.
							%if &_h2.>1 %then %do; +eta2_&_i._&_h1._%eval(&_h2.-1) %end; 
							; 
						%end;
					%end; 
				%end;
				%if &&_anchor&_i^=Y %then %do;	
					%do _h1=0 %to &&_max1_&_i;
						%do _h2=1 %to &&_max2_&_i;	
							estimate "&&_item2_&_i.|&_h2. (&&_item1_&_i = &_h1.)" -eta2_&_i._&_h1._&_h2. 
							%if &_h2.>1 %then %do; +eta2_&_i._&_h1._%eval(&_h2.-1) %end; 
							; 
						%end;
					%end;
				%end;
			%end;
		%end;
		/***********************************************/
		/* estimate statement time 2 (for other items) */
		/***********************************************/
		%do _i=1 %to &_nitems1.; 
			%if &&_split&_i^=Y %then %do;
				%if &&_anchor&_i^=Y %then %do;
					%do _h2=1 %to &&_max2_&_i; 
						estimate "&&_item2_&_i|&_h2." -eta2_&_i._&_h2. 
						%if &_h2.>1 %then %do; +eta2_&_i._%eval(&_h2.-1) %end;
						; 
					%end; 
				%end;
				%if &&_anchor&_i=Y %then %do;
					%do _h2=1 %to &&_max2_&_i; 
						estimate "&&_item2_&_i|&_h2." -eta1_&_i._&_h2. 
						%if &_h2.>1 %then %do; +eta1_&_i._%eval(&_h2.-1) %end;
						; 
					%end; 
				%end;
			%end;
		%end;
	run;
	/*****************************************/
	/* end of item parameter estimation part */
	/*****************************************/
	data _pe_orig; set _pe; run;
	/* log likelihood */
	data &out._logl; 
		set _logl; 
	run;
	/* population parameters */
	data &out._poppar;
		set _pe;
		if substr(parameter,1,3)='eta' then delete;
		drop df tvalue probt alpha gradient;
	run;
	/* data with item parameters at the time 1 */
	data _eta1;
		set _pe;
		if substr(parameter,1,4)='eta1';
		keep parameter estimate;
	run;
	data _ipar1;
		set _eta1;
		item_name='                                                 ';
		if char(parameter,7)='_' then item_no=substr(parameter,6,1)*1; else item_no=substr(parameter,6,2)*1; 
		temp=scan(parameter,1,'|');
		item_name=substr(temp,6,14);
		if char(parameter,7)='_' then score=substr(item_name,3,1)*1; else score=substr(item_name,4,1)*1;
	run;
	/* data with item parameters at the time 2 */
	data _eta2;
		set _pe;
		if substr(parameter,1,4)='eta2';
		keep parameter estimate;
	run;
	/* '_ipar2' (SPLIT=N, ANCHOR=Y) */
	data _ipar2_NY;
		set _ipar1;
		%do _i=1 %to &_nitems1;
			%if &&_split&_i=Y %then %do; if item_no=&_i then delete; %end;
			%if &&_anchor&_i^=Y %then %do; if item_no=&_i then delete; %end;
		%end;
		item_name_='                                                 ';
	run;
	data _ipar2_NY;
		set _ipar2_NY;
		%do _i=1 %to &_nitems2.;
			if item_no=&_i then do; max=&&_max2_&_i; item_name="&&_item2_&_i"; end;
		%end;
		if substr(parameter,1,4)='eta1' then substr(parameter,1,4)='eta2';
	run;
	/* '_ipar2' (SPLIT=N, ANCHOR=N)*/
	data _ipar2_NN;
		set _pe;
		if substr(parameter,1,4)='eta2';
		item_name='                                                 ';
		if char(parameter,7)='_' then item_no=substr(parameter,6,1)*1; 
		else item_no=substr(parameter,6,2)*1; 
		keep parameter estimate item_no item_name;
	run;
	data _ipar2_NN;
		set _ipar2_NN;
		%do _i=1 %to &_nitems1;
			%if &&_split&_i=Y %then %do; if item_no=&_i then delete; %end;
			%if &&_anchor&_i=Y %then %do; if item_no=&_i then delete; %end;
		%end;
	run;
	data _ipar2_NN;
		set _ipar2_NN;
		item_name_='                                                 ';
		if char(parameter,7)='_' then score=substr(parameter,8,1)*1; 
		else score=substr(parameter,9,1)*1; 
	run;
	data _ipar2_NN;
		set _ipar2_NN;
		%do _i=1 %to &_nitems2.;
			if item_no=&_i then do; max=&&_max2_&_i; item_name="&&_item2_&_i"; end;
		%end;
	run;
	/* '_ipar2' (SPLIT=Y)*/
	data _ipar2_Y;
		set _pe;
		if substr(parameter,1,4)='eta2';
		item_name='                                                 ';
		if char(parameter,7)='_' then item_no=substr(parameter,6,1)*1; 
		else item_no=substr(parameter,6,2)*1;
		keep parameter estimate item_no;
	run;
	data _ipar2_Y;
		set _ipar2_Y;
		item_names='                                                 ';
		%do _i=1 %to &_nitems2.;
			if item_no=&_i then do; item_names="&&_item2_&_i"; max=&&_max2_&_i; end;
		%end;
	run;
	data _ipar2_Y;
		set _ipar2_Y;
		%do _i=1 %to &_nitems1;
			%if &&_split&_i^=Y %then %do; if item_no=&_i then delete; %end;
		%end;
	run;
	data _ipar2_Y;
		set _ipar2_Y;
		if char(parameter,7)='_' then item_no_=substr(parameter,6,3); 
		else item_no_=substr(parameter,6,4); 
		if char(parameter,9)='_' then score=substr(parameter,10,1)*1; 
		else score=substr(parameter,11,1)*1; 
	run;
	/* attach artifical item numbers to the split items*/
	data _ipar2_Y;
		set _ipar2_Y;
		%do _i=1 %to 9;
			if substr(item_no_,1,1)="&_i" and substr(item_no_,2,1)='_' then do; 
				item_no=100*&_i+substr(item_no_,3,1);
			end;
		%end;
		%do _i=10 %to &_nitems1;
			if substr(item_no_,1,2)="&_i" and substr(item_no_,3,1)='_' then do; 
				item_no=100*&_i+substr(item_no_,4,1);
			end;
		%end;
		split='Y';
		if substr(item_no_,2,1)='_' then do; old_item_no=substr(item_no_,1,1)*1; splitval=substr(item_no_,3,1)*1; end;
		if substr(item_no_,3,1)='_' then do; old_item_no=substr(item_no_,1,2)*1; splitval=substr(item_no_,4,1)*1; end;
	run;
	/* data sets for later use */
	data _ipar1;
		set _ipar1;
		parameter='eta'||strip(item_no)||'_'||strip(score);
		%do _i=1 %to &_nitems1;
			if item_no=&_i then do; max=&&_max1_&_i; item_name="&&_item1_&_i"; end;
		%end;
		drop temp;
	run; 
	/* TIME2: recode item numbers for unsplit and split items to get consecutive integers */
	data _ipar2;
		set _ipar2_nn _ipar2_ny _ipar2_y;
		old_item_no=item_no;
		drop item_name_ temp;
	run;
	proc sql noprint;
		select count(distinct(item_no))
		into :___nitems_
		from _ipar2;
	quit;
	%let ___nitems_=&___nitems_.;
	proc sql noprint;
		select unique(item_no) 
		into :___item_no_1-:___item_no_&___nitems_ 
		from _ipar2;
	quit;
	%do _i=1 %to &___nitems_.;
		%let ___item_no_&_i=&&___item_no_&_i.;
	%end;
	data _ipar2;
		set _ipar2;
		%do _i=1 %to &___nitems_.;
			if item_no=&&___item_no_&_i then item_no=&_i;
		%end;
	run;
	data _ipar2;
		set _ipar2;
		parameter='eta'||strip(item_no)||'_'||strip(score);
	run;
	/* 'item_name' for split items */
	data _ipar2;
		set _ipar2;
		if item_name='' and substr(item_no_,2,1)='_' then item_name=strip(item_names)||strip(substr(item_no_,2,2));
		if item_name='' and substr(item_no_,3,1)='_' then item_name=strip(item_names)||strip(substr(item_no_,3,2));
	run;

	/* output dataset with thresholds */
	data &out._thresholds; 
		set _ae(keep=label estimate standarderror lower upper); 
	run;
	/* range of estimated thresholds */
	proc sql noprint;
		select round(min(estimate)-1), round(max(estimate)+1)
		into :xmin, :xmax
		from &out._thresholds;
	quit;
	%let xmin=&xmin.;
	%let xmax=&xmax.;

	/************************************/
	/* start of part for plotting ICC's */
	/************************************/

	%if %upcase(&icc.)^=NO %then %do; 
		data _zeros;
			%do _i=1 %to &_nitems1; item_no=&_i; score=0; estimate=0; output; %end;
		run;
		data _icc1; 
			set _ipar1 _zeros; 
		run;
		data _icc1; 
			set _icc1; 
			do _theta=%eval(&xmin.) to %eval(&xmax.) by 0.01; output; end;
		run;
		proc sql;
			create table _plot as select 
			parameter,
			item_no,
			score,
			estimate,
			_theta,
			exp(score*_theta+estimate)/(sum(exp(score*_theta+estimate))) as _prob
			from _icc1
			group by item_no, _theta
			order by item_no, score, _theta;
		quit;
		/**************************************/
		/* plot time 1 ICC's                  */
		/**************************************/
		ods listing;
			%do _i=1 %to &_nitems1;
				%if &&_anchor&_i=Y %then %do; %if &&_split&_i^=Y %then %do; 
					title "&&_item1_&_i and &&_item2_&_i";
				%end; %end;
				%if &&_anchor&_i^=Y %then %do; 
					title "&&_item1_&_i";
				%end;
				%if &&_split&_i=Y %then %do; 
					title "&&_item1_&_i";
				%end;
				axis1 order=(%eval(&xmin.) to %eval(&xmax.) by 1) length=15 cm value=(H=2) minor=NONE label=(H=2 'Latent variable');
				axis2 label=(H=2 A=90 'Probability') length=10 cm value=(H=2) minor=NONE;
				symbol1 v=none i=join l=1 w=1 c=black w=3 r=%eval(&&_max1_&_i+1);
				proc gplot data=_plot (where=(item_no=&_i.));
					plot _prob*_theta=score / haxis=axis1 vaxis=axis2 nolegend;
				run; quit;
			%end;
		ods listing close;
		/**************************************/
		/* plot time 2 ICC's for split items  */
		/**************************************/
		data _zeros2; 
			%do _i=1 %to &_nitems1; 
				%do _sv=0 %to &&_max1_&_i;
					old_item_no=&_i*1; splitval=&_sv*1; score=0; estimate=0; split="&&_split&_i"; output; 
				%end;
			%end;
		run;
		data _icc2_spl; 
			set _ipar2_Y _zeros2(where=(split='Y')); 
		run;
		data _icc2_spl; 
			set _icc2_spl; 
			do _theta=%eval(&xmin.) to %eval(&xmax.) by 0.01; output; end;
		run;
		proc sql;
			create table _plot_spl as select 
			parameter,
			old_item_no,
			splitval,
			score,
			estimate,
			_theta,
			exp(score*_theta+estimate)/(sum(exp(score*_theta+estimate))) as _prob
			from _icc2_spl
			group by old_item_no, splitval, _theta
			order by old_item_no, splitval, score, _theta;
		quit;
		ods listing;
			%do _i=1 %to &_nitems1;
				%if &&_split&_i=Y %then %do; 
					%do _h1=0 %to &&_max1_&_i;
						title "&&_item2_&_i. (&&_item1_&_i = &_h1.)";
						axis1 order=(%eval(&xmin.) to %eval(&xmax.) by 1) value=(H=2) minor=NONE label=(H=2 'Latent variable');
						axis2 label=(H=2 A=90 'Probability') value=(H=2) minor=NONE;
						symbol1 v=none i=join l=1 w=1 c=black w=3 r=%eval(&&_max1_&_i+1);
						proc gplot data=_plot_spl;
							plot _prob*_theta=score / haxis=axis1 vaxis=axis2 nolegend;
							where (old_item_no=&_i.) and (splitval=&_h1.);
						run; quit;
					%end;
				%end;
			%end;
		ods listing close;
	%end; 
	/*******************************************/
	/* plot time 2 ICC's for unanchored items  */
	/*******************************************/
	data _zeros2_unan; 
		%do _i=1 %to &_nitems1; 
			item_no=&_i*1; score=0; estimate=0; split="&&_split&_i"; anch="&&_anchor&_i"; output; 
		%end;
	run;
	data _icc2_unan; 
		set _ipar2_nn _zeros2_unan(where=(split^='Y' and anch^='Y')); 
	run;
	data _icc2_unan; 
		set _icc2_unan; 
		do _theta=%eval(&xmin.) to %eval(&xmax.) by 0.01; output; end;
	run;
	proc sql;
		create table _plot as select 
		parameter,
		item_no,
		score,
		estimate,
		_theta,
		exp(score*_theta+estimate)/(sum(exp(score*_theta+estimate))) as _prob
		from _icc2_unan
		group by item_no, _theta
		order by item_no, score, _theta;
	quit;
	ods listing;
		%do _i=1 %to &_nitems1;
			%if &&_anchor&_i^=Y %then %do; 
				%if &&_split&_i^=Y %then %do; 
					title "&&_item2_&_i";
					axis1 order=(%eval(&xmin.) to %eval(&xmax.) by 1) length=15 cm value=(H=2) minor=NONE label=(H=2 'Latent variable');
					axis2 label=(H=2 A=90 'Probability') length=10 cm value=(H=2) minor=NONE;
					symbol1 v=none i=join l=1 w=1 c=black w=3 r=%eval(&&_max2_&_i+1);
					proc gplot data=_plot (where=(item_no=&_i.));
						plot _prob*_theta=score / haxis=axis1 vaxis=axis2 nolegend;
					run; quit;
				%end; 
			%end;
		%end;
	ods listing close;

	/**********************************/
	/* end of part for plotting ICC's */
	/**********************************/

	/******************************/
	/* delete temporary data sets */
	/******************************/
	ods output Datasets.Members=_mem;
	proc datasets; run; quit;

	data _mem; set _mem; _underscore=substr(name,1,1); run;
	data _mem; set _mem(where=(_underscore='_')); run;
	data _null_; set _mem end=final; if final then call symput('_nd',trim(left(_N_))); run;

	%put clean-up deleting &_nd temporary data sets;
	proc sql noprint;
		select distinct(name)
		into :_var1-:_var&_nd.
		from _mem;
	quit;
	proc datasets;
		delete %do _d=1 %to &_nd; &&_var&_d %end;;
	run; quit;
	/**************************/
	/* print output data sets */
	/**************************/
	ods listing;
	options nocenter;
	title "&out: LRASCH MML estimation loglikelihood";
	proc print data=&out._logl round noobs; run;
	title "&out: LRASCH MML estimated thresholds";
	proc print data=&out._thresholds (keep=Label Estimate StandardError lower upper ) noobs; 
		format Estimate StandardError lower upper 10.2; 
	run;
	title "&out: LRASCH change in latent mean (mu), variances, and latent correlation (rho)";
	proc print data=&out._poppar noobs; 
		format Estimate StandardError lower upper 10.3; 
	run;
	options notes;
	/****************/
	/* end of macro */
	/****************/
%mend lrasch_mml;
