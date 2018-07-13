# lrasch

A SAS macro that can be used to fit longitudinal Rasch models. Inlude the macro using
```
%let homep=https://publicifsv.sund.ku.dk/~kach/macro/lrasch_mml.sas;
filename lr URL "&homep";
%include lr;
```
data
```
proc import datafile="p:\public_html\macro\longi"
			dbms=sav 
			out=longi 
			replace;
run;
data longi;
	set longi;
	rename item1response1=i1_1; rename item1response2=i1_2;
	rename item2response1=i2_1; rename item2response2=i2_2;
	rename item3response1=i3_1; rename item3response2=i3_2;
	rename item4response1=i4_1; rename item4response2=i4_2;
	rename item5response1=i5_1; rename item5response2=i5_2;
	rename item6response1=i6_1; rename item6response2=i6_2;
	rename item7response1=i7_1; rename item7response2=i7_2;
	rename item8response1=i8_1; rename item8response2=i8_2;
	rename item9response1=i9_1; rename item9response2=i9_2;
	rename item10response1=i10_1; rename item10response2=i10_2;
run;
```
specify variable names and response formats (item SC01 at baseline and items SC03 and SC06 at follow-up considered to have only three response options)
```
data inames;
input item_no item_names1 $ item_names2 $ item_text $ max1 max2;
datalines;
1 i1_1 i1_2 x 1 1
2 i2_1 i2_2 x 1 1
3 i3_1 i3_2 x 1 1
4 i4_1 i4_2 x 1 1
5 i5_1 i5_2 x 1 1
6 i6_1 i6_2  x 1 1
7 i7_1 i7_2 x 1 1
8 i8_1 i8_2 x 1 1
9 i9_1 i9_2 x 1 1
10 i10_1 i10_2 x 1 1
;
run;
```
fit the longitudinal Rasch model - plot ICC's
```
%lrasch_mml(DATA=longi, 
            ITEM_NAMES=inames, 
            ANCHOR=ALL, 
            SPLIT=NONE, 
            OUT=LRM,
			ICC=YES);
```
fit longitudinal Rasch model - allow item parameter drift for item 1 
%lrasch_mml(DATA=longi, 
            ITEM_NAMES=inames, 
            ANCHOR=ALL, 
            SPLIT=i1_1, 
            OUT=DIF1);
```
compare the two models using a likelihood ratio test
```
data lrt; 
	merge LRM_logl(rename=(value=logl_0)) DIF1_logl(rename=(value=logl));
	where Descr='-2 Log Likelihood';
	lrt=logl_0-logl;
	df=3;
	p=1-probchi(lrt,df);
run;
proc print noobs; 
run;
```
fit longitudinal Rasch model - allow for local dependence across time points 
%lrasch_mml(DATA=data.scqol, 
            ITEM_NAMES=inames, 
            ANCHOR=ALL, 
            SPLIT=i1_1, 
            OUT=LD1);
```
compare the two models using a likelihood ratio test
```
data lrt;
    merge LRM_logl(rename=(value=logl_0)) LD1_logl(rename=(value=logl));
    where Descr='-2 Log Likelihood';
    lrt=logl_0-logl;
    df=12;
    p=1-probchi(lrt,df);
run;
proc print noobs; 
run;
```

