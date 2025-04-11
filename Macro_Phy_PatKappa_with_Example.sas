/** ~*~~*~~*~~*~~*~~*~~*~~*~~*~~*~~*~~*~~*~~*~~*~~*~~*~~*~~*~~*~~*~~*~~*~~*~~*~~*~~*~~*~** 
### SAS macro to calculate the Kappa statistic and its variance under clustered data
 		Code to create the results in Examples section of the manuscript;
### "Kappa statistic for clustered physician-patients polytomous data";
### for Computational Statistics and Data Analysis;
### by Zhao Yang, Ming Zhou;
### Version Date: 04Feb2014
** ~*~~*~~*~~*~~*~~*~~*~~*~~*~~*~~*~~*~~*~~*~~*~~*~~*~~*~~*~~*~~*~~*~~*~~*~~*~~*~~*~~*~**/;
/** ~*~~*~~*~~*~~*~~*~~*~~*~~*~~*~~*~~*~~*~~*~~*~~*~~*~~*~~*~~*~~*~~*~~*~~*~~*~~*~~*~~*~** 
##### Revision History

** ~*~~*~~*~~*~~*~~*~~*~~*~~*~~*~~*~~*~~*~~*~~*~~*~~*~~*~~*~~*~~*~~*~~*~~*~~*~~*~~*~~*~**/;
%macro Phy_PatKappa (dsin = , /* the input dataset */
					 clusV = , /* the cluster variable*/
					 unitV = , /* the unit variable*/
					 Resp1V =, /* Response variable for physicians, row variable, 
					 			need to be numeric integer and > 0 */
					 Resp2V =  /* Response variable for patients, column variable,
								need to be numeric integer and > 0 */);
	proc sql noprint;
		create table clustDS as
			select distinct &clusV from &dsin
			order by &clusV; quit;
	data clustDS;
		set clustDS; clusterID + 1;
	proc sort data = &dsin;
		by &clusV &unitV;
	data resp_;
		merge &dsin clustDS;
		by &clusV;
	proc sort data = resp_;
		by clusterID &unitV;
	data F_Ana_DS (drop = &clusV &unitV rename = (clusterID = cluster));
		set resp_;
		by clusterID &unitV;
		if first.clusterID then subject = 0;
			subject + 1; 
	data F_Ana_DS;
		set F_Ana_DS (rename = (&Resp1V = resp1__ &Resp2V = resp2__));
		newres = max(resp1__, resp2__);
	run;

	proc sql noprint;
		select max(newres) into: maxr
		from F_Ana_DS; quit;
	%let maxr2 = %eval(&maxr*&maxr);

	** put the code in another way;
	data F_Ana_DS (drop = i);
		set F_Ana_DS;
		%do i = 1 %to &maxr2;
			bs&i = .;
		%end;
		array basisV{*} bs1 -- bs&maxr2; ** resp1__ physician, resp2__ patients;
		do i = 1 to dim(basisV);
			if resp1__ = i - (resp2__ - 1)*&maxr then basisV{resp1__ + (resp2__ - 1)*&maxr} = 1;
			if basisV{i} = . then basisV{i} = 0;
		end;
	run;

	*** calculation ussing the Kappa under the independence (non-clusterd situation);
	data fakestr;
		do resp1__ = 1 to &maxr;
			do resp2__ = 1 to &maxr;
				count = 0;	output;
			end;
		end;
	proc sort data = fakestr;
		by resp1__ resp2__;
	run;

	proc summary data = F_Ana_DS;
		class resp1__ resp2__;
		output out = interm_D (where = (_TYPE_ = 3) rename = (_FREQ_ = count) );
	proc sort data = interm_D (drop = _TYPE_) out = interm_D;
		by resp1__ resp2__;
	data F_Ana_DSI;
		merge fakestr interm_D;
		by resp1__ resp2__;
	run;

	title2 "*** Part 1: Kappa statistics assuming independence *** ";	
	proc freq data = F_Ana_DSI;
		label resp1__ = "Physician"  resp2__ = "Patient";
		tables resp1__ * resp2__ /agree;
		weight count/zeros;
	run;
    /*****************************************************************************/;

	proc sql noprint;
		select max(cluster) into: nclus
		from F_Ana_DS; quit;

	data intermDS;
		set F_Ana_DS (keep = cluster bs1 -- bs&maxr2);
	run;
	proc sql noprint;
		create table clussize as
			select cluster, count(*) as clussize, 
					calculated clussize**2 as clussize2,
					calculated clussize*(calculated clussize - 1) as clussize3
			from intermDS
			group by cluster;
	quit;

	%do i = 1 %to &nclus;
		data clusds&i;
			set intermDS (where = (cluster = &i));
		proc sql noprint;
			create table clusds_n&i as /* this contains the total count for each cluster  */
				select distinct cluster, count(*) as clusize,
						%do j = 1 %to &maxr2 - 1; 
							sum(bs&j) as bs_s&j, 
						%end; 
						sum(bs&maxr2) as bs_s&maxr2
				from clusds&i;
			create table clusds_p&i as /* this contains the percentage for each cluster */
				select distinct cluster, count(*) as clusize,
						%do j = 1 %to &maxr2 - 1; 
							sum(bs&j)/calculated clusize as bs_s_p&j, 
						%end; 
						sum(bs&maxr2)/calculated clusize as bs_s_p&maxr2
				from clusds&i;
		quit;
	%end;

	data overall;  ** this is to calculate the overall percentage vector;
		set %do i = 1 %to &nclus; 
				clusds_n&i  
			%end;;
	proc sql noprint;
		create table overall_Sn as
			select sum(clusize) as total,
				%do j = 1 %to &maxr2 - 1; 
					sum(bs_s&j)/calculated total as bs_O_p&j, 
				%end; 
				sum(bs_s&maxr2)/calculated total as bs_O_p&maxr2
		from overall;
	quit;


/**************************************************************************************************************/;
/**************************************************************************************************************/;
/****************************** this starts the calculation of the statistic **********************************/;
	proc iml;
		use clussize;
		read all ;
		clusziV = clussize;
		clusziV2 = clussize2;
		clusziV3 = clussize3;
		close clussize;

/******************************************************************************************************/;
	*** this part is for calculating the rho_w ***;
	%do i = 1 %to &nclus;
		use clusds&i; 
		read all ;
		Clus_Level_M&i = bs1 %do j = 2 %to &maxr2; || bs&j %end; ; ** this is the unit level for each cluster;
		close clusds&i;

		use clusds_p&i;
		read all ;
		Clus_Overall_M&i = bs_s_p1 %do j = 2 %to &maxr2; || bs_s_p&j %end; ; ** this is the overall level for each cluster;
		clu_size = nrow(Clus_Level_M&i);
		close clusds_p&i;
		call symput( "clu_size&i", left( char(clu_size) ) );   
		
		clust_&i = SHAPE (Clus_Overall_M&i, &maxr, &maxr);
		clust_S&i = clust_&i[, +];

		use overall_Sn;
		read all ;
		Overall_M = bs_O_p1 %do j = 2 %to &maxr2; || bs_O_p&j %end; ; ** this is the overall ;
		bigN = total;
		close overall_Sn;
		call symput( "bigN", left( char(bigN) ) );   

		Overallp = SHAPE (Overall_M, &maxr, &maxr);
		Overallp_S = Overallp[, +];

/********* this is the basis to calculate the MSW ******************************************/;

		%do s = 1 %to &&clu_size&i;
			unitL = SHAPE (Clus_Level_M&i[&s,], &maxr, &maxr);
			unitL_S = unitL[, +]; 
			clu_mat_&i._&s = (unitL_S - clust_S&i)` * (unitL_S - clust_S&i);
		%end;

		clu_mat_&i = clu_mat_&i._1 
			%do s = 2 %to &&clu_size&i;
				+ clu_mat_&i._&s
			%end; ;

/******************************************************************************************************/;

/********* this is the basis to calculate the MSB ******************************************/;

		Overall_mat_&i = &&clu_size&i * (clust_S&i - Overallp_S)` * (clust_S&i - Overallp_S);

/******************************************************************************************************/;
	%end;

/********* this is the basis to calculate the rho_w ******************************************/;

	SSW = clu_mat_1 %do i = 2 %to &nclus; + clu_mat_&i	%end; ;
	SSB = Overall_mat_1 %do i = 2 %to &nclus; + Overall_mat_&i	%end; ;

	MSB = SSB/(&nclus - 1);
	MSW = SSW/(&bigN - &nclus);
	n0 = (&bigN - clusziV2[+]/&bigN)/(&nclus - 1);

	rho_w = TRACE(MSB - MSW)/TRACE( MSB + (n0 - 1) * MSW );  ** this is the value of rho_w;

/******************************************************************************************************/;

/********************* Define some vectors and matrices ******************************************/;

	iden = I(&maxr);
	L0 = iden[,1]` %do i = 2 %to &maxr; || iden[,&i]` %end;	;
    vec1 = J(&maxr, 1) ;
    L  = iden@vec1;
    U  = vec1@iden;
    B  = l0` || L || U ;
    vec2 = {1}||J(&maxr, 1, 0)`||J(&maxr, 1, 0)` ;

/******************************************************************************************************/;

	V0 = diag(Overall_M) - Overall_M` * Overall_M; ** this is for V0;

/********************* this part is for V1 ******************************************/;
	perc_v = SHAPE (Overall_M, &maxr, &maxr);
	mu_y_v = perc_v[, +];
	%do i = 1 %to &maxr;
		rowM_&i = perc_v[&i,]` * perc_v[&i,]/mu_y_v[&i];
	%end;

	part2 = BLOCK( rowM_1 %do i = 2 %to &maxr; ,rowM_&i %end; );

	V1 = rho_w * (part2 - Overall_M` * Overall_M);
/******************************************************************************************************/;

    po = l0 * Overall_M`;                               
    p_plus = L` * Overall_M`;
    pplus_ = U` * Overall_M`;
    pe = p_plus` * pplus_;
    kappa = (po - pe)/(1 - pe);

    V_Z_bar = V0/&bigN + (clusziV3[+]/&bigN**2)*V1 ;
    V_zeta = B` * V_Z_bar * B ;

	vec3 = {0}||pplus_`||p_plus`; 
	vec4 = vec2//vec3;

	V_pope = vec4 * V_zeta * vec4`;
    kappa_deriv = 1/(1 - pe) || (po - 1)/(1 - pe)**2;
    kappa_ase = sqrt( kappa_deriv * V_pope * kappa_deriv`) ;

	lower = kappa - 1.96 # kappa_ase;
	upper = kappa + 1.96 # kappa_ase;

	******* this part is for pooling the final results information;
		
	title2 "*** Part 2: Kappa statistics for clustered Physician-patients data *** ";	
	print PO PE rho_w [label="ICC"],,
			kappa kappa_ase[label="Standard Error"],, 
			 lower[label="95% CI Lower limit"] upper[label="95% CI Upper limit"];
	quit;
%mend;

/* here is the example for physician-patients dichotomous data **/;

data responseA;
 	input Cluster a b c d;
cards;
1  1  0  4  2
2  1  1  2  2
3  3  1  3  0
4  2  1  2  1
5  2  0  1  1
6  0  3  3  0
7  0  3  0  2
8  0  3  2  3
9  0  1  2  1
10  0  0  4  3
11  1  0  2  0
12  1  1  4  0
13  1  0  3  3
14  2  0  2  2
15  2  0  1  2
16  1  0  3  0
17  2  0  1  3
18  1  0  2  2
19  0  3  0  1
20  2  1  0  0
21  1  1  0  3
22  0  2  2  0
23  1  0  1  2
24  0  0  2  0
25  0  0  2  1
26  0  0  0  2
27  0  0  0  2
28  1  0  1  1
29  0  0  1  0
;
quit;
data responseA1 (keep = Cluster Subject Trt2 Resp2 Trt1 Resp1);
	set responseA;
	if a ne 0 then do;
		do Subject = 1 to a;	
			Trt1 = 1; Resp1 = 1;
			Trt2 = 2; Resp2 = 1;
			output;
		end;
	end;
run;

data responseA2 (keep = Cluster Subject Trt2 Resp2 Trt1 Resp1);
	set responseA;
	if b ne 0 then do;
		do Subject = 1 to b;	
			Trt1 = 1; Resp1 = 1;
			Trt2 = 2; Resp2 = 2;
			output;
		end;
	end;
run;

data responseA3 (keep = Cluster Subject Trt2 Resp2 Trt1 Resp1);
	set responseA;
	if c ne 0 then do;
		do Subject = 1 to c;	
			Trt1 = 1; Resp1 = 2;
			Trt2 = 2; Resp2 = 1;
			output;
		end;
	end;
run;

data responseA4 (keep = Cluster Subject Trt2 Resp2 Trt1 Resp1);
	set responseA;
	if d ne 0 then do;
		do Subject = 1 to d;	
			Trt1 = 1; Resp1 = 2;
			Trt2 = 2; Resp2 = 2;
			output;
		end;
	end;
run;

data response;
	set responseA1 responseA2 responseA3 responseA4;
run;

proc sort data = response;
	by Cluster;
run;

data new (keep = Cluster subjid Resp1 Resp2);
	set response;
	by Cluster;
	if first.Cluster then subjid = 0;
		subjid + 1;
run;

/*
proc print data = new noobs;
var Cluster subjid Resp1 Resp2;
run; */

%Phy_PatKappa (dsin = new, /* the input dataset */
				clusV = Cluster, /* the cluster variable*/
				unitV = subjid, /* the unit variable*/
				Resp1V = Resp1, /* Response variable for physicians, row variable
					 			need to be numeric integer and > 0 */
				Resp2V = Resp2 /* Response variable for patient, column variable 
								need to be numeric integer and > 0 */);

