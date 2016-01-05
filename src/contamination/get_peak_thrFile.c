 /*  File: get_signif_maf_thrFile.c // computes maf distribution over variant positions from vcf files
 * Authors: designed by Irina Abnizova (ia1)
 *

  Last edited: 3 feb-reads thr from a file; 29 Jan 2014  - get autom thr, starting with default ones 1,5
  still old fields of vcf input
  , 21 Nov 2013, for two thresholds

  outputs vector from function????? histAF  NO

  just mu, std and then Z value
 *-------------------------------------------------------------------
 * Description: computes maf distrib over variant positions from just extracted 6 fileld

 * Exported functions:
 * HISTORY:

computes AF (alternative to Ref allele frequency) for each error base call from vcf

applies minimum depth threshold, both for Ref and Alternative alleles fwd rev (total min depth will be
4*MIN_DEPTH

input1 thrfile  (thRR thAA)
input2: extracted vcf,


output1: info AF.txt=

chr  pos  %AF Q30frac (for this position)
6    1546 100 0.89
8    4469 100 0.86
22   4532 100 0.78
22   4559  97 0.95

output2: distrAF.txt=histogram of AF with significance for each bin
% count signif
8 0 0
9 0 0
10 385 1
11 780 3


usage:./get_signif_maf name.thr extracted.vcf .infAF .distrAF

*/



#include <stdio.h>
//#include "conio.h"
#include <time.h>
#include <math.h>
#include <string.h>

// ******** predefined user constants, 22 october

#define NRID 6      // Number of variant ids=fields: pos, chromosome, Depth,DP4,(no PV4+maf) in extracted_vcf files
#define NPAR 5      // max Number of PARameters and arguments
#define Nthr 2      //  Number of thresholds


// ******** declarations of  functions
int   GetAf (int []);// computes AF percent for a variant position
float   GetMu( int data[]);
float   GetStd( int data[]);
//int GetMinPercAF(int data[]);
//int GetSignifPercAF(int data[], int data1[]);
int GetPeakPercAF(int data[]);

int main    (int argc, char *argcv[]) {

    // flags
    int firstSixAreOK  = 1;
   // int lastSevenAreRight = 1;
    int canWriteAF        = 1;
    int canWriteDistrib      = 1;

    FILE *extract_vcfFile, *afFile, *distribFile, *thrFile; //file handles

    int n,i,perc;//, cntZero=0, counter=0;
    //char chrom[50];
    int count_afterF=0;// after filtering with min depths
    int count_beforeF=0;// before filtering with min depths
    int DP4[4];//
    int D,pos,AF,Z;// Depth, genome pos of error
    int sumD=0,sumDP4=0;//across all pos!
    int sumDbefore=0;
    //int thRR,thAA;// to compute!=optimal: given here
    int thR,thA;// from file

    float mu,std;//mean and std of alternative frequency, AF
    float avDbefore;
    float avD,avDP4,Q30frac;
    float sDP4;
    float perc_left;
    float x,conf;//confidence (0,1)

    int histAF[100];// to store mafs percentages
    int ZZ[100];//signif of each histo bin, integerised
    int THR[2];// define from input1 thr file


    if(argc < NPAR)//four input_output files are submitted
    {
        printf("not enough of parms |input_output files\n");
        printf("usage:./get_signifAF name.thr input2 output1 output2\n");
        return -1;
    }


    thrFile=fopen(argcv[1],"r");
		    if (thrFile == NULL) {
		      printf("cannot open first input _threshold file %s\n", argcv[1]);
		      return -1;
    }

    extract_vcfFile=fopen(argcv[2],"r");
    if (extract_vcfFile == NULL) {
      printf("cannot open first input _vcf.txt file %s\n", argcv[2]);
      return -1;
    }

    afFile=fopen(argcv[3],"w");
    if (afFile == NULL) {
    printf("cannot open first output1 infoAF.txt file %s for writing\n", argcv[3]);
    return -1;
    }

    distribFile=fopen(argcv[4],"w");
    if (distribFile == NULL) {
    printf("cannot open second output file distAF %s for writing\n", argcv[4]);
    return -1;
    }

    printf("get_peak_thr\n");


    // initiate zero vector for histograme
			 for(i=0;i<100;i++)//
	         {
	              histAF[i]=0;
		     }
    // initiate zero vector for threshold
			 for(i=0;i<1;i++)//
			 {
			 	  THR[i]=0;
		     }

		// scan thr file	and assign current precomputed thresholds for filtering
    while( (n = fscanf(thrFile,"%d %d", &thR, &thA )) >= 0)
         // until the end of the input threshold  file
    {
		      if (n != Nthr) // incorrect format, NCOLS=number of columns in pipeCT input
		      {
		      printf ("corrupted input thrFile format\n");
		      return -1;
		      }

		THR[0]=thR;
		THR[1]=thA;
	}


	//6 fields of input file	// no PV4
//CSV reading...fine!!!
    while( (n = fscanf(extract_vcfFile,"%d,%d,%d,%d,%d,%d", &pos, &D, &DP4[0], &DP4[1], &DP4[2], &DP4[3])) >= 0 && firstSixAreOK == 1 && canWriteAF == 1)
    // read the Read Title
    {
        if( n != NRID )     // incorrect format
        {
            firstSixAreOK = 0;
            break;
        }

           count_beforeF++;
           sumDbefore=sumDbefore+D;

           // f1 f2 Josie : filter for DP4 separately for ref and alternative alleles

            if ( DP4[0] >= THR[0] && DP4[1]>= THR[0] && DP4[2]> THR[1] && DP4[3]> THR[1])// 0 ref are ok
           {
                AF = GetAf(DP4);
                histAF[AF-1]++;
		        // compute average run depth default (Q13), D across pos,sumD
				// compute average run depth after Q30, sum(DP4) across pos,sumDP4
				        count_afterF++;
				        sumD=sumD+D;
				        sDP4=DP4[0]+DP4[1]+DP4[2]+DP4[3];
				        sumDP4=sumDP4+sDP4;
				        Q30frac=sDP4/D;
                 //printf("frac filtered out Q= %.2f\n", Q30frac);

                // output check: 2 columns position and its AF
                if( fprintf(afFile,"%d %d %.2f\n", pos, AF, Q30frac) < 0 )
                 {
                 canWriteAF = 0;
                 break;
                 }

             }// end filter DP4



        if (canWriteAF == 0) {break;}

        // removing space symbols before carriage return
        do  {
            n = fgetc (extract_vcfFile);
        }
        while ((char)n != '\n');

    } //END of outer while loop

    // to find significat peak of AF histo!
        mu=GetMu(histAF);
        std=GetStd(histAF);

      for(i=0;i<49;i++)//only 49, because we need only lower freq, esp for diploid organisms
      {
		x=(histAF[i]-mu)/std;

		if (x>=0)
        Z= (long)( x +0.5);

        if (x<0)
        Z= (long)( x -0.5);

        //printf("%d\n", Z);
        ZZ[i]=Z;
		 if( fprintf(distribFile,"%d %d %.2f\n",i, histAF[i],x) <= 0 ) {
         //if( fprintf(distribFile,"ae = %3.5f\tAF_avr = %6.2f\tsigma = %6.2f\n",
         // (float)cntZero/(float)counter, avrAF, sigma ) <= 0 ) {
        canWriteDistrib = 0;
	     }

       }//for

      // find all percents with non-zero histAF, and then MIN among them
       //perc=GetMinPercAF(histAF);
        // perc=GetSignifPercAF(histAF,ZZ);

       // find signid peak in low freq histogram as a mixture
         perc=GetPeakPercAF(histAF);
       // conition on mixture percent: if more than 25, it is likely to be diploid allele  :  17 Nov
       if (perc < 25)
       {
		   printf("percent_mix= %d\n", perc);
	   }
		if (perc > 25 || perc == 25)
       {
		   printf("not likely there is a low percent mixture here at this min depth\n");
        }

       // compute avarege depth :default and after filtering
      avD=sumD/count_afterF;
      avDP4=sumDP4/count_afterF;
       avDbefore=sumDbefore/count_beforeF;
       perc_left=(100.0*count_afterF/count_beforeF);

       conf=(1-1/avDP4)*0.8;//(xx[perc]+1)*(1-1/avDP4);// 0.8 is MY constant

	   printf("Confidence of mixture value= %.2f\n", conf);


       //printf("average Depth default after min_depth filter= %.2f\n", avD);
         printf("average actual Depth after Q25-filter and min_depth filter= %.2f\n", avDP4);
           printf("average Depth before filtering= %.2f\n", avDbefore);
    printf("percent data left after min_depth filtering= %.2f\n", perc_left);


    fclose(distribFile);
    fclose(extract_vcfFile);
    fclose(afFile);
    //free(outFileName);


    if( firstSixAreOK  == 0 ||
        canWriteAF == 0 || canWriteDistrib ==0) {
        printf ("Error during execution. Details: \n");
        printf ("\tfirstSixAreOK %d\n",  firstSixAreOK);
        printf ("\tcanWriteAF %d\n",        canWriteAF);
        printf ("\tcanWriteDistrib %d\n",      canWriteDistrib);
        printf ("Execution aborted\n");
        return -1;
    }

    printf(" done.\n");
    return 0;
}//main

/***************************************************** Functions******************************/


////////////////////////////////////////////////////
// Calculates AF=alternative allele frequebcy from data=DP4 counts
// input - array (4-vector) of DP4 counts
//output -one float number from (0,1)
////////////////////////////////////////////////////
int GetAf (int data[])
{
    float Isum,AF1;
    int i,AF;

    Isum = data[0];
    for (i=1; i<4; i++)
    {
        Isum += data[i];
    }
    AF1=(data[2]+data[3])/Isum;
    AF=(long) (100*AF1+0.5);
    //(long) (x+0.5)
    return AF;
}
// GetMu.c calculates mu for first 50 non-zero elements of 100-vector
float GetMu( int data[])
 {
	 float mu,std,Isum,Isum2;
	 int i;
     int cnz=0;

     Isum = data[0];
     Isum2 = data[0]*data[0];
	 for (i=1; i<49; i++)
	 {
		 if (data[i] > 0)
		 {
		 cnz++;
	     Isum += data[i];
	     Isum2 += data[i]*data[i];
	     }
      }
      mu=Isum/cnz;
      //std=Isum2/cnz-mu*mu;
      return mu;
}

// stand deviation
float GetStd( int data[])
 {
	 float mu,std,Isum,Isum2;
	 int i;
     int cnz=0;

     Isum = data[0];
     Isum2 = data[0]*data[0];
	 for (i=1; i<49; i++)
	 {
		 if (data[i] > 0)
		 {
		 cnz++;
	     Isum += data[i];
	     Isum2 += data[i]*data[i];
	     }
      }
      mu=Isum/cnz;
      std=sqrt(Isum2/cnz-mu*mu);
      return std;
}

// min peak percents with non-zero histAF
int GetPeakPercAF(int data[])
{
//consider only first 49

// find minimal (first) non-zero peak
	//data is histAF, its index is bin=percentage


	int i;
	int perc=0;

	for (i=1; i<48; i++)
	 {
	  if ( ((data[i]-data[i-1]) >0 ) && ((data[i+1]-data[i]) <0))//here peak in [i]
	  {
	   perc=i;
	   break;
      }
     }
     return perc;
 }






















