 /*  File: get_mixture_hap13.c // computes and analyses AAF( maf) distribution over variant positions from vcfq files
 * Authors: designed by Irina Abnizova (ia1)edited by Steve Leonard (srl)
 *

  Last edited:
  13 June-output sk25,sk75
  12 june-- final mixture file; modified as output name.mix:  mix_percent=17
possible_mix_percent= 17
confidence_nonzeroMix=0.64
AvActDepth=52 (avDP4)
min_depthR=2
min_depthA=2

  30 May last 75 implement
  29 May  mixture file instead of infoAF file
  29 May tuning params sk25,sk75, n25 n49 n51 n75
  23 May: filtered on 'bad' regions/variants: if their D>mu+std
  7 May-functions to get skew25
  6 May - do specially for human 20x-diploid, and haploid thr_sk=0.8

  25 April: incorporate info about (100-mix)% to increase proly of detection
  24 April-output in infoAF 4 fields for each variant position, last is its actual coverage sDP4
  22 April: max sginif or first max if not signif  + all 100 AAF to see 'complimentary'
  10 April:   if ( DP4[0] >= THR[0] && DP4[1]>= THR[0] && DP4[2]>= THR[1] && DP4[3]>= THR[1])// 0 ref alt are ok
              thr_sk
  10 March 2014 for skewness to first 25% AAF
  3 feb-reads thr from a file; 29 Jan 2014  - get autom thr, starting with default ones 1,5
  still old fields of vcf input
  21 Nov 2013, for two thresholds

 *-------------------------------------------------------------------
 * Description: computes maf distrib over variant positions from just extracted 6 fileld

 * Exported functions:
 * HISTORY:

skewness to first 25% of AAF  ----------introduced 10 March
computes AF (alternative to Ref allele frequency) for each error base call from vcf

applies minimum depth threshold, both for Ref and Alternative alleles fwd rev (total min depth will be
4*MIN_DEPTH

input1 thrfile  (thRR thAA)
input2: extracted vcf,


output1: mix_conf AF.txt=
mixture confidence
12        0.87

output2: .distrAF=histogram of AF with significance for each bin
% count
8 0
9 0
10 385
11 780


usage:./get_mixture rams.thr extracted.vcf .mix .distrAF

*/

#include <stdio.h>
//#include "conio.h"
#include <time.h>
#include <math.h>
#include <string.h>

// ******** predefined user constants, 22 october

#define NRID 6      // Number of variant ids=fields: pos, chromosome, Depth,DP4,(no PV4+maf) in esigtracted_vcf files
#define NPAR 5      // masig Number of PARameters and arguments
#define Nthr 4     //  Number of thresholds and muD1 stdD1 12 June


// ******** declarations of  functions
int   GetAf (int []);// computes AF percent for a variant position
float   GetMu( int data[]);
float   GetStd( int data[]);
float   GetSkew25( int data[],int n25, int n49);
float   GetSkew75( int data[],int n51, int n75);
float GetLikely25(int perc, float sk25);

float GetConfidence( float avDP4,float likely25,float likely75, float sk25, float sk75,float thr25,float thr75);
int MakeDecisionMix(int perc, float sk25, float thr25,float sk75,float thr75);// update percentage mix

int GetMaxPeak25(int data[], int n25);
float ComplementaryPeaks75(int data[], int perc, float sk75);// sh be around=3 etc

//=====================================================MAIN

int main    (int argc, char *argcv[]) {

    // flags
    int firstSixAreOK  = 1;
    int canWriteMix        = 1;
    int canWriteDistrib      = 1;

    FILE *extract_vcfFile, *mixFile, *distribFile, *thrFile; //file handles

        int n,i;

    //from input vcfq
	    int DP4[4];//
	    int D,pos;// Depth, genome pos of error

	    //stats
	    int avDP4,stdDP4;// of actual sDP4 depth after RA & bad regions filtering
	    int sumDP4=0,sumDP4_2=0;//across all pos for mu std
	    int sumDbefore=0;
	    float avDbefore;
		float Q30frac;
		float sDP4;
		float perc_left;

      // from rams(now four values thrR thrA muD stdD
       int muD,stdD;
       int thR,thA;


       // settings for skewness and their likelihoods
       int n25,n49;
       int n51,n75;
       float thr25,thr75;

      //what we count
       int count_afterF=0;// after filtering with min depths
	   int count_beforeF=0;// before filtering with min depths

       int perc,perc_update;
       int AF;
       float conf;//confidence (0,1)
       float sk25,sk75;//skewness to first 25% of AAF  ----------introduced 10 March
       float likely75,likely25;

       //arrays
       int histAF[100];// to store mafs percentages

       int THR[2];// compute from input1 thr file
       int msD[2];


    if(argc < NPAR)//four input_output files are submitted
    {
        printf("not enough of parms |input_output files\n");
        printf("usage:./get_mixture name.thr input2 output1 output2\n");
        return -1;
    }


    thrFile=fopen(argcv[1],"r");
		    if (thrFile == NULL) {
		      printf("cannot open first input _threshold file %s\n", argcv[1]);
		      return -1;
    }

    extract_vcfFile=fopen(argcv[2],"r");
    if (extract_vcfFile == NULL) {
      printf("cannot open first input _.vcfq file %s\n", argcv[2]);
      return -1;
    }

    mixFile=fopen(argcv[3],"w");
    if (mixFile == NULL) {
    printf("cannot open first output1 mix_conf file %s for writing\n", argcv[3]);
    return -1;
    }

    distribFile=fopen(argcv[4],"w");
    if (distribFile == NULL) {
    printf("cannot open second output file distAF %s for writing\n", argcv[4]);
    return -1;
    }

    printf("get_mixture_haploid\n");

    // ===SETTINGS     parameters for skewness measurements and peak consideretion Diploid! (more noise around 50%)
            n25=25;
            n49=40;
            n51=70;
            n75=83;
       // give value for thr25=    0.55 for human  diploid; 1.0, 0.8 for pathogens-haploid
            thr25=0.8;
            thr75=0.8;//complementary skewness threshold

    // initiate zero vector for histograme
			 for(i=0;i<100;i++)//
	         {
	              histAF[i]=0;
		     }
    // initiate zero vector for threshold
			 for(i=0;i<1;i++)//
			 {
			 	  THR[i]=0;
			 	  msD[i]=0;
		     }

		// scan thr file	and assign current precomputed thresholds for filtering
    while( (n = fscanf(thrFile,"%d %d %d %d", &thR, &thA,&muD, &stdD)) >= 0)
         // until the end of the input threshold  file
    {
		      if (n != Nthr) // incorrect format input
		      {
		      printf ("corrupted input thrFile format\n");
		      return -1;
		      }

		THR[0]=thR;
		THR[1]=thA;
		msD[0]=muD;
        msD[1]=stdD;
	}

	       // printf("recommended precomputed Ref threshold= %d\n", thR);
	       // printf("recommended precomputed Alternative threshold= %d\n", thA);
	       // printf("actual depth of input vcfq after filtering RA 1,1= %d\n", muD);
            //printf("std of actual depth including possible bad regions= %d\n", stdD);

	//6 fields of input file	// no PV4: doit to compute AAF=AF

    while( (n = fscanf(extract_vcfFile,"%d,%d,%d,%d,%d,%d", &pos, &D, &DP4[0], &DP4[1], &DP4[2], &DP4[3])) >= 0 && firstSixAreOK == 1 && canWriteMix == 1)
    // read the Read Title
    {
        if( n != NRID )     // incorrect format
        {
            firstSixAreOK = 0;
            break;
        }
           count_beforeF++;
           sumDbefore=sumDbefore+D;

         // f1 Josie : filter for DP4 separately for ref and alternative alleles+ filter abnormal Depth

            if ( DP4[0] >= THR[0] && DP4[1]>= THR[0] && DP4[2]>= THR[1] && DP4[3]>= THR[1] && D < (msD[0]+msD[1]))// 0 ref alt are ok
           {
                AF = GetAf(DP4);
                histAF[AF-1]++;

		        // for EACH POSITION, compute depths after thrRA filtering:
                       // average run depth after Q30, sum(DP4) across pos,=  sumDP4
                        count_afterF++;
				        sDP4=DP4[0]+DP4[1]+DP4[2]+DP4[3];
				        sumDP4=sumDP4+sDP4;
				        sumDP4_2=sumDP4_2+sDP4*sDP4;//for std
				        Q30frac=sDP4/D;// how actaul quality Depth differs from default
             }// end min_depth filter of DP4

       // removing space symbols before carriage return
        do  {
            n = fgetc (extract_vcfFile);
        }
        while ((char)n != '\n');

    } //END of outer while loop across vcfq file   AF is computed!

    //====================== output2    write AF
      for(i=0;i<101;i++)
      {
		  if( fprintf(distribFile,"%d %d\n",i, histAF[i]) <= 0 ) {
          canWriteDistrib = 0;
	      }

       }//  for i=100 cycle
//=======================================analysis of AF distribution to find mixture sample

         sk25=GetSkew25(histAF,n25,n49);
         sk75=GetSkew75(histAF,n51,n75);

// find  what percent AAF gives first max peak among first 25% AAF

             perc=GetMaxPeak25(histAF,n25);//peaks in first 25% AAF
             likely25= GetLikely25(perc,sk25);
             //printf("likelihood of having low AAF peak= %.2f\n", likely25);

      // stats after filtering :compute average depth :default and after filtering


             avDP4=ceil(sumDP4/count_afterF);//actual average cov  after Filt
	         stdDP4=ceil(sqrt(sumDP4_2/count_afterF));

	         avDbefore=sumDbefore/count_beforeF;
	         perc_left=(100.0*count_afterF/count_beforeF);


	         //printf("average default Q13 Depth before min_depth filtering= %.2f\n", avDbefore);
	         //printf("average actual Depth (avDP4) after min depth and bad regions filter= %d\n", avDP4);
	         //printf("std actual Depth (stdDP4) min depth and bad regions filter= %d\n", stdDP4);
	         //printf("percent data left after min depth and bad regions filter = %.2f\n", perc_left);

	         // output percent of bad regions!

     // impact of (100-%mix) AAF : factor probability
            likely75=ComplementaryPeaks75(histAF, perc,sk75);
            //printf("likelihood of having complementary AAF peak= %.2f\n", likely75);

     // make decision about on mixture percent: update it to zero if needed

	        perc_update=MakeDecisionMix(perc, sk25, thr25,sk75,thr75);

      // confidence of mixture estimating

           conf=GetConfidence( avDP4,likely25,likely75,sk25,thr25,sk75,thr75);
           printf("updated percent mixture= %d\n", perc_update);
         //printf("possible percent mixture= %d\n", perc);
           printf("Confidence of non-zero mixture = %.2f\n", conf);
	     //printf("skewness of first 25 perc AF = %.2f\n", sk25);
	     //printf("skewness of last 25 perc AF = %.2f\n", sk75);

	   // MAIN output
	  //if (fprintf(mixFile,"mix_percent=%d\tpossible_mix_percent= %d\tconfidence_nonzeroMix=%.2f\tAvActDepth=%d\tmin_depthR=%d\tmin_depthA=%d\n", perc_update,perc,conf,avDP4,thR,thA)<= 0 )
	  if (fprintf(mixFile,"mix_percent=%d\npossible_mix_percent= %d\nconfidence_nonzeroMix=%.2f\nAvActDepth=%d\nmin_depthR=%d\nmin_depthA=%d\nskew25=%.2f\nskew75=%.2f\n", perc_update,perc,conf,avDP4,thR,thA,sk25,sk75)<= 0 )

	  {
	                canWriteMix = 0;
	                //return -1;
          }
    fclose(distribFile);
    fclose(extract_vcfFile);
    fclose(mixFile);

    if( firstSixAreOK  == 0 ||
        canWriteMix == 0 || canWriteDistrib ==0) {
        printf ("Error during execution. Details: \n");
        printf ("\tfirstSixAreOK %d\n",  firstSixAreOK);
        printf ("\tcanWriteMix %d\n",        canWriteMix);
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
    for (i=1; i<4; i++)//there are four conts for each base
    {
        Isum += data[i];
    }
    AF1=(data[2]+data[3])/Isum;
    AF=(long) (100*AF1+0.5);
    //(long) (sig+0.5)
    return AF;
}
// GetMu.c calculates mu for all non-zero elements of 100-vector
float GetMu( int data[])
 {
	 float mu,Isum;
	 int i;
     int cnz=0;

     Isum = data[0];
     //Isum2 = data[0]*data[0];
	 for (i=1; i<101; i++)
	 {
		 if (data[i] > 0)
		 {
		 cnz++;
	     Isum += data[i];
	     //Isum2 += data[i]*data[i];
	     }
      }
      mu=Isum/cnz;
      //std=Isum2/cnz-mu*mu;
      return mu;
}

// ========================stand deviation
float GetStd( int data[])
 {
	 float mu,std,Isum,Isum2;
	 int i;
     int cnz=0;

     Isum = data[0];
     Isum2 = data[0]*data[0];
	 for (i=1; i<101; i++)
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


//----------------------max from the first 25
 int GetMaxPeak25(int data[], int n25)
 {

  // input data=AAF
	 int Imax;
     int i,perc,k;
     int peak[n25], per_peak[n25];

     // initialise perc

     perc=0.0;

// find max peak among peaks, first 25 percent only, data is maf(AAF) vector

//what if no peaks?
     peak[0]=0;


     // find peaks
     k=0;
     for (i=1; i<n25; i++)
     {
		 if ( ((data[i]-data[i-1]) >0 ) && ((data[i+1]-data[i]) <=0))// && ( Imax < data[i]))//here peak in [i]
         {

             peak[k] = data[i];
             per_peak[k]=i;
             //printf("peak = %d\n", peak[k]);
             //printf("per_peak = %d\n", per_peak[k]);
             k=k+1;

         }
	 }
// find max peak/s if any exist

  if (peak[0]>0)
  {
	 Imax = peak[1];

	 	      for (i=1; i<k; i++)
	 	      {
	 	          if( Imax < peak[i])
	 	              Imax = peak[i];
	 	              //perc=i;
	           }

       //printf("max peak = %d\n", Imax);
   // find percent giving Imax, first one if several are max
       for (i=0; i<k; i++)
  	   {
  	          if( peak[i]==Imax)
  	            { perc=per_peak[i];

  	              break;}
       }
   }// if (peak[1]>0)-exist any peak <25%

       return perc;

}
// compute likely25 based on existence of peaks, which is  perc=0, and sk25

float GetLikely25(int perc, float sk25)
{
     float likely25;

    if (perc<=0) likely25=0.0;
    //display('not likely any low AAF mixture');
     if (perc>0)
     {
       if (sk25<=1) likely25=sk25*sk25;
       if (sk25>1) likely25=1+(1-1/sk25);
    }

    return likely25;
}

//===================================
float   GetSkew25( int data[],int n25, int n49)
 {

	 float sk25;
     int i,k;
     int less25;
     int bw25_50;
           // sum up from 0 to n25 %
           less25=data[0];
	 	   for(i=0;i<n25;i++)
	 	   {
	 	   		less25+=data[i];
	 	   }

	 	   bw25_50=data[n25];
	 	   for(i=n25;i<n49;i++)
	 	   {
	 	   	    bw25_50+=data[i];
	 	   }

	   //compute 25% AF skewness, sk25

	 	   sk25=0.00;// if no bases for bw25_50

	 	   if (bw25_50 >0 )
	 	   {
	 	   sk25=1.0*less25/bw25_50;// if proportion of less 25% AF is high (more 1?), more likely to be a mixture
	        }
	  return sk25;
}
//  ---------------------float   GetSkew( int data[],int n51, int n75);
float   GetSkew75( int data[],int n51, int n75)
 {

	 float sk75;
     int i,k;
     int more75;
     int bw75_50;
     int last;

     last=5;
           // sum up from n75 to almost 100%
           more75=data[n75];
	 	   for(i=n75;i<(100-last);i++)
	 	   {
	 	   		more75+=data[i];
	 	   }
           // sum 50 to 75
	 	   bw75_50=data[n51];
	 	   for(i=n51;i<n75;i++)
	 	   {
	 	   	    bw75_50+=data[i];
	 	   }

	   //compute 75% AF skewness, sk75

	 	   sk75=0.00;// if no bases for bw75_50

	 	   if (bw75_50 >0 )
	 	   {
	 	   sk75=1.0*more75/bw75_50;// if proportion of less 75% AF is high (more 1?), more likely to be a mixture
	        }
	  return sk75;
}

//===================== confidence of mixture estimation
float GetConfidence( float avDP4, float likely25,float likely75, float sk25, float sk75,float thr25,float thr75)

{
        float conf;
     // confidence of mixture estimating

       conf=(1-1/avDP4)*0.6;

       if ((sk25 >= thr25) && (sk75 >= thr75))// larger any how!
       {conf=(1-1/avDP4)*0.7;//(sigsig[perc]+1)*(1-1/avDP4);// 0.8 is MY arbitrary constant
       }

         //likely75=ComplementaryPeaks75(histAF, perc); and likely25 impacts

         conf=conf + likely75*likely25;

         if (conf>1)
         {conf=1.0;}

  return conf;
}
//float GetConfidence( float avDP4, float likely75, float sk25, float sk75,float thr25,float thr75);

// make decision if perc noise or mixture: update it to zero if needed
//if perc is more than n25 (here =25), it is likely to be diploid allele noise  :  17 Nov+ amount AF<25%  skewed (too high)
//int MakeDecisionMix(int perc, float sk25, float thr25)
int MakeDecisionMix(int perc, float sk25, float thr25,float sk75,float thr75)// update percentage mix

{
       int perc_update;
       perc_update=0.0;
       if ((sk25 >= thr25) && (sk75 >= thr75))// thre os mixture=perc ;thr25 was 0.95 for path, 0.62 for human
       {
		   perc_update=perc;
		   //printf("percent mixture= %d\n", perc);
	   }


        return perc_update;
}

//===================================
//=================================================finds AAF complementary peaks, (100-perc)AAF if any

float ComplementaryPeaks75(int data[], int perc,float sk75)
{
      float likely75;
      int i,k,n,around; ;// frequency may vary up to this amount of %
      int peak[25], per_peak[25];


// find max peak among peaks, last 75% (because as mixture was first 25 percent only)
     // find peaks
    k=0;
    for (i=75; i<101; i++)
    {
        if ( ((data[i]-data[i-1]) >0 ) && ((data[i+1]-data[i]) <=0))
        {
			 peak[k] = data[i];
             per_peak[k]=i;
             k=k+1;
         }
     }
	 // see if these peaks are complementary to detected perc=% mixture <25%
	 n=0;
     if (k<=0) likely75=0;// no peaks
     if (k>0)// exist any peaks after 75
	 {
           for (i=0; i<k; i++)
	       {
			   if ( (per_peak[i] > (100-perc-around)) && (per_peak[i] < (100-perc+around)))  n=n+1;
		    }
      }// if k>0

      if (n>0)// exist complementary peaks after 75
      {
       if (sk75<=1) likely75=sk75*sk75;
       if (sk75>1) likely75=1+(1-1/sk75);
      }
       return likely75;
}
