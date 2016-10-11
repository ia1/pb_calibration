/*  File: compare_samples.c // compares multiple AAF( maf) distributions and mixture probabilities
 * Authors: designed by Irina Abnizova (ia1)
 *

  Last edited:
 10 october : Fail/pass per tag based on combined info about mixture likelihood and outlier's AAF distribution in a lane across tags
 7 october- extra output: Final Table per tag
 6 October- Normalise all histos up to percentages
 3 october 2016- created the first version

 *-------------------------------------------------------------------
 * Description: compares AAF distrib over variant positions across tags, and finds outliers
   compares with possible percentage of contaminated mixture (three modes) and its confidence

 * Exported functions:
 * HISTORY:

*/

#include <stdio.h>
#include <errno.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <string.h>
#include <getopt.h>

// ******** predefined user constants,
#define NPAR 4                // Number of arguments int1 int2 out1
#define NF 6                  // number of fields in mix file
#define N_perc  101          // Number of columns =percentages AAF for histogram per tag
#define MAX_LINE_NUM 999  // Max line=max tags

// ******** declarations of  functions
//void usage(int code); // usage
float   GetMu( int data[]);
float   GetStd( int data[]);
int SmoothWin2(int val1,int val2);
void Smoothing_iterate( int histS[], int hist[] );
void Normalise_percentage( int hist[] );
int SmoothWin3(int val1,int val2,int val3);
int GetMaxIndex (float data[], int i1, int i2);
float GetMax(float data[], int i1, int i2);
float Decision(float Lik[], float mix[], float thr_lik);

////////////////////////////////////////////////////
// main
////////////////////////////////////////////////////
//int main(int argc, char *argv[]) {
int main(int argc, char *argcv[])
{


    int rrow;
    int tag;
    int nbins;
    int wid;

    char buffer[10024] ;// declared 10x larger char array! was [1024] before
	char *record,*line1;
	int i=0,j=0;

    int mat[1000][101];
    int hist[101];
    int histS[101];
    int matS[1000][101];// smoothed matrix

    int sum_hist[101];// sum of hists across (smoothed) tag histos
    int av_hist[101];// av of hists across (smoothed) tag histos

    int diff_up37[1000];// differences b/w 'model'=averaged and tags shoulder as sum 0:37 bins; per tag
    int diff_after63[1000];// differences b/w 'model'=averaged and tags shoulder as sum 63:95 bins
    int diffs[1000];// sum up differences b/w 'model'=averaged
    int FP[1000];// Fail or pass (0/1) per tag


    float mu;// for diffs
    float sd;

     //---------------results to compute:
	int count_tags; // tag counts
	float mixture; // mixture
	float max_conf;

	    //---------------------------2D arrays to read into memory to operate
	    float mix_Lik[1000][6]; // mix likelihood (or confidences) three modes, per tag
	    float mix_conf[1000][2];// mix ,  confidences max mode, per tag

	    float mix[3];
	    float Lik[3];

    // file names
    char *mix_file, *distrib_file, *fintab_file;

    // file handles
    FILE *mixFile, *distribFile, *fintabFile;

    //params
    float thr_lik;//how low sh be likelihood to ignore it; was 0.25, now 0.3

    // to read vcfa/vcfq
    static const int line_size = 8192; // maximum line size
    char line[line_size];

    // open files for read/write

      // ----------------------------------------check arguments/parameters
	    if(argc < NPAR)//four input_output files are submitted
	    {
	        printf("not enough of input_output files\n");
	        printf("usage:./compare in1 in2\n");
	        return -1;
	    }

    mixFile = fopen(argcv[1],"r");// make read
    if (mixFile == NULL) {
        fprintf(stderr, "cannot open mixture_file %s: %s\n", mix_file, strerror(errno));
        exit(EXIT_FAILURE);
    }
    distribFile = fopen(argcv[2],"r");// make read
    if (distribFile == NULL) {
        fprintf(stderr, "cannot open distribution_file %s: %s\n", distrib_file, strerror(errno));
        exit(EXIT_FAILURE);
    }
     fintabFile = fopen(argcv[3],"w");//
	    if (fintabFile == NULL) {
	        fprintf(stderr, "cannot open fintab_file %s: %s\n", fintab_file, strerror(errno));
	        exit(EXIT_FAILURE);
    }


      fprintf(stderr, "compare_tags\n");


	      // param values
	      thr_lik=0.3;
	      nbins=100;
	      wid=1;
	      //bins[0]=wid;

        tag=0;
    ///1 input========================= read from distrFile into 2D histAAF array

while((line1=fgets(buffer,sizeof(buffer),distribFile))!=NULL)
   {
     record = strtok(line1,",");
     while(record != NULL)
     {
     //printf("record : %s", record) ;    //here I can put the record into the array .
     mat[i][j++] = atoi(record) ;
     record = strtok(NULL,",");
     }
     ++i;
     rrow=i;
   }
 printf("\nnumber of rows in long histo file = %d\n",rrow);


 //input2================================== read mix - 6 fields into mix_Lik array
    tag=0;
    while (fgets(line, line_size, mixFile))
    {
        int k = sscanf(line, "%f,%f,%f,%f,%f,%f", &mix[0], &mix[1], &mix[2], &Lik[0], &Lik[1], &Lik[2]);
        if (k != NF)
        {
            // number of fields read not correct
            //fprintf (stderr, "corrupt mix file %s\n", line);
             fprintf(stderr, "skipping malformed mix line %s", line);
            //exit(EXIT_FAILURE);
            continue;//skip
        }

        mix_Lik[tag][0]=mix[0];
        mix_Lik[tag][1]=mix[1];
        mix_Lik[tag][2]=mix[2];
        mix_Lik[tag][3]=Lik[0];
        mix_Lik[tag][4]=Lik[1];
        mix_Lik[tag][5]=Lik[2];


        mixture = Decision(Lik, mix, thr_lik);// !!!!
        max_conf=GetMax(Lik,0,3);

	    fprintf(stderr, "final decided mixture = %.2f\n",  mixture);
        mix_conf[tag][0]=mixture;
        mix_conf[tag][1]=max_conf;
         fprintf(stderr, "max conf = %.2f\n",  mix_conf[tag][1]);

         tag++;
      }

    printf("\nnum of lines-tags in mix_LIk file %d\n",tag);


//-----------starting to smooth histos
    for (j=0; j<tag; j++)
    {
          // initiate histS for each tagj
	      for (i=0; i<101; i++)
	      {
		         hist[i]=mat[j][i];
		         histS[i]=mat[j][i];// initialised
		  }

		  //----------------------Normalise both

		  Normalise_percentage( hist );
		  Normalise_percentage( histS);

		  // smooth hist
             Smoothing_iterate( histS, hist);

          for (i=0; i<101; i++)
          {
              matS[j][i]=histS[i];
		  }

  }


   //-----------starting to summing up n(already Normalised to percentage)

     //initialise sum_hist as zeros
     for (i=0; i<101; i++)
   	      sum_hist[i]=0;

       for (j=0; j<tag; j++)
       {
   	      for (i=0; i<101; i++)
   	      {
           sum_hist[i]=sum_hist[i]+matS[j][i];
           av_hist[i]=(int) (sum_hist[i]/tag+0.5);

		  }
       }

    printf("differences with the model\n");

   //diff_up37[1000]
   // initialise diff for each tag
   for (j=0; j<tag; j++)
   {
        diff_up37[j]=0;
        diff_after63[j]=0;
   }
//----------diff 37
        for (j=0; j<tag; j++)
		{
		   	      for (i=0; i<37; i++)
   	              {
			      diff_up37[j]=diff_up37[j]+matS[j][i]-av_hist[i];// shoulders are wider for contaminated
			      }
         //printf("diff 37 = %d\n",diff_up37[j]);
        }
//------------------diff after63

        for (j=0; j<tag; j++)
		{
		   	      for (i=63; i<96; i++)
   	              {
			      diff_after63[j]=diff_after63[j]+matS[j][i]-av_hist[i];// shoulders are wider for contaminated
			      }
         //printf("diff 63 = %d\n",diff_after63[j]);
        }

// sum up and find mu and std for each of diff
        for (j=0; j<tag; j++){
                 diffs[j]=diff_after63[j]+diff_up37[j];
                //printf("diff = %d\n",diffs[j]);
			 }

        mu=GetMu( diffs);
        sd=GetStd( diffs);

        //printf("mu diff = %.2f\n",mu);
        //printf("std diff = %.2f\n",sd);

        // combine information about tags mixture confidence and its local distrAAF (outlier?...)

         for (j=0; j<tag; j++){
			 FP[j]=1;
			 if ((diffs[j] > mu+sd) && ( mix_conf[j][1] > thr_lik))
			 FP[j]=0;
			 printf("fail or pass = %d\n",FP[j]);
		 }


			//  ====================================oUTPUT1

  for (i=0;i<tag;i++) {

		    if (fprintf(fintabFile,"%d %.2f %.2f %d %d\n",i, mix_conf[i][0],mix_conf[i][1],diffs[i],FP[i]) <= 0) {
                fprintf(stderr, "error writing final table file: %s\n", strerror(errno));
            exit(EXIT_FAILURE);
        }
    }

    //FILE *mixFile, *distribFile, *fintabFile;

    fclose(distribFile);
    fclose(fintabFile);
    fclose(mixFile);

    fprintf(stderr, "done compare tags AAF distributions \n");

return 0;
}

//==============================================functions

///==============================================
// first initiate histS[]=0  //num_peS = num_peI;
void Normalise_percentage( int hist[] )
{

   int i;
   int s;

   // sum up histo
   s=0;
    for (i=0;i<101;i++)
    	s=s+hist[i];

    	if (s > 0)
    	{
			for (i=0;i<101;i++)
			hist[i]=(int) ((1000.0*hist[i]/s)*0.1 +0.5);
	    }
}
///==============================================
// first initiate histS[]=0  //num_peS = num_peI;
void Smoothing_iterate( int histS[], int hist[] )
{

    //diff = 1;
    int k,i;
    //int num_peaks;

    k=0;
    //while ((diff>0))
    for (k=0;k<4;k++)// 4 smoothings
    {
        //float histS[100];

        //---------------smooth given hist
        histS[0] = SmoothWin2(hist[0],hist[1]);
        for (i=1;i<100;i++)
        {
            histS[i] = SmoothWin3(hist[i-1],hist[i],hist[i+1]);
        }
        histS[100] = SmoothWin2(hist[100-1],hist[100]);


        // ------------re-assign current hist as histS
        for (i=0;i<101;i++)
        {
            hist[i] = histS[i];
        }
        //k++;// number of smoothing here
    } //end while; for

    //printf("number of smoothings= %d\n", k);
}

///----------smoothing in a 2 bin window
int SmoothWin2(int val1,int val2)
{
    int hi2;

    hi2 = (int)((val1 + val2)/2.0+0.5);
    return hi2;
}
int SmoothWin3(int val1,int val2,int val3)
{
    int hi3;

    hi3 =(int)( (val1 + val2 + val3) / 3.0+0.5);

    return hi3;
}
////////////////////////////////////////////////////
// find position of the max value in a sub-interval (i1<=i< i2) of an array
////////////////////////////////////////////////////
int GetMaxIndex (float data[], int i1, int i2)
{
  float data_max;
  int i, ind;

  data_max = data[i1];
  ind = i1;
  for (i=i1; i<i2; i++) {
      if( data_max < data[i]) {
          data_max = data[i];
          ind = i;
      }
  }

  return ind;
}

////////////////////////////////////////////////////
// find max value in a sub-interval (i1<=i< i2) of an array data[i]
////////////////////////////////////////////////////
float GetMax(float data[], int i1, int i2)
{
  float data_max;
  int ind, i;

  data_max = data[i1];
  ind = i1;
  for (i=i1; i<i2; i++) {
      if( data_max < data[i]) {
          data_max = data[i];
          ind = i;
      }
  }

  return data_max;
}

////////////////////////////////////////////////////
// decision based on confidences modes 0,1,2
// we assume that the larger is skewness of a mode (pronounced humph), the more likely mixture belongs there
////////////////////////////////////////////////////

float Decision(float Lik[], float mix[], float thr_lik)
{
    float mixx;
    int nm;

     // if there is no peak at mode0, then take its middle
	    if (mix[0]==0.0)
			mix[0]=6.0;

    // if each Lik[i] <=thr_lik, mixture=0
    if (Lik[0] <= thr_lik && Lik[1] <= thr_lik && Lik[2] <= thr_lik){
        mixx = 0.0;
    }
    else {
        nm = GetMaxIndex(Lik, 0,3);
        mixx = 1.0*mix[nm];
    }

    return mixx;
}

// GetMu.c calculates mu for all non-zero elements of 100-vector
float GetMu( int data[])
 {
	 float mu,Isum;
	 int i;
     int cnz=0;

     Isum = data[0];
     for (i=1; i<101; i++)
	 {
		 if (data[i] > 0)
		 {
		 cnz++;
	     Isum += data[i];
	     }
      }
      mu=Isum/cnz;

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

