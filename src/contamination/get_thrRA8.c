/*  File: get_rams.c // filters 1 1 vcfq file and computes av std Depth  depending on Depth of vcfq file

 * Authors: designed by Irina Abnizova (ia1) and edited by Steve Leonard(srl)
 *
  Last edited:
  17 April 2015-now 8 fields in vcfq file: '%POS,%REF,%ALT,%INFO/DP,%INFO/DP4\n' $vcf > $vcfq

  11 June: one thr file. output is e.g.  2 2 41 113
  9 June 2014 :n1 from pipe of three:filters 1 1 vcfq file and computes av std Depth  depending on Depth of vcfq file


  28 May : not output filtered vcfq
  22May  get_filters 1 1 .vcfq  RA_ms.filters

  21 May  filters bad varians: abnormal Depth
  30 April: adds Filteref 1,1 vcfq + distr_1_1
  29 April
  // computes suitable thresholds bases on avDP4=actual (not cov before as earlier!!!
  9 April: less strong threshold for AAF

  4 March, to make less strict thr for deep avDepth

  3 feb-put thrR thrA into file thrD ;29 Jan 2014  - get autom thr, starting with default ones 1,5
  still old fields of vcf input
  21 Nov 2013, for two thresholds

  *-------------------------------------------------------------------
 * Description: computes computes threshold file depending on Depth of vcfq file

 * Exported functions:
 * HISTORY:

computes minimum depth threshold, both for Ref and Alternative alleles fwd rev (total min depth will be
4*MIN_DEPTH


input: MIN_DEPTH_R(integer)=1[default]  MIN_DEPTH_A(integer)=1[default]   extracted=vcfq,

output1:
thrR  thrA
3      1
output2:  .distAF_1_1
output3:  .vcfqF_1_1  filtered 1,1
usage:./get_thrRA MIN_DEPTH_R(integer) MIN_DEPTH_A(integer) name.vcf  ra_i_j.thrD  name.vcfqF_1_1

*/


#include <stdio.h>
#include <time.h>
#include <math.h>
#include <string.h>

// ******** predefined user constants,

#define NRID 8      // Nunber of variant ids=fields: pos,Ref, Alt, Depth,DP4,(no PV4+maf) in extracted_vcf files
#define NPAR 5      // Number of PARameters and arguments


// *********    Global variables min_depths for Ref and Alternative alleles
int thR = 0;
int thA = 0;

// ******** declarations of  functions

int main    (int argc, char *argcv[]) {

    // flags
    int firstEightAreOK  = 1;
    //int canWriteDistrib      = 1;
    //int canWriteFilteredVCF      = 1;

    FILE *extract_vcfFile, *thrFile; //file handles

    int n,i;//
    int count_afterF=0;// after filtering with min depths
    int count_beforeF=0;// before filtering with min depths
    int DP4[4];//
    int D,pos;// Depth, genome pos of error
    int ref,alt;//nucleotides from vcf file

    //stats
    int sumDP4=0,sumDP4_2=0;//across all pos for mu std
    int sumDbefore=0;
    float avDbefore;
	float Q30frac;
	float sDP4;
	float perc_left;

    //to compute
    int avDP4,stdDP4,avDP41;// of actual sDP4 depth after 1 1 filtering
    int thRR, thAA;// recommended based on actual depth avDP4 and stdDP4
    float ratio;

      if(argc < NPAR)//four input_output files are submitted
    {
        printf("not enough of parms |input_output files\n");
        printf("usage:./get_rams thrR thrA input2 output1\n");
        return -1;
    }

    if (sscanf(argcv[1],"%d",&thR) == EOF) {
                printf("Failed to convert the first argument %s to integer\n", argcv[1]);
            }
            //printf(" min depth for Ref counts, fwd and rev =%d\n",thR);

     if (sscanf(argcv[2],"%d",&thA) == EOF) {
                printf("Failed to convert the second argument %s to integer\n", argcv[2]);
            }
           // printf(" min depth for Alternative allele counts, fwd and rev =%d\n",thA);

    extract_vcfFile=fopen(argcv[3],"r");
    if (extract_vcfFile == NULL) {
      printf("cannot open first input _vcf.txt file %s\n", argcv[3]);
      return -1;
    }
 // three outputs
     thrFile=fopen(argcv[4],"w");
            if (thrFile == NULL) {
            printf("cannot open third output file thrD %s\n", argcv[4]);
            return -1;
    }


        printf("get_rams\n");

    //6 fields of input file,CSV reading
    while( (n = fscanf(extract_vcfFile,"%d,%d,%d,%d,%d,%d,%d,%d", &pos, &ref, &alt,&D, &DP4[0], &DP4[1], &DP4[2], &DP4[3])) >= 0 && firstEightAreOK == 1)// && canWriteAF == 1)
    // read the Read Title
    {
        if( n != NRID )     // incorrect format
        {
            firstEightAreOK = 0;
            break;
        }

           count_beforeF++;
           sumDbefore=sumDbefore+D;



           // f1 Josie : filter for DP4 separately for ref and alternative alleles

           if ( DP4[0] >= thR && DP4[1]>= thR && DP4[2]>= thA && DP4[3]>= thA)// 0 ref are ok
           {
			     // for EACH POSITION, compute actual depths sDP4 after thrRA filtering:

                                // average run depth after Q25, sum(DP4) across pos,=  sumDP4
                                        count_afterF++;
                                        sDP4=DP4[0]+DP4[1]+DP4[2]+DP4[3];
                                        sumDP4=sumDP4+sDP4;//for av
                                        sumDP4_2=sumDP4_2+sDP4*sDP4;//for std

                                        Q30frac=sDP4/D;// how actaul quality Depth differs from default

            }// end filter DP4

        // removing space symbols before carriage return
        do  {
            n = fgetc (extract_vcfFile);
        }
        while ((char)n != '\n');

    } //END of outer while loop for all vcfq file



    // compute average depth : default and actual after 1,1 and q25 filtering

      // stats after 1,1 filtering

       avDP4=ceil(sumDP4/count_afterF);//actual average cov  after Filt
       stdDP4=ceil(sqrt(sumDP4_2/count_afterF));

       avDbefore=sumDbefore/count_beforeF;
       perc_left=(100.0*count_afterF/count_beforeF);


      // printf("average default Q13 Depth before min_depth filtering= %.2f\n", avDbefore);
       //printf("average actual Depth (avDP4) after Q25-filter and min_depth filter= %d\n", avDP4);
      // printf("std actual Depth (avDP4) after Q25-filter and min_depth filter= %d\n", stdDP4);
      // printf("percent data left after q25 & min_depth filtering 1 1 = %.2f\n", perc_left);

 // adjust actual depth to Bad regions and recommend thr RA

        avDP41=avDP4;
        ratio=(sqrt(sumDP4_2/count_afterF))/(sumDP4/count_afterF);//stdDP4/avDP4;
        printf("cv actual depth = %.2f\n", ratio);
        if ( ratio > 1.5) //bad regions
           {avDP41=ceil(0.8*avDP4);
           printf("moderate bad regions here\n");
           printf("adjusted(avDP4) = %d\n", avDP41);
           }
         if (ratio > 3) //very bad regions
		            {avDP41=ceil(0.7*avDP4);
		            printf("extreme bad regions here\n");
		            printf("adjusted(avDP4) = %d\n", avDP41);
           }
      // printf("adjusted(avDP4) = %d\n", avDP41);
           // recommendation for thrRR AA
           if ( avDP41 > 0 && avDP41 < 4 )
		     	     { thRR = 100;thAA=100;//stupidly high?
		     	     printf("better not to bother to look for mixture: too low actual coverage=%d\n",avDP4);
		     	     }
		     	     if ( avDP41 >= 4 && avDP41 < 10  )
		     		      { thRR = 1;thAA=1;
		     	     }
		     	     if ( avDP41 >= 10 && avDP41 < 70  )
		     		 { thRR = 2;thAA=2;
		     	     }
		     	     if ( avDP41 >= 70 && avDP41 < 90  )
		     		 { thRR = 3;thAA=2;
		     	     }
		     	     if ( avDP41 >= 90 )
		     		 { thRR = 3;thAA=3;
		  	     }

 // fill in the output file
        fprintf(thrFile,"%d %d %d %d\n",thRR,thAA,avDP4,stdDP4);


    fclose(extract_vcfFile);
    fclose(thrFile);

       // checking write/read
        if( firstEightAreOK  == 0)// || canWriteFilteredVCF == 0)
        {
		        printf ("Error during execution. Details: \n");
		        printf ("\tfirstEightAreOK %d\n",  firstEightAreOK);
		        //printf ("\tcanWriteDistrib %d\n",      canWriteDistrib);
		         //printf ("\tcanWriteFilteredVCF %d\n",   canWriteFilteredVCF);
		        printf ("Execution aborted\n");
		        return -1;
        }

    printf(" done get_rams.\n");
    return 0;
}//main
