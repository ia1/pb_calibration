
//ia1@sf2-farm-srv2:/nfs/sf39/tmp_ia1$
//ia1@sf2-farm-srv2:/nfs/sf39/tmp_ia1$ more get_thr_min.c//

 /*  File: get_thr_min.c // computes maf distribution over variant positions from vcf files
 * Authors: designed by Irina Abnizova (ia1) and edited by Steve Leonard(srl)
 *
  Last edited:
  4 March, to make less strict thr for deep avDepth
  22 feb-min of peaks outputs
  3 feb-put thrR thrA into file thrD ;29 Jan 2014  - get autom thr, starting with default ones 1,5
  still old fields of vcf input
  21 Nov 2013, for two thresholds

  *-------------------------------------------------------------------
 * Description: computes maf distrib over variant positions from just extracted 6 fileld

 * Exported functions:
 * HISTORY:

computes AF (alternative to Ref allele frequency) for each error base call from vcf

applies minimum depth threshold, both for Ref and Alternative alleles fwd rev (total min depth will be
4*MIN_DEPTH


input: extracted vcf,


output3:
thrR  thrA
0      1

usage:./get_thr_min MIN_DEPTH_R(integer) MIN_DEPTH_A(integer) extracted.vcf  output.thrD

*/


#include <stdio.h>
#include <time.h>
#include <math.h>
#include <string.h>

// ******** predefined user constants,

#define NRID 6      // Number of variant ids=fields: pos, chromosome, Depth,DP4,(no PV4+maf) in extracted_vcf files
#define NPAR 5      // max Number of PARameters and arguments+thrD file


// *********    Global variables min_depths for Ref and Alternative alleles
int thR = 0;
int thA = 0;

// ******** declarations of  functions
//int   GetAf (int []);// computes AF percent for a variant position

int main    (int argc, char *argcv[]) {

    // flags
    int firstSixAreOK  = 1;
   // int lastSevenAreRight = 1;
    //int canWriteAF        = 1;
    //int canWriteDistrib      = 1;

    FILE *extract_vcfFile, *thrFile; //file handles

    int n,i,perc;//, cntZero=0, counter=0;
    //char chrom[50];
    int count_afterF=0;// after filtering with min depths
    int count_beforeF=0;// before filtering with min depths
    int DP4[4];//
    int D,pos;// Depth, genome pos of error
    int sumD=0,sumDP4=0;//across all pos!
    int sumDbefore=0;
    int thRR,thAA;// to compute!
    //int thr_R_A[2];
    //int AF;

    //float mu,std;//mean and std of alternative frequency, AF

    float avDbefore;
    float avD,avDP4,Q30frac;
    float sDP4;
    float perc_left;

      if(argc < NPAR)//four input_output files are submitted
    {
        printf("not enough of parms |input_output files\n");
        printf("usage:./get_thr_min thrR thrA input2 output1\n");
        return -1;
    }

    if (sscanf(argcv[1],"%d",&thR) == EOF) {
                printf("Failed to convert the first argument %s to integer\n", argcv[1]);
            }
            printf(" min depth for Ref counts, fwd and rev =%d\n",thR);

     if (sscanf(argcv[2],"%d",&thA) == EOF) {
                printf("Failed to convert the second argument %s to integer\n", argcv[2]);
            }
            printf(" min depth for Alternative allele counts, fwd and rev =%d\n",thA);

    extract_vcfFile=fopen(argcv[3],"r");
    if (extract_vcfFile == NULL) {
      printf("cannot open first input _vcf.txt file %s\n", argcv[3]);
      return -1;
    }

     thrFile=fopen(argcv[4],"w");
            if (thrFile == NULL) {
            printf("cannot open third output file thrD %s for writing\n", argcv[4]);
            return -1;
    }

    printf("get_thr_min\n");

    //6 fields of input file,CSV reading
    while( (n = fscanf(extract_vcfFile,"%d,%d,%d,%d,%d,%d", &pos, &D, &DP4[0], &DP4[1], &DP4[2], &DP4[3])) >= 0 && firstSixAreOK == 1)// && canWriteAF == 1)
    // read the Read Title
    {
        if( n != NRID )     // incorrect format
        {
            firstSixAreOK = 0;
            break;
        }

           count_beforeF++;
           sumDbefore=sumDbefore+D;

           // f1 Josie : filter for DP4 separately for ref and alternative alleles

           if ( DP4[0] >= thR && DP4[1]>= thR && DP4[2]> thA && DP4[3]> thA)// 0 ref are ok
           {

                        // compute average run depth default (Q13), D across pos,sumD
                                // compute average run depth after Q30, sum(DP4) across pos,sumDP4
                                        count_afterF++;
                                        sumD=sumD+D;
                                        sDP4=DP4[0]+DP4[1]+DP4[2]+DP4[3];
                                        sumDP4=sumDP4+sDP4;
                                        Q30frac=sDP4/D;
                 //printf("frac filtered out Q= %.2f\n", Q30frac);

            }// end filter DP4

        // removing space symbols before carriage return
        do  {
            n = fgetc (extract_vcfFile);
        }
        while ((char)n != '\n');

    } //END of outer while loop

    // compute avarege depth :default and after filtering
       avD=sumD/count_afterF;
       avDP4=sumDP4/count_afterF;
       avDbefore=sumDbefore/count_beforeF;
       perc_left=(100.0*count_afterF/count_beforeF);

      // printf("average default Depth before filtering= %.2f\n", avDbefore);
       //printf("average Depth default after min_depth filter= %.2f\n", avD);
       printf("average actual Depth after Q25-filter and min_depth filter= %.2f\n", avDP4);

       //printf("percent data left after min_depth filtering= %.2f\n", perc_left);

    // computes suitable thresholds bases on avD

      for (i=0; i<100; i++)
           {
                         if ( avDbefore > i*10 && avDbefore <= (i+1)*10 )
                         {
                         //thAA= i;
                         thRR = (long)( avDbefore*0.07 +0.5);// or i+2 ???
                         thAA= thRR-3;
                         }
           }
           if (thRR<3){
			   thAA=0;
		   }

        printf("recommended Ref thresholds= %d\n", thRR);
        printf("recommended Alternative threshold= %d\n", thAA);

        fprintf(thrFile,"%d %d\n",thRR,thAA);


    fclose(extract_vcfFile);

    fclose(thrFile);

       if( firstSixAreOK  == 0 )//||
        //canWriteAF == 0 || canWriteDistrib ==0)
        {
        printf ("Error during execution. Details: \n");
        printf ("\tfirstSixAreOK %d\n",  firstSixAreOK);
        printf ("Execution aborted\n");
        return -1;
    }

    printf(" done get_thr_min.\n");
    return 0;
}//main

