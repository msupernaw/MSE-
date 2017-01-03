/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/*
 * File:   MSE.hpp
 * Author: matthewsupernaw
 *
 * Created on January 19, 2016, 10:56 AM
 * Updated on December 27, 2016 by ZTA
 *
 */


/* Box Muller Adapted from GSL Libraries  - A.W. Seaver*/

/* GNU General Public License http://www.gnu.org/licenses/gpl.html
 * randist/gauss.c
 *
 * Copyright (C) 1996, 1997, 1998, 1999, 2000, 2006, 2007 James Theiler, Brian Gough
 * Copyright (C) 2006 Charles Karney
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 */


/* Box Muller Adapted from GSL Libraries  - A.W. Seaver*/


/* Polar (Box-Muller) method; See Knuth v2, 3rd ed, p122 */


#ifndef MSE_HPP
#define MSE_HPP

#include "util.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <random>
#include <iostream>
#include <memory>

#define   MAXBUF    2048
#define   FILBUF     256
#define   TXTBUF     128

#define   MXITER     100
#define   NRECMAX    500
#define   MCLEN       25

#define   XTOL   1.0E-15
#define   XXLW         6
#define   XXNM         2
#define   XXDF         2

#define   XXSD     10000

#define   XLBOUND  0.0001
#define   XUBOUND  5.0

#define   TAGSIZE  16

#define   FORMFEED  '\x0C'

#define   NOPT        10
#define   NPHASE       8
#define   NSELMOD      6
#define   NSTRING     21

#define  AGEPROCV      9

#define   M_SQRT1_2 (1.0 / 1.41421356237309504880168872421)


namespace noaa {

    template<class T = double>
    class MSE {
        std::mt19937 generator; //// mt19937 is a standard mersenne_twister_engine
        std::normal_distribution<T> distribution;
        std::uniform_real_distribution<T> uniform_distribution;
        std::string input;
        std::string output;

        /* Global Variables */

        short assessment_frequency;     // MRS addition to v3.4.3
        double Blim;                    // ZTA addition

        FILE *fp1,*fp2, *fp3, *fp4, *fp5, *fp6;

        long NFYear;
        long NXYear;
        long NYears;

        long MXYear;
        long KYears;

        long MaxYears;

        long MaxAge;
        long NAgePlus;

        long NSurvey;

        long NIter;

        long ModelFlag;

        long ISeed;

        long DiscFlag;

        long RecFlag;

        long NRec1;
        long NRec2;

        double RecAlpha;
        double RecBeta;
        double RecGamma;
        double CvRec;
        double CvPop;

        double RecScale;
        double SSBScale;

        double *RecObs1;
        double *RecObs2;
        double *ObsRec;

        double SSBCutPt;


        char CaseID[TXTBUF];

        double *InitPop;
        double *TF;
        double *TM;
        double *Recruits;

        long *HarvestSpec;

        double *FullF;
        double *Landings;
        double *Discards;
        double *Biomass;
        double *CatchBio;
        double *SSB;
        double *HarvestValue;

        double **StockNumbers;
        double **StockBiomass;

        double **CatchNumbers;
        double **CatchBiomass;

        double **SpawnNumbers;
        double **SpawnBiomass;

        double **NatMort;
        double **FMort;
        double **Mature;

        double **StockWeight;
        double **SpawnWeight;
        double **CatchWeight;


        double **FSelec;
        double **DiscFrac;

        double *AgeComp;
        double *CatchAge;

        double *CatchCV;
        double *CatchDraw;

        double **CatchSamples;

        double *SurveyAge;

        double **SurveyQ;
        double **SurveyCV;
        double **SurveyDraw;

        double ***SurveySamples;
        double ***SurveySelec;


        double **StockEstim;
        double **SSBEstim;
        double **FEstim;

        long KIter;

        char WorkPath[FILBUF];

        /* MSR */

        double Fmsy;
        double Bmsy;

        double SSBEst;
        double FEst;
        double FTarg;
        double FPrev;
        double TargetLand;
        double ProjLand;
        double LandEst;

        long DeltaFlag;

        double DeltaLand;
        double PrevLand;
        double PrevProj;


        double biasSSB;
        double biasFEst;
        double biasLand;

        double cvSSB;
        double cvFEst;
        double cvLand;


        /* Template Files */

        char AgeProPath[FILBUF];

        char ModelPath[FILBUF];

        long NBoot;

        long NIndex;

        long KAges;

        struct si
        {
            short StartAge;
            short EndAge;
            short Type;
            short SurveyIndex;
        };

        struct si *IndexData;

        short Fold[2];

        double *PRVec;
        double *StockEst;

        double **NMortMat;
        double **IndexValues;

        double *AvgFVPA;

        double CatchFill;

        long *AgeEst;

        short OptionFlags[NOPT];

        short NRange[2];

        double *NMortVec;
        double **RSelec;


        double *StockVec;
        double *FMortVec;

        long NRecAge;

        double NMortEst;

        /* VPA */


        double **StockVPA;
        double **SSBVPA;
        double **FVPA;

        /* AgePro */

        long NAgeProSim;

        double cvAgeProErr;
        double biasAgeProErr;

        double **AgeProCV;


        /* ASAP */

        short NBlock;

        short *Fleetblock;
        short *FleetSelType;
        short *IndexSelType;

        short PhaseAsap[NPHASE];

        short AvgF[2];

        short NAAFlag;
        short SRRFlag;

        long KThin;
        long KStart;
        long KEnd;

        double *Fleetsamp;
        double *Fleetcv;
        double *Discsamp;
        double *Disccv;
        double *Recruitcv;

        double *LambdaIndex;
        double *Lambdaq;
        double *Lambdaqdev;


        double LambdaFleet;
        double LambdaDisc;
        double LambdaFMult;
        double LambdaFDev;
        double LambdaN1;
        double LambdaRec;
        double LambdaSteep;
        double LambdaVirgin;

        double cvFMult;
        double cvFDev;
        double cvN1;
        double cvSteep;
        double cvVirgin;

        double *cvq;
        double *cvqdev;

        double *InitNAA;
        double *Initq;

        double InitFMult;
        double InitVirgin;
        double InitSteep;
        double InitFMax;
        double IndexDatcv;

        double **StockAsap;

        double SSBAsap;


        struct seldat
        {
            double value;
            double lambda;
            double cv;
            short phase;
        };

        struct seldat *Fleetsel;
        struct seldat *Indexsel;

        /*  ASPIC */

        char AspicPath[FILBUF];
        char AspicpPath[FILBUF];


        double *AspicF;
        double *AspicB;
        double *AspicY;

        char *PString[NSTRING];

        double cvAspicErr;
        double biasAspicErr;


        void RNOPT(long* o) {

        }

        void RNSET(long* seed) {
            generator.seed(*seed);
        }

        void RNGET(long* seed) {
            *seed = this->ISeed;
        }

        static T Erf(const T &x) {
            // constants
            T a1 = T(0.254829592);
            T a2 = T(-0.284496736);
            T a3 = T(1.421413741);
            T a4 = T(-1.453152027);
            T a5 = T(1.061405429);
            T p = T(0.3275911);

            // Save the sign of x
            T sign = T(1);
            if (x < T(0))
                sign = T(-1);
            T xx = std::fabs(x);

            // A&S formula 7.1.26
            T t = T(1.0) / (T(1.0) + p * xx);
            T y = T(1.0) - T((((((a5 * t + a4) * t) + a3) * t + a2) * t + a1) * t * std::exp(T(-1.0) * xx * xx));

            return sign*y;
        }

        const T pnorm(const T &x, const T &m, const T &s) {
            //            return T(0.5)*(T(1.0) + Erf((x - mean) / std::sqrt(T(2.0) * sd * sd)));
            static const T inv_sqrt_2pi = 0.3989422804014327;
            T a = (x - m) / s;

            return inv_sqrt_2pi / s * std::exp(-T(0.5) * a * a);
        }

        T DNORDF(T* x) {
//            return pnorm(*x, this->distribution.mean(), this->distribution.stddev());
             return 0.5 * erfc(-*x * M_SQRT1_2);
        }

        T DRNUNF() {
            return uniform_distribution(this->generator);
        }

        void DRNNOR(long* n, T* array) {
            for (int i = 0; i < (*n); i++) {
                array[i] = this->distribution(this->generator);
            }
        }

        double gsl_box_muller()
        {
          double x, y, r2;

          do
            {
              /* choose x,y in uniform square (-1,-1) to (+1,+1) */
              x = -1.0 + 2.0 * DRNUNF();
              y = -1.0 + 2.0 * DRNUNF();

              /* see if it is in the unit circle */
              r2 = x * x + y * y;
            }
          while (r2 > 1.0 || r2 == 0);

          /* Box-Muller transform */
          return (y * sqrt (-2.0 * log (r2) / r2));
        }

        std::string exec(const char* cmd, int& error) {
            std::shared_ptr<FILE> pipe(_popen(cmd, "r"), _pclose);
            if (!pipe) {
                error = 0;
                return "";
            }
            char buffer[128];
            std::string result = "";
            while (!feof(pipe.get())) {
                if (fgets(buffer, 128, pipe.get()) != NULL) {
                    std::cout << buffer << std::flush;
                }
                //                    result += buffer;
            }
            return result;
        }

        long CreateConsoleProcess(char *cmd) {
            std::cout << cmd << "\n";
            int error = 1;
            std::string out = this->exec(cmd, error);
            std::cout << out.c_str() << std::endl;
            return error;
        }

    public:

        MSE(int argc, char **argv) {
            this->assessment_frequency = 1;

            long iopt = 5;

            char *fn;
            char *c;
            char fname[FILBUF];
            char dstrng[15];
            char tstrng[15];
            long i;
            long iter;
            long kx;
            short flag;
            double tpop;
            double sd, zscore, xx;
            double ssbx;


            if (argc < 2)
            {
                fprintf(stderr,"usage: mse41 filename\n");
                exit(1);
            }

            fn = *++argv;

            strcpy(fname,fn);

            DateString(dstrng);
            TimeString(tstrng);

            strcpy(WorkPath,fname);

            if ((c = strrchr(WorkPath,'\\')) != NULL)
                *c = '\0';

            ReadInputDataFile(fname);

            /* Open Report File */

            c = strrchr(fname,'.');

            *c = '\0';

            strcat(fname,".rep");

            if ((fp2 = fopen(fname,"w")) == NULL)
            {
                fprintf(stderr,"Unable to Open Report File: %s\n",fname);
                exit(1);
            }

            /* Open Auxiliary Files */

            *c = '\0';

            strcat(fname,"_xx1.csv");

            if ((fp3 = fopen(fname,"w")) == NULL)
            {
                fprintf(stderr,"Unable to Open Auxiliary File: %s\n",fname);
                exit(1);
            }


            *c = '\0';

            strcat(fname,"_xx2.csv");

            if ((fp4 = fopen(fname,"w")) == NULL)
            {
                fprintf(stderr,"Unable to Open Auxiliary File: %s\n",fname);
                exit(1);
            }


            fprintf(fp2,"Management Strategy Evaluation Version 4.1\n\n");
            fprintf(fp2,"Date & Time of Run: %s  %s\n\n",dstrng,tstrng);
            fprintf(fp2,"Case ID: %s\n\n",CaseID);

            if (ModelFlag == 0)
            {
                fprintf(fp2,"General Projection Model\n\n");

                fprintf(fp4,"Year,");
                fprintf(fp4,"F-Full,");
                fprintf(fp4,"Removals,");
                fprintf(fp4,"SSB (N-1),");
                fprintf(fp4,"SSB Estim,");
                fprintf(fp4,"F-Target,");
                fprintf(fp4,"F-Estimate (N+1),\n");
            }
            else if (ModelFlag == 1)
            {
                fprintf(fp2,"Model - VPA / AgePro\n\n");

                ReadVPATemplateFile(fn);
                ReadAgeProTemplateFile(fn);

                fprintf(fp4,"Year,");
                fprintf(fp4,"F-Full,");
                fprintf(fp4,"Removals,");
                fprintf(fp4,"SSB (VPA),");
                fprintf(fp4,"Target-F,");
                fprintf(fp4,"AgePro Catch,");
                fprintf(fp4,"Target Catch,");
                fprintf(fp4,"Proj Catch,");
                fprintf(fp4,"F-Calc (N+1),\n");
            }
            else if (ModelFlag == 2)
            {
                fprintf(fp2,"Model - ASAP / AgePro\n\n");

                ReadAsapTemplateFile(fn);
                ReadAgeProTemplateFile(fn);

                fprintf(fp4,"Year,");
                fprintf(fp4,"F-Full,");
                fprintf(fp4,"Removals,");
                fprintf(fp4,"SSB (VPA),");
                fprintf(fp4,"Target-F,");
                fprintf(fp4,"AgePro Catch,");
                fprintf(fp4,"Target Catch,");
                fprintf(fp4,"Proj Catch,");
                fprintf(fp4,"F-Calc (N+1),\n");
            }
            else if (ModelFlag == 3)
            {
                fprintf(fp2,"Model - Aspic / Aspic Projection\n\n");

                ReadAspicTemplateFile(fn);

                fprintf(fp4,"Year,");
                fprintf(fp4,"F-Full,");
                fprintf(fp4,"Removals,");
                fprintf(fp4,"SSB (Aspic),");
                fprintf(fp4,"Target-F,");
                fprintf(fp4,"Aspic-P Catch,");
                fprintf(fp4,"Target Catch,");
                fprintf(fp4,"Proj Catch,");
                fprintf(fp4,"F-Calc (N+1),\n");
            }

            /* Expand Data to MaxYears + 1 */

            ExpandData();


            this->assessment_frequency = this->assessment_frequency < 1 ? 1 : this->assessment_frequency;
            std::cout << "assessment frequency = " << this->assessment_frequency << "\n";
            std::cout << "done.\n";


            /* Initialize Random Number Generator */

            RNOPT(&iopt);

            RNSET(&ISeed);

            /* Normalize Initial Composition at Age */

            tpop = 0.0;
            for (i = 0; i < MaxAge; i++)
                tpop += InitPop[i];

            for (i = 0; i < MaxAge; i++)
                InitPop[i] = InitPop[i] / tpop;

            for (iter = 0; iter < NIter; iter++)
            {
                printf("Iteration - %d\n\n",iter+1);
                fprintf(fp2,"Iteration - %d\n\n",iter+1);
                fprintf(fp3,"Iteration, %ld,\n",iter+1);
                fprintf(fp4,"Iteration, %ld,\n",iter+1);

                /* Apply Error to Initial Population (Total) */

                sd = sqrt(log(CvPop * CvPop + 1.0));

                zscore = gsl_box_muller();

                xx = tpop * exp(sd * zscore - sd * sd * 0.5);

                ZeroMatrix(StockNumbers,MaxAge,MaxYears+1);

                for (i = 0; i < MaxAge; i++)
                    StockNumbers[i][0] = xx * InitPop[i];

                fprintf(fp2,"Initial Population = %18.0f\n\n",xx);

                Recruits[0] = StockNumbers[0][0];

                /* Generate Calculated Population & Catch over Base Years */

                for (i = 0; i < NYears; i++)
                {
                    /* Calculate Jan-1 Biomass in Current Year */

                    CalcStockBiomass(i);

                    /* Apply Mortality in Current Year */

                    CalcCatch(i);

                    /* Insert Recruitment */

                    ApplyRecruitment(i);
                }

                CalcStockBiomass(NYears);

                if (ModelFlag == 0)
                {
                    fprintf(fp2,"Year     ");
                    fprintf(fp2,"F-Full      ");
                    fprintf(fp2,"Removals     ");
                    fprintf(fp2,"SSB (N-1)    ");
                    fprintf(fp2,"SSB Estim    ");
                    fprintf(fp2,"F-Target   ");
                    fprintf(fp2,"F-Estim (N+1)\n\n");

                    FPrev = FullF[NYears-1];
                }
                else if (ModelFlag == 1)
                {
                    fprintf(fp2,"Year     ");
                    fprintf(fp2,"F-Full      ");
                    fprintf(fp2,"Removals     ");
                    fprintf(fp2,"SSB (VPA)    ");
                    fprintf(fp2,"Target-F    ");
                    fprintf(fp2,"AgePro Catch ");
                    fprintf(fp2,"Target Catch ");
                    fprintf(fp2,"Proj Catch   ");
                    fprintf(fp2,"F-Calc (N+1)\n\n");
                }
                else if (ModelFlag == 2)
                {
                    fprintf(fp2,"Year     ");
                    fprintf(fp2,"F-Full      ");
                    fprintf(fp2,"Removals     ");
                    fprintf(fp2,"SSB (ASAP)   ");
                    fprintf(fp2,"Target-F    ");
                    fprintf(fp2,"AgePro Catch ");
                    fprintf(fp2,"Target Catch ");
                    fprintf(fp2,"Proj Catch   ");
                    fprintf(fp2,"F-Calc (N+1)\n\n");
                }
                else if (ModelFlag == 3)
                {
                    fprintf(fp2,"Year     ");
                    fprintf(fp2,"F-Full      ");
                    fprintf(fp2,"Removals     ");
                    fprintf(fp2,"SSB (Aspic)   ");
                    fprintf(fp2,"Target-F    ");
                    fprintf(fp2,"Aspic-P Catch ");
                    fprintf(fp2,"Target Catch ");
                    fprintf(fp2,"Proj Catch   ");
                    fprintf(fp2,"F-Calc (N+1)\n\n");
                }


                while (NYears <= MaxYears)
                {
                    if (this->assessment_frequency > 1 && (NYears % this->assessment_frequency) != 0)
                    {
                        NYears++;

                        continue;
                    }

                    if (ModelFlag == 0)
                    {

                        TargetLand = Landings[NYears-1] + Discards[NYears-1];

                        ssbx = SSB[NYears-2];

                        sd = sqrt(log(cvSSB * cvSSB + 1.0));
                        zscore = gsl_box_muller();
                        xx = exp(sd * zscore);
                        SSBEst = ssbx * (1.0 + biasSSB) * xx;

                        FTarg  = MSRule(SSBEst);

                        zscore = gsl_box_muller();
                        xx = exp(sd * zscore);
                        FEst = FTarg * (1.0 + biasFEst) * xx;

                        flag = 0;

                        if (DeltaFlag && NYears > KYears)
                        {
                            xx = (FEst - FPrev) / FPrev;
                            if (xx > DeltaLand)
                            {
                                FEst = FPrev * (1.0 + DeltaLand);
                                flag = 1;
                            }
                            else if (fabs(xx) > DeltaLand)
                            {
                                FEst = FPrev * (1.0 - DeltaLand);
                                flag = 1;
                            }
                        }


                        FullF[NYears] = FEst;

                        HarvestSpec[NYears] = 0;

                        CalcStockBiomass(NYears);

                        CalcCatch(NYears);

                        ApplyRecruitment(NYears);

                        fprintf(fp2,"%ld  ",NFYear+NYears-1);
                        fprintf(fp2,"%10.5f ",FullF[NYears-1]);
                        fprintf(fp2,"%12.3f ",TargetLand);
                        fprintf(fp2,"%12.3f ",ssbx);
                        fprintf(fp2,"%12.3f ",SSBEst);
                        fprintf(fp2,"%10.5f ",FTarg);
                        fprintf(fp2,"%10.5f ",FEst);


                        fprintf(fp4,"%ld,",NFYear+NYears-1);
                        fprintf(fp4,"%10.5f,",FullF[NYears-1]);
                        fprintf(fp4,"%12.3f,",TargetLand);
                        fprintf(fp4,"%12.3f,",ssbx);
                        fprintf(fp4,"%12.3f,",SSBEst);
                        fprintf(fp4,"%10.5f,",FTarg);
                        fprintf(fp4,"%10.5f,",FEst);

                        if (flag)
                        {
                            fprintf(fp2," *\n");
                            fprintf(fp4,"*,\n");
                        }
                        else
                        {
                            fprintf(fp2,"\n");
                            fprintf(fp4,"\n");
                        }


                        FPrev = FEst;

                        printf("Last Year in Projection = %ld\n\n",NFYear+NYears-1);
                    }
                    else if (ModelFlag == 1)
                    {
                        TargetLand = Landings[NYears-1] + Discards[NYears-1];

                        strcpy(fname,fn);

                        DrawCatchSamples();
                        DrawSurveySamples();
                        WriteVPAInputFile(fname);
                        RemoveOutputFiles(fname);

                        kx = LaunchVPA(fname);
                        if (kx == 0)
                        {
                            fprintf(stderr,"Error Launching VPA !\n");
                            exit(1);
                        }

                        ScanVPAResults(fname);

                        SSBEst = 0.0;
                        for (i = 0; i < NAgePlus; i++)
                            SSBEst += SSBVPA[i][NYears-1] / 1000.0;

                        FTarg  = MSRule(SSBEst);


                        WriteAgeProInputFile(fname);
                        kx = LaunchAgePro(fname);
                        if (kx == 0)
                        {
                            fprintf(stderr,"Error Launching AgePro !\n");
                            exit(1);
                        }
                        ScanAgeProResults(fname);

                        flag = 0;

                        if (NYears == KYears)
                        {
                            LandEst = TargetLand;
                            PrevProj = TargetLand;
                        }
                        else
                        {
                            LandEst = PrevProj;

                            if (DeltaFlag)
                            {
                                xx = (PrevProj - PrevLand) / PrevLand;
                                if (xx > DeltaLand)
                                {
                                    LandEst = PrevLand * (1.0 + DeltaLand);
                                    flag = 1;
                                }
                                else if (fabs(xx) > DeltaLand)
                                {
                                    LandEst = PrevLand * (1.0 - DeltaLand);
                                    flag = 1;
                                }
                            }
                        }

                        sd = sqrt(log(cvLand * cvLand + 1.0));
                        zscore = gsl_box_muller();
                        xx = exp(sd * zscore);
                        LandEst = LandEst * (1.0 + biasLand) * xx;

                        if (DiscFlag)
                            HarvestSpec[NYears] = 2;
                        else
                            HarvestSpec[NYears] = 1;

                        HarvestValue[NYears] = LandEst;

                        CalcStockBiomass(NYears);

                        CalcCatch(NYears);

                        ApplyRecruitment(NYears);

                        fprintf(fp2,"%ld  ",NFYear+NYears-1);
                        fprintf(fp2,"%10.5f ",FullF[NYears-1]);
                        fprintf(fp2,"%12.3f ",TargetLand);
                        fprintf(fp2,"%12.3f ",SSBEst);
                        fprintf(fp2,"%10.5f ",FTarg);
                        fprintf(fp2,"%12.3f ",ProjLand);
                        fprintf(fp2,"%12.3f ",PrevProj);
                        fprintf(fp2,"%12.3f ",LandEst);
                        fprintf(fp2,"%10.5f ",FullF[NYears]);


                        fprintf(fp4,"%ld,",NFYear+NYears-1);
                        fprintf(fp4,"%10.5f,",FullF[NYears-1]);
                        fprintf(fp4,"%12.3f,",TargetLand);
                        fprintf(fp4,"%12.3f,",SSBEst);
                        fprintf(fp4,"%10.5f,",FTarg);
                        fprintf(fp4,"%12.3f,",ProjLand);
                        fprintf(fp4,"%12.3f,",PrevProj);
                        fprintf(fp4,"%12.3f,",LandEst);
                        fprintf(fp4,"%10.5f,",FullF[NYears]);

                        if (flag)
                        {
                            fprintf(fp2," *\n");
                            fprintf(fp4,"*,\n");
                        }
                        else
                        {
                            fprintf(fp2,"\n");
                            fprintf(fp4,"\n");
                        }


                        PrevLand = LandEst;

                        PrevProj = ProjLand;

                    }
                    else if (ModelFlag == 2)
                    {
                        TargetLand = Landings[NYears-1] + Discards[NYears-1];

                        strcpy(fname,fn);

                        DrawCatchSamples();
                        DrawSurveySamples();
                        RemoveOutputFiles(fname);
                        WriteAsapInputFile(fname);

                        kx = LaunchAsap(fname);
                        if (kx == 0)
                        {
                            fprintf(stderr,"Error Launching ASAP !\n");
                            exit(1);
                        }
                        ScanAsapResults(fname);

                        SSBEst = SSBAsap;

                        FTarg  = MSRule(SSBEst);


                        WriteAgeProInputFile(fname);
                        kx = LaunchAgePro(fname);
                        if (kx == 0)
                        {
                            fprintf(stderr,"Error Launching AgePro !\n");
                            exit(1);
                        }
                        ScanAgeProResults(fname);

                        flag = 0;

                        if (NYears == KYears)
                        {
                            LandEst = TargetLand;
                            PrevProj = TargetLand;
                        }
                        else
                        {
                            LandEst = PrevProj;

                            if (DeltaFlag)
                            {
                                xx = (PrevProj - PrevLand) / PrevLand;
                                if (xx > DeltaLand)
                                {
                                    LandEst = PrevLand * (1.0 + DeltaLand);
                                    flag = 1;
                                }
                                else if (fabs(xx) > DeltaLand)
                                {
                                    LandEst = PrevLand * (1.0 - DeltaLand);
                                    flag = 1;
                                }
                            }
                        }

                        sd = sqrt(log(cvLand * cvLand + 1.0));
                        zscore = gsl_box_muller();
                        xx = exp(sd * zscore);
                        LandEst = LandEst * (1.0 + biasLand) * xx;

                        if (DiscFlag)
                            HarvestSpec[NYears] = 2;
                        else
                            HarvestSpec[NYears] = 1;

                        HarvestValue[NYears] = LandEst;

                        CalcStockBiomass(NYears);

                        CalcCatch(NYears);

                        ApplyRecruitment(NYears);

                        fprintf(fp2,"%ld  ",NFYear+NYears-1);
                        fprintf(fp2,"%10.5f ",FullF[NYears-1]);
                        fprintf(fp2,"%12.3f ",TargetLand);
                        fprintf(fp2,"%12.3f ",SSBEst);
                        fprintf(fp2,"%10.5f ",FTarg);
                        fprintf(fp2,"%12.3f ",ProjLand);
                        fprintf(fp2,"%12.3f ",PrevProj);
                        fprintf(fp2,"%12.3f ",LandEst);
                        fprintf(fp2,"%10.5f ",FullF[NYears]);


                        fprintf(fp4,"%ld,",NFYear+NYears-1);
                        fprintf(fp4,"%10.5f,",FullF[NYears-1]);
                        fprintf(fp4,"%12.3f,",TargetLand);
                        fprintf(fp4,"%12.3f,",SSBEst);
                        fprintf(fp4,"%10.5f,",FTarg);
                        fprintf(fp4,"%12.3f,",ProjLand);
                        fprintf(fp4,"%12.3f,",PrevProj);
                        fprintf(fp4,"%12.3f,",LandEst);
                        fprintf(fp4,"%10.5f,",FullF[NYears]);

                        if (flag)
                        {
                            fprintf(fp2," *\n");
                            fprintf(fp4,"*,\n");
                        }
                        else
                        {
                            fprintf(fp2,"\n");
                            fprintf(fp4,"\n");
                        }


                        PrevLand = LandEst;

                        PrevProj = ProjLand;

                    }
                    else if (ModelFlag == 3)
                    {
                        TargetLand = Landings[NYears-1] + Discards[NYears-1];

                        strcpy(fname,fn);

                        DrawCatchSamples();
                        DrawSurveySamples();
                        RemoveOutputFiles(fname);
                        WriteAspicInputFile(fname);

                        kx = LaunchAspic(fname);
                        if (kx == 0)
                        {
                            fprintf(stderr,"Error Launching Aspic\n");
                            exit(1);
                        }

                        CaptureAspicResults(fname);

                        SSBEst = AspicB[NYears-1] / 1000.0;

                        FTarg  = MSRule(SSBEst);

                        WriteAspicPFile(fname);

                        kx = LaunchAspicProj(fname);

                        if (kx == 0)
                        {
                            fprintf(stderr,"Error Launching Aspic Projection\n");
                            exit(1);
                        }

                        ScanAspicProj(fname);

                        flag = 0;

                        if (NYears == KYears)
                        {
                            LandEst = TargetLand;
                            PrevProj = TargetLand;
                        }
                        else
                        {
                            LandEst = PrevProj;

                            if (DeltaFlag)
                            {
                                xx = (PrevProj - PrevLand) / PrevLand;
                                if (xx > DeltaLand)
                                {
                                    LandEst = PrevLand * (1.0 + DeltaLand);
                                    flag = 1;
                                }
                                else if (fabs(xx) > DeltaLand)
                                {
                                    LandEst = PrevLand * (1.0 - DeltaLand);
                                    flag = 1;
                                }
                            }
                        }

                        sd = sqrt(log(cvLand * cvLand + 1.0));
                        zscore = gsl_box_muller();
                        xx = exp(sd * zscore);
                        LandEst = LandEst * (1.0 + biasLand) * xx;

                        if (DiscFlag)
                            HarvestSpec[NYears] = 2;
                        else
                            HarvestSpec[NYears] = 1;

                        HarvestValue[NYears] = LandEst;

                        CalcStockBiomass(NYears);

                        CalcCatch(NYears);

                        ApplyRecruitment(NYears);

                        fprintf(fp2,"%ld  ",NFYear+NYears-1);
                        fprintf(fp2,"%10.5f ",FullF[NYears-1]);
                        fprintf(fp2,"%12.3f ",TargetLand);
                        fprintf(fp2,"%12.3f ",SSBEst);
                        fprintf(fp2,"%10.5f ",FTarg);
                        fprintf(fp2,"%12.3f ",ProjLand);
                        fprintf(fp2,"%12.3f ",PrevProj);
                        fprintf(fp2,"%12.3f ",LandEst);
                        fprintf(fp2,"%10.5f ",FullF[NYears]);


                        fprintf(fp4,"%ld,",NFYear+NYears-1);
                        fprintf(fp4,"%10.5f,",FullF[NYears-1]);
                        fprintf(fp4,"%12.3f,",TargetLand);
                        fprintf(fp4,"%12.3f,",SSBEst);
                        fprintf(fp4,"%10.5f,",FTarg);
                        fprintf(fp4,"%12.3f,",ProjLand);
                        fprintf(fp4,"%12.3f,",PrevProj);
                        fprintf(fp4,"%12.3f,",LandEst);
                        fprintf(fp4,"%10.5f,",FullF[NYears]);

                        if (flag)
                        {
                            fprintf(fp2," *\n");
                            fprintf(fp4,"*,\n");
                        }
                        else
                        {
                            fprintf(fp2,"\n");
                            fprintf(fp4,"\n");
                        }


                        PrevLand = LandEst;

                        PrevProj = ProjLand;

                    }

                    NYears++;

                }


                fprintf(fp2,"\n");
                BasicReport();
                AuxiliaryReport();

                NYears = KYears;
            }

            fclose(fp2);
            fclose(fp3);
            fclose(fp4);

            printf("---- Run Completed ----\n");

        }

        void ReadInputDataFile(char *fn)
        {
            char buffer[MAXBUF];
            char *c;
            char *t;
            long i,j,k;


            if ((fp1 = fopen(fn,"r")) == NULL)
            {
                fprintf(stderr,"Unable to Open Input File: %s\n",fn);
                exit(1);
            }

            fgets(buffer,MAXBUF-1,fp1);
            if (!strstr(buffer,"MSE VERSION 4.1"))
            {
                fprintf(stderr,"Not Current PopSim Input Data File\n");
                exit(1);
            }

            DiscFlag = 0;

            while (!feof(fp1))
            {
                fgets(buffer,MAXBUF-1,fp1);

                if (strstr(buffer,"GENERAL"))
                {
                    fgets(buffer,MAXBUF,fp1);
                    t = strtok(buffer," \t\r\n");
                    NFYear = atol(t);
                    t = strtok(NULL," \t\r\n");
                    NXYear = atol(t);
                    t = strtok(NULL," \t\r\n");
                    MXYear = atol(t);
                    t = strtok(NULL," \t\r\n");
                    MaxAge = atol(t);
                    t = strtok(NULL," \t\r\n");
                    NAgePlus = atol(t);
                    t = strtok(NULL," \t\r\n");
                    NSurvey = atol(t);
                    t = strtok(NULL," \t\r\n");
                    NIter = atol(t);
                    t = strtok(NULL," \t\r\n");
                    ISeed = atol(t);
                    t = strtok(NULL," \t\r\n");
                    ModelFlag = atol(t);

                    AllocData();
                }
                else if (strstr(buffer, "ASSESSMENT_FREQUENCY"))
                {
                    fgets(buffer, MAXBUF - 1, fp1);
                    t = strtok(buffer, " \t\r\n");
                    this->assessment_frequency = atof(t);
                }
                else if (strstr(buffer,"CASEID"))
                {
                    fgets(buffer,TXTBUF,fp1);
                    strcpy(CaseID,buffer);

                    if ((c = strrchr(CaseID,'\n')) !=NULL)
                        *c = '\0';
                }
                else if (strstr(buffer,"INITPOP"))
                {
                    fgets(buffer,MAXBUF,fp1);
                    t = strtok(buffer," \t\r\n");
                    CvPop = atof(t);
                    for (i =0; i < MaxAge; i++)
                    {
                        fgets(buffer,MAXBUF,fp1);
                        t = strtok(buffer," \t\r\n");
                        InitPop[i] = atof(t);
                    }
                }
                else if (strstr(buffer,"NATMORT"))
                {
                    for (j = 0; j < MaxYears; j++)
                    {
                        fgets(buffer,MAXBUF,fp1);
                        t = strtok(buffer," \t\r\n");
                        for (i = 0; i < MaxAge; i++)
                        {
                            NatMort[i][j] = atof(t);
                            t = strtok(NULL," \t\r\n");
                        }
                    }
                }
                else if (strstr(buffer,"MATURE"))
                {
                    for (j = 0; j < MaxYears; j++)
                    {
                        fgets(buffer,MAXBUF,fp1);
                        t = strtok(buffer," \t\r\n");
                        for (i = 0; i < MaxAge; i++)
                        {
                            Mature[i][j] = atof(t);
                            t = strtok(NULL," \t\r\n");
                        }
                    }
                }
                else if (strstr(buffer,"TFRAC"))
                {
                    for (i = 0; i < MaxYears; i++)
                    {
                        fgets(buffer,MAXBUF,fp1);
                        t = strtok(buffer," \t\r\n");
                        TF[i] = atof(t);
                        t = strtok(NULL," \t\r\n");
                        TM[i] = atof(t);
                    }
                }
                else if (strstr(buffer,"STOCKWT"))
                {
                    for (j = 0; j < MaxYears; j++)
                    {
                        fgets(buffer,MAXBUF,fp1);
                        t = strtok(buffer," \t\r\n");
                        for (i = 0; i < MaxAge; i++)
                        {
                            StockWeight[i][j] = atof(t);
                            t = strtok(NULL," \t\r\n");
                        }
                    }
                }
                else if (strstr(buffer,"SPAWNWT"))
                {
                    for (j = 0; j < MaxYears; j++)
                    {
                        fgets(buffer,MAXBUF,fp1);
                        t = strtok(buffer," \t\r\n");
                        for (i = 0; i < MaxAge; i++)
                        {
                            SpawnWeight[i][j] = atof(t);
                            t = strtok(NULL," \t\r\n");
                        }
                    }
                }
                else if (strstr(buffer,"CATCHWT"))
                {
                    for (j = 0; j < MaxYears; j++)
                    {
                        fgets(buffer,MAXBUF,fp1);
                        t = strtok(buffer," \t\r\n");
                        for (i = 0; i < MaxAge; i++)
                        {
                            CatchWeight[i][j] = atof(t);
                            t = strtok(NULL," \t\r\n");
                        }
                    }
                }
                else if (strstr(buffer,"DISCFRAC"))
                {
                    DiscFlag = 1;

                    for (j = 0; j < MaxYears; j++)
                    {
                        fgets(buffer,MAXBUF,fp1);
                        t = strtok(buffer," \t\r\n");
                        for (i = 0; i < MaxAge; i++)
                        {
                            DiscFrac[i][j] = atof(t);
                            t = strtok(NULL," \t\r\n");
                        }
                    }
                }
                else if (strstr(buffer,"FSELEC"))
                {
                    for (j = 0; j < MaxYears; j++)
                    {
                        fgets(buffer,MAXBUF,fp1);
                        t = strtok(buffer," \t\r\n");
                        for (i = 0; i < MaxAge; i++)
                        {
                            FSelec[i][j] = atof(t);
                            t = strtok(NULL," \t\r\n");
                        }
                    }
                }
                else if (strstr(buffer,"HARVEST"))
                {
                    for (i = 0; i < NYears; i++)
                    {
                        fgets(buffer,MAXBUF,fp1);
                        t = strtok(buffer," \t\r\n");
                        t = strtok(NULL," \t\r\n");
                        HarvestSpec[i] = atol(t);
                        t = strtok(NULL," \t\r\n");
                        HarvestValue[i] = atof(t);

                        if (HarvestSpec[i] == 0)
                            FullF[i] = HarvestValue[i];
                    }
                }
                else if (strstr(buffer,"RECRUITS"))
                {
                    fgets(buffer,MAXBUF,fp1);
                    t = strtok(buffer," \t\r\n");
                    RecFlag = atol(t);
                    t = strtok(NULL," \t\r\n");
                    RecScale = atof(t);
                    t = strtok(NULL," \t\r\n");
                    SSBScale = atof(t);

                    switch (RecFlag)
                    {
                        case 0:
                            fgets(buffer,MAXBUF,fp1);
                            t = strtok(buffer," \t\r\n");
                            CvRec = atof(t);

                            RecObs1 = (double *) calloc(MaxYears+2,sizeof(double));
                            for (i = 1; i < MaxYears; i++)
                            {
                                fgets(buffer,MAXBUF,fp1);
                                t = strtok(buffer," \t\r\n");
                                RecObs1[i] = atof(t);
                            }

                            RecObs1[MaxYears]   = RecObs1[MaxYears-1];
                            RecObs1[MaxYears+1] = RecObs1[MaxYears-1];

                            break;
                        case 1:
                        case 2:
                            fgets(buffer,MAXBUF,fp1);
                            t = strtok(buffer," \t\r\n");
                            RecAlpha = atof(t);
                            t = strtok(NULL," \t\r\n");
                            RecBeta = atof(t);
                            t = strtok(NULL," \t\r\n");
                            CvRec = atof(t);
                            break;
                        case 3:
                            fgets(buffer,MAXBUF,fp1);
                            t = strtok(buffer," \t\r\n");
                            RecAlpha = atof(t);
                            t = strtok(NULL," \t\r\n");
                            RecBeta = atof(t);
                            t = strtok(NULL," \t\r\n");
                            RecGamma = atof(t);
                            t = strtok(NULL," \t\r\n");
                            CvRec = atof(t);
                            break;
                        case 4:
                            fgets(buffer,MAXBUF,fp1);
                            t = strtok(buffer," \t\r\n");
                            NRec1 = atol(t);

                            RecObs1 = (double *) calloc(NRec1,sizeof(double));
                            for (i = 0; i < NRec1; i++)
                            {
                                fgets(buffer,MAXBUF,fp1);
                                t = strtok(buffer," \t\r\n");
                                RecObs1[i] = atof(t);
                            }
                            break;
                        case 5:
                            fgets(buffer,MAXBUF,fp1);
                            t = strtok(buffer," \t\r\n");
                            NRec1 = atol(t);
                            t = strtok(NULL," \t\r\n");
                            NRec2 = atol(t);
                            t = strtok(NULL," \t\r\n");
                            SSBCutPt = atof(t);

                            RecObs1 = (double *) calloc(NRec1,sizeof(double));
                            for (i = 0; i < NRec1; i++)
                            {
                                fgets(buffer,MAXBUF,fp1);
                                t = strtok(buffer," \t\r\n");
                                RecObs1[i] = atof(t);
                            }
                            fgets(buffer,MAXBUF,fp1);
                            RecObs2 = (double *) calloc(NRec2,sizeof(double));
                            for (i = 0; i < NRec2; i++)
                            {
                                fgets(buffer,MAXBUF,fp1);
                                t = strtok(buffer," \t\r\n");
                                RecObs2[i] = atof(t);
                            }
                            break;
                    }
                }
                else if (strstr(buffer,"CATCH_SAMPLES"))
                {
                    for (i = 0; i < MaxYears; i++)
                    {
                        fgets(buffer,MAXBUF,fp1);
                        t = strtok(buffer," \t\r\n");
                        CatchCV[i] = atof(t);
                        t = strtok(NULL," \t\r\n");
                        CatchDraw[i] = atof(t);
                    }
                }
                else if (strstr(buffer,"SURVEY_SAMPLES"))
                {
                    for (i = 0; i < NSurvey; i++)
                    {
                        for (j = 0; j < MaxYears; j++)
                        {
                            fgets(buffer,MAXBUF,fp1);
                            t = strtok(buffer," \t\r\n");
                            SurveyCV[j][i] = atof(t);
                            t = strtok(NULL," \t\r\n");
                            SurveyQ[j][i] = atof(t);
                            t = strtok(NULL," \t\r\n");
                            SurveyDraw[j][i] = atof(t);
                        }
                    }
                    for (k = 0; k < NSurvey; k++)
                    {
                        for (j = 0; j < MaxYears; j++)
                        {
                            fgets(buffer,MAXBUF,fp1);
                            t = strtok(buffer," \t\r\n");
                            for (i = 0; i < MaxAge; i++)
                            {
                                SurveySelec[i][j][k] = atof(t);
                                t = strtok(NULL," \t\r\n");
                            }
                        }
                    }

                }
                else if (strstr(buffer,"MSR"))
                {
                    fgets(buffer,MAXBUF,fp1);
                    t = strtok(buffer," \t\r\n");
                    Fmsy = atof(t);
                    t = strtok(NULL," \t\r\n");
                    Bmsy = atof(t);
                    t = strtok(NULL, " \t\r\n");
                    Blim = atof(t);

                    fgets(buffer,MAXBUF,fp1);
                    t = strtok(buffer," \t\r\n");
                    cvFEst = atof(t);
                    t = strtok(NULL," \t\r\n");
                    biasFEst = atof(t);

                    fgets(buffer,MAXBUF,fp1);
                    t = strtok(buffer," \t\r\n");
                    cvSSB = atof(t);
                    t = strtok(NULL," \t\r\n");
                    biasSSB = atof(t);

                    fgets(buffer,MAXBUF,fp1);
                    t = strtok(buffer," \t\r\n");
                    cvLand = atof(t);
                    t = strtok(NULL," \t\r\n");
                    biasLand = atof(t);

                    fgets(buffer,MAXBUF,fp1);
                    t = strtok(buffer," \t\r\n");
                    DeltaFlag = atol(t);
                    t = strtok(NULL," \t\r\n");
                    DeltaLand = atof(t);
                }

            }

            fclose(fp1);

        }

        void AllocData() {

            NYears = NXYear - NFYear + 1;
            MaxYears = MXYear - NFYear + 1;

            KYears = NYears;

            InitPop      = (double *) calloc(MaxAge,sizeof(double));
            TF           = (double *) calloc(MaxYears+1,sizeof(double));
            TM           = (double *) calloc(MaxYears+1,sizeof(double));
            Recruits     = (double *) calloc(MaxYears+2,sizeof(double));
            FullF        = (double *) calloc(MaxYears+1,sizeof(double));
            Landings     = (double *) calloc(MaxYears+1,sizeof(double));
            Discards     = (double *) calloc(MaxYears+1,sizeof(double));
            Biomass      = (double *) calloc(MaxYears+1,sizeof(double));
            SSB          = (double *) calloc(MaxYears+1,sizeof(double));
            CatchCV      = (double *) calloc(MaxYears+1,sizeof(double));
            CatchDraw    = (double *) calloc(MaxYears+1,sizeof(double));
            HarvestValue = (double *) calloc(MaxYears+1,sizeof(double));

            ObsRec = (double *) calloc(NRECMAX,sizeof(double));

            HarvestSpec = (long *) calloc(MaxYears+1,sizeof(long));

            StockNumbers = AllocMatrix(MaxAge,MaxYears+2);
            StockBiomass = AllocMatrix(MaxAge,MaxYears+2);

            CatchNumbers = AllocMatrix(MaxAge,MaxYears+1);
            CatchBiomass = AllocMatrix(MaxAge,MaxYears+1);

            SpawnNumbers = AllocMatrix(MaxAge,MaxYears+1);
            SpawnBiomass = AllocMatrix(MaxAge,MaxYears+1);

            NatMort      = AllocMatrix(MaxAge,MaxYears+1);
            FMort        = AllocMatrix(MaxAge,MaxYears+1);
            Mature       = AllocMatrix(MaxAge,MaxYears+1);

            StockWeight  = AllocMatrix(MaxAge,MaxYears+1);
            SpawnWeight  = AllocMatrix(MaxAge,MaxYears+1);
            CatchWeight  = AllocMatrix(MaxAge,MaxYears+1);


            FSelec       = AllocMatrix(MaxAge,MaxYears+1);
            DiscFrac     = AllocMatrix(MaxAge,MaxYears+1);

            CatchAge = (double *) calloc(MaxAge,sizeof(double));
            AgeComp  = (double *) calloc(MaxAge,sizeof(double));


            CatchSamples = AllocMatrix(MaxAge,MaxYears+1);

            SurveyAge = (double *) calloc(MaxAge,sizeof(double));

            SurveyCV   = AllocMatrix(MaxYears+1,NSurvey);
            SurveyDraw = AllocMatrix(MaxYears+1,NSurvey);
            SurveyQ    = AllocMatrix(MaxYears+1,NSurvey);

            SurveySelec   = AllocMat3(MaxAge,MaxYears+1,NSurvey);
            SurveySamples = AllocMat3(MaxAge,MaxYears+1,NSurvey);

            StockEstim = AllocMatrix(NAgePlus,MaxYears+1);
            SSBEstim   = AllocMatrix(NAgePlus,MaxYears+1);
            FEstim     = AllocMatrix(NAgePlus,MaxYears+1);

            KAges = NAgePlus;

            ZeroMatrix(DiscFrac,MaxAge,MaxYears);

            if (ModelFlag == 1)
            {
                StockVPA = AllocMatrix(NAgePlus,MaxYears+1);
                SSBVPA   = AllocMatrix(NAgePlus,MaxYears+1);
                FVPA     = AllocMatrix(NAgePlus,MaxYears+1);

                AgeProCV = AllocMatrix(NAgePlus,AGEPROCV);
            }
            else if (ModelFlag  == 2)
            {
                AgeProCV = AllocMatrix(NAgePlus,AGEPROCV);
            }
        }

        void CalcStockBiomass(long k)
        {
            long i;

            Biomass[k] = 0.0;

            for (i = 0; i < MaxAge; i++)
            {
                StockBiomass[i][k] = StockNumbers[i][k] * StockWeight[i][k];
                Biomass[k] += StockBiomass[i][k] / 1000.0;
            }
        }

        void CalcCatch(long k)
        {

            if (HarvestSpec[k])
                CalcHarvestF(k);
            else
                ApplyHarvestF(k);
        }

        void ApplyHarvestF(long k)
        {
            long i;
            double m;
            double f;
            double z;
            double zx;


            for (i = 0; i < MaxAge; i++)
            {
                StockNumbers[i][k+1]= 0.0;
                SpawnNumbers[i][k]= 0.0;
            }

            for (i = 0; i < MaxAge-1; i++)
            {

                m = NatMort[i][k];

                f = FullF[k] * FSelec[i][k];

                z = m + f;

                StockNumbers[i+1][k+1] = StockNumbers[i][k] * exp(-z);
            }


            Landings[k] = 0.0;
            Discards[k] = 0.0;

            SSB[k] = 0.0;

            for (i = 0; i < MaxAge; i++)
            {

                m = NatMort[i][k];

                f = FullF[k] * FSelec[i][k];

                z = m + f;

                FMort[i][k] = f;

                CatchNumbers[i][k] = StockNumbers[i][k] * f * (1.0 - exp(-z)) / z;

                CatchBiomass[i][k] = CatchNumbers[i][k] * CatchWeight[i][k];

                Landings[k] += CatchBiomass[i][k] * (1.0 - DiscFrac[i][k]) / 1000.0;

                Discards[k] += CatchBiomass[i][k] * DiscFrac[i][k] / 1000.0;

                zx = m * TM[k] + f * TF[k];

                SpawnNumbers[i][k] = Mature[i][k] * StockNumbers[i][k] * exp(-zx);

                SpawnBiomass[i][k] = SpawnNumbers[i][k] * SpawnWeight[i][k];

                SSB[k] += SpawnBiomass[i][k] / 1000.0;

            }

        }

        void CalcHarvestF(long k)
        {
            long iter;
            double err, xx;
            double hcalc, htarget;
            double xlo, xhi;

            /* Initial Guess */

            xx = 0.5;

            xlo = XLBOUND;
            xhi = XUBOUND;

            /* Trial & Error */

            iter = 0;

            htarget = HarvestValue[k];

            while (iter < MXITER)
            {

                FullF[k] = xx;

                hcalc = 0.0;

                ApplyHarvestF(k);

                if (HarvestSpec[k] == 1)
                    hcalc = Landings[k];
                else
                    hcalc = Landings[k] + Discards[k];

                err = (hcalc - htarget) / htarget;

                if (fabs(err) < XTOL)
                    return;

                if (hcalc < htarget)
                {
                    xlo = xx;
                    xx = (xhi + xx) * 0.5;

                }
                else
                {
                    xhi = xx;
                    xx = (xx + xlo) * 0.5;
                }

                iter++;
            }
            if (iter == MXITER)
                fprintf(stderr,"Maximum Iterations Exceeded in Harvest\n");

        }

        void ApplyRecruitment(long k)
        {
            double ssb;

            switch(RecFlag)
            {
                case 0:
                    Recruits[k+1] = UserRecruitment(k);
                    break;
                case 1:
                    ssb = SSB[k];
                    Recruits[k+1] = BevertonHolt(ssb);
                    break;
                case 2:
                    ssb = SSB[k];
                    Recruits[k+1] = Ricker(ssb);
                    break;
                case 3:
                    ssb = SSB[k];
                    Recruits[k+1] = Shepherd(ssb);
                    break;
                case 4:
                    Recruits[k+1] = EmpiricalRecruits1();
                    break;
                case 5:
                    ssb = SSB[k];
                    Recruits[k+1] = EmpiricalRecruits2(ssb);
                    break;
            }

            StockNumbers[0][k+1] = Recruits[k+1];
        }

        double UserRecruitment(long k)
        {

            double sd;
            double zscore;
            double xx, xr;

            xx = RecObs1[k+1];

            sd = sqrt(log(CvRec * CvRec + 1.0));

            zscore = gsl_box_muller();

            xx = xx * exp(sd * zscore - sd * sd * 0.5);

            xr = xx * RecScale;

            return(xr);
        }

        double BevertonHolt(double ssb)
        {
            double alpha, beta;
            double x, xx, xr;
            double sd;
            double zscore;

            ssb = ssb / SSBScale;

            alpha = RecAlpha;
            beta  = RecBeta;

            x = alpha * ssb / (beta + ssb);

            sd = sqrt(log(CvRec * CvRec + 1.0));

            do
            {

                /* Get Zscore from Box-muller Transformation*/

                zscore = gsl_box_muller();


                xx = x * exp(zscore * sd - sd * sd * 0.5);

            } while (xx < 0.001);

            xr = xx * RecScale;

            return(xr);
        }

        double Ricker(double ssb)
        {
            double alpha, beta;
            double x, xx, xr;
            double sd;
            double zscore;

            ssb = ssb / SSBScale;

            alpha = RecAlpha;
            beta  = RecBeta;

            x = alpha * ssb * exp(-beta*ssb);

            sd = sqrt(log(CvRec * CvRec + 1.0));

            do
            {

                /* Get Zscore from Box muller Transformation */

                zscore = gsl_box_muller();


                xx = x * exp(zscore * sd - sd * sd * 0.5);

            } while (xx < 0.001);

            xr = xx * RecScale;

            return(xr);

        }

        double Shepherd(double ssb)
        {
            double alpha, beta, kparm;
            double x, xx, xr;
            double sd;
            double zscore;

            ssb = ssb / SSBScale;

            alpha = RecAlpha;
            beta  = RecBeta;
            kparm = RecGamma;

            xx = ssb / kparm;

            x = alpha * ssb / (1.0 + pow(xx,beta));

            sd = sqrt(log(CvRec * CvRec + 1.0));

            do
            {

                /* Get Zscore from Box-muller Transformation */

                zscore = gsl_box_muller();

                xx = x * exp(zscore * sd - sd * sd * 0.5);

            } while (xx < 0.001);

            xr = xx * RecScale;

            return(xr);

        }

        double EmpiricalRecruits1()
        {
            long i, k;
            double xr;


            xr = 0.0;

            k = NRec1;

            for (i = 0; i < k; i++)
                ObsRec[i] = RecObs1[i];

            xr = EmpiricalCDF(k,ObsRec);


            xr = xr * RecScale;


            return(xr);
        }

        double EmpiricalRecruits2(double ssb)
        {
            long i, k;
            double xr;


            ssb = ssb / SSBScale;


            if (ssb > SSBCutPt)
            {

                k = NRec2;

                for (i = 0; i < k; i++)
                    ObsRec[i] = RecObs2[i];
            }
            else
            {
                k = NRec1;

                for (i = 0; i < k; i++)
                    ObsRec[i] = RecObs1[i];
            }


            xr = 0.0;


            xr = EmpiricalCDF(k,ObsRec);

            xr = xr * RecScale;


            return(xr);

        }

        double EmpiricalCDF(long n,double *x)
        {
            long k;
            double xn, xk;
            double xr, xx, rnum;

            /* Vector Must be sorted in ascending order */

            BubbleUpSort(n,x);

            xn = (double) (n - 1);

            /* Get Uniform Pseudorandom Number */

            rnum = DRNUNF();

            xk = floor(rnum * xn);

            k = (long) xk;

            xx = xk / xn;


            xr = x[k] + (x[k+1] - x[k]) * (rnum - xx) * xn;

            return xr;

        }

        void AuxiliaryReport()
        {
            long i,j;
            double xt;


            fprintf(fp3,"Year,");
            for (i = 0; i < MaxYears; i++)
                fprintf(fp3,"%ld,",NFYear+i);
            fprintf(fp3,"\n");

            fprintf(fp3,"Stock Numbers,");
            for (j = 0; j < MaxYears; j++)
            {
                xt = 0.0;

                for (i = 0; i < MaxAge; i++)
                    xt += StockNumbers[i][j];

                fprintf(fp3,"%-18.8E, ",xt);
            }
            fprintf(fp3,"\n");

            fprintf(fp3,"Catch Numbers,");
            for (j = 0; j < MaxYears; j++)
            {
                xt = 0.0;

                for (i = 0; i < MaxAge; i++)
                    xt += CatchNumbers[i][j];

                fprintf(fp3,"%-18.8E, ",xt);
            }
            fprintf(fp3,"\n");

            fprintf(fp3,"Stock Biomass (MT),");
            for (j = 0; j < MaxYears; j++)
                fprintf(fp3,"%-18.8E, ",Biomass[j]);
            fprintf(fp3,"\n");

            fprintf(fp3,"Spawning Stock Biomass (MT),");
            for (j = 0; j < MaxYears; j++)
                fprintf(fp3,"%-18.8E, ",SSB[j]);
            fprintf(fp3,"\n");

            fprintf(fp3,"Recruits,");
            for (j = 0; j < MaxYears; j++)
                fprintf(fp3,"%-18.8E, ",Recruits[j]);
            fprintf(fp3,"\n");

            fprintf(fp3,"Landings (MT),");
            for (j = 0; j < MaxYears; j++)
                fprintf(fp3,"%-18.8E, ",Landings[j]);
            fprintf(fp3,"\n");

            fprintf(fp3,"Discards (MT),");
            for (j = 0; j < MaxYears; j++)
                fprintf(fp3,"%-18.8E, ",Discards[j]);
            fprintf(fp3,"\n");


            fprintf(fp3,"Fully Recruited F,");
            for (j = 0; j < MaxYears; j++)
                fprintf(fp3,"%-18.8E, ",FullF[j]);
            fprintf(fp3,"\n");
        }

        void BasicReport()
        {
            long i,j,k;
            double xt;

            /* Stock Numbers */

            fprintf(fp2,"Stock Numbers\n");
            fprintf(fp2,"Year                  ");

            for (i = 0; i < MaxYears; i++)
                fprintf(fp2,"%4ld            ",NFYear+i);
            fprintf(fp2,"\n\n");

            for (i = 0; i < NAgePlus-1; i++)
            {
                fprintf(fp2,"Age %3ld    ",i+1);
                for (j = 0; j < MaxYears; j++)
                    fprintf(fp2,"%15.0f ",StockNumbers[i][j]);
                fprintf(fp2,"\n");
            }

            fprintf(fp2,"Age %3ld+   ",NAgePlus);
            for (j = 0; j < MaxYears; j++)
            {
                xt = 0.0;
                for (k = NAgePlus-1; k < MaxAge; k++)
                    xt += StockNumbers[k][j];

                fprintf(fp2,"%15.0f ",xt);
            }
            fprintf(fp2,"\n");
            fprintf(fp2,"\n");


            fprintf(fp2,"Total      ");
            for (j = 0; j < MaxYears; j++)
            {
                xt = 0.0;
                for (k = 0; k < MaxAge; k++)
                    xt += StockNumbers[k][j];

                fprintf(fp2,"%15.0f ",xt);
            }
            fprintf(fp2,"\n");
            fprintf(fp2,"\n");

            /* Catch Numbers */

            fprintf(fp2,"Catch Numbers\n");
            fprintf(fp2,"Year                  ");

            for (i = 0; i < MaxYears; i++)
                fprintf(fp2,"%4ld            ",NFYear+i);
            fprintf(fp2,"\n\n");

            for (i = 0; i < NAgePlus-1; i++)
            {
                fprintf(fp2,"Age %3ld    ",i+1);
                for (j = 0; j < MaxYears; j++)
                    fprintf(fp2,"%15.0f ",CatchNumbers[i][j]);
                fprintf(fp2,"\n");
            }

            fprintf(fp2,"Age %3ld+   ",NAgePlus);
            for (j = 0; j < MaxYears; j++)
            {
                xt = 0.0;
                for (k = NAgePlus-1; k < MaxAge; k++)
                    xt += CatchNumbers[k][j];

                fprintf(fp2,"%15.0f ",xt);
            }
            fprintf(fp2,"\n");
            fprintf(fp2,"\n");


            fprintf(fp2,"Total      ");
            for (j = 0; j < MaxYears; j++)
            {
                xt = 0.0;
                for (k = 0; k < MaxAge; k++)
                    xt += CatchNumbers[k][j];

                fprintf(fp2,"%15.0f ",xt);
            }
            fprintf(fp2,"\n");
            fprintf(fp2,"\n");

            /* Spawning Stock Numbers */

            fprintf(fp2,"Spawning Stock Numbers\n");
            fprintf(fp2,"Year                  ");

            for (i = 0; i < MaxYears; i++)
                fprintf(fp2,"%4ld            ",NFYear+i);
            fprintf(fp2,"\n\n");

            for (i = 0; i < NAgePlus-1; i++)
            {
                fprintf(fp2,"Age %3ld    ",i+1);
                for (j = 0; j < MaxYears; j++)
                    fprintf(fp2,"%15.0f ",SpawnNumbers[i][j]);
                fprintf(fp2,"\n");
            }

            fprintf(fp2,"Age %3ld+   ",NAgePlus);
            for (j = 0; j < MaxYears; j++)
            {
                xt = 0.0;
                for (k = NAgePlus-1; k < MaxAge; k++)
                    xt += SpawnNumbers[k][j];

                fprintf(fp2,"%15.0f ",xt);
            }
            fprintf(fp2,"\n");
            fprintf(fp2,"\n");


            fprintf(fp2,"Total      ");
            for (j = 0; j < MaxYears; j++)
            {
                xt = 0.0;
                for (k = 0; k < MaxAge; k++)
                    xt += SpawnNumbers[k][j];

                fprintf(fp2,"%15.0f ",xt);
            }
            fprintf(fp2,"\n");
            fprintf(fp2,"\n");

            /* Landings & Discards */

            fprintf(fp2,"Year                  ");

            for (i = 0; i < MaxYears; i++)
                fprintf(fp2,"%4ld            ",NFYear+i);
            fprintf(fp2,"\n\n");

            fprintf(fp2,"Landings (MT)");
            for (j = 0; j < MaxYears; j++)
                fprintf(fp2,"%15.3f ",Landings[j]);
            fprintf(fp2,"\n");

            fprintf(fp2,"Discards (MT)");
            for (j = 0; j < MaxYears; j++)
                fprintf(fp2,"%15.3f ",Discards[j]);
            fprintf(fp2,"\n\n");

            fprintf(fp2,"Total (MT)   ");
            for (j = 0; j < MaxYears; j++)
                fprintf(fp2,"%15.3f ",Landings[j]+Discards[j]);
            fprintf(fp2,"\n\n");

            /* Biomass */

            fprintf(fp2,"Year                  ");

            for (i = 0; i < MaxYears; i++)
                fprintf(fp2,"%4ld            ",NFYear+i);
            fprintf(fp2,"\n\n");

            fprintf(fp2,"Stock Bio (MT)");
            for (j = 0; j < MaxYears; j++)
                fprintf(fp2,"%15.3f ",Biomass[j]);
            fprintf(fp2,"\n\n");

            fprintf(fp2,"Year                  ");

            for (i = 0; i < MaxYears; i++)
                fprintf(fp2,"%4ld            ",NFYear+i);
            fprintf(fp2,"\n\n");
            fprintf(fp2,"SSB (MT)     ");
            for (j = 0; j < MaxYears; j++)
                fprintf(fp2,"%15.3f ",SSB[j]);
            fprintf(fp2,"\n");

            fprintf(fp2,"Recruits     ");
            for (j = 0; j < MaxYears; j++)
                fprintf(fp2,"%15.0f ",Recruits[j]);
            fprintf(fp2,"\n\n");


            fprintf(fp2,"Fishing Mortality\n\n");

            fprintf(fp2,"Year                  ");

            for (i = 0; i < MaxYears; i++)
                fprintf(fp2,"%4ld            ",NFYear+i);
            fprintf(fp2,"\n\n");

            fprintf(fp2,"Full-F       ");
            for (j = 0; j < MaxYears; j++)
                fprintf(fp2,"%15.6f ",FullF[j]);
            fprintf(fp2,"\n\n");

            fprintf(fp2,"\n");

        }

        void DrawCatchSamples()
        {
            long i, j, n;
            long ka;
            double x, xx, xxx;
            double cv, sd, zscore;


            for (j = 0; j < NYears; j++)
            {

                x = 0.0;

                for (i = 0; i < MaxAge; i++)
                    x += CatchNumbers[i][j];

                cv = CatchCV[j];
                sd = sqrt(log(cv * cv + 1.0));

                do
                {

                    /* Get Zscore from Box muller Transformation */

                    zscore = gsl_box_muller();


                    xx = x * exp(zscore * sd - sd * sd * 0.5);

                } while (xx < 0.001);



                for (i = 0; i < MaxAge; i++)
                    CatchAge[i] = CatchNumbers[i][j] * xx / x;


                AgeComp[0] = CatchAge[0] / xx;

                for (i = 1; i < MaxAge; i++)
                    AgeComp[i] = AgeComp[i-1] + CatchAge[i] / xx;

                for (i = 0; i < MaxAge; i++)
                    CatchSamples[i][j] = 0.0;

                n = (long) floor(CatchDraw[j]);

                for (i = 0; i < n; i++)
                {
                    ka = DrawSample();
                    CatchSamples[ka][j]++;
                }

                xxx = 0.0;
                for (i = 0; i < MaxAge; i++)
                {
                    CatchSamples[i][j] = CatchSamples[i][j] * xx / (double) n;
                    xxx += CatchSamples[i][j];
                }

            }
        }

        long DrawSample()
        {
            long i, k;
            double rnum;

            rnum = DRNUNF();

            for (i = 0; i < MaxAge; i++)
            {
                if (rnum < AgeComp[i])
                {
                    k = i;
                    break;
                }
            }

            return(k);
        }

        void DrawSurveySamples()
        {
            long i, j, k, n;
            long ka;
            double x, xx;
            double cv, sd, zscore;


            for (j = 0; j < NYears+1; j++)
            {
                for (k = 0; k < NSurvey; k++)
                {

                    for (i = 0; i < MaxAge; i++)
                        SurveyAge[i] = StockNumbers[i][j] * SurveyQ[j][k] * SurveySelec[i][j][k];

                    x = 0.0;

                    for (i = 0; i < MaxAge; i++)
                        x += SurveyAge[i];

                    if (x > 0.0)
                    {

                        cv = SurveyCV[j][k];
                        sd = sqrt(log(cv * cv + 1.0));

                        do
                        {

                            /* Get Zscore from Box muller Transformation */

                            zscore = gsl_box_muller();


                            xx = x * exp(zscore * sd - sd * sd * 0.5);

                        } while (xx < 0.001);



                        for (i = 0; i < MaxAge; i++)
                            SurveyAge[i] = SurveyAge[i] * xx / x;


                        AgeComp[0] = SurveyAge[0] / xx;

                        for (i = 1; i < MaxAge; i++)
                            AgeComp[i] = AgeComp[i-1] + SurveyAge[i] / xx;

                        for (i = 0; i < MaxAge; i++)
                            SurveySamples[i][j][k] = 0.0;

                        n = (long) floor(SurveyDraw[j][k]);

                        for (i = 0; i < n; i++)
                        {
                            ka = DrawSample();
                            SurveySamples[ka][j][k]++;
                        }

                        for (i = 0; i < MaxAge; i++)
                            SurveySamples[i][j][k] = SurveySamples[i][j][k] * xx / (double) n;
                    }
                    else
                    {
                        for (i = 0; i < MaxAge; i++)
                            SurveySamples[i][j][k] = 0.0;
                    }

                }
            }
        }

        double MSRule(double B)
        {
            double FTarg;

            if (B < Blim)
                FTarg = 0.0;
            else if (B < Bmsy)
                FTarg = Fmsy * B / Bmsy;
            else
                FTarg = Fmsy;

            return(FTarg);

        }

        void ReadVPATemplateFile(char *fn)
        {
            char xname[FILBUF];
            char buffer[MAXBUF];
            char *c;
            char *t;
            long i, j;


            strcpy(xname,fn);

            if ((c = strrchr(xname,'.')) != NULL)
                *c = '\0';

            strcat(xname,".tx1");

            if ((fp1 = fopen(xname,"r")) == NULL)
            {
                fprintf(stderr,"Unable to Open Template File: %s\n",xname);
                exit(1);
            }


            PRVec    = (double *) calloc(KAges,sizeof(double));
            AgeEst   = (long *) calloc(KAges,sizeof(short));
            StockEst = (double *) calloc(KAges,sizeof(double));

            AvgFVPA  = (double *) calloc(MaxYears+1,sizeof(double));

            NMortMat = AllocMatrix(KAges,MaxYears+1);

            CatchFill  = 0.1;

            while (!feof(fp1))
            {
                fgets(buffer,MAXBUF,fp1);

                if (strstr(buffer,"PATH"))
                {
                    fgets(buffer,MAXBUF,fp1);
                    if ((c = strrchr(buffer,'\n')) != NULL)
                        *c = '\0';
                    strcpy(ModelPath,buffer);
                }
                else if (strstr(buffer,"INDEX"))
                {
                    fgets(buffer,MAXBUF,fp1);
                    t = strtok(buffer," \t\r\n");
                    NIndex = atol(t);

                    IndexValues = AllocMatrix(MaxYears+1,NIndex);

                    IndexData = (struct si *) calloc(NIndex,(sizeof(struct si)));
                    for (i = 0; i < NIndex; i++)
                    {
                        fgets(buffer,MAXBUF,fp1);
                        t = strtok(buffer," \t\r\n");
                        IndexData[i].StartAge = atoi(t);
                        t = strtok(NULL," \t\r\n");
                        IndexData[i].EndAge = atoi(t);
                        t = strtok(NULL," \t\r\n");
                        IndexData[i].Type = atoi(t);
                        t = strtok(NULL," \t\r\n");
                        IndexData[i].SurveyIndex = atoi(t);
                    }
                }
                else if (strstr(buffer,"FOLD"))
                {
                    fgets(buffer,MAXBUF,fp1);
                    t = strtok(buffer," \t\r\n");
                    Fold[0] = atoi(t);
                    t = strtok(NULL," \t\r\n");
                    Fold[1] = atoi(t);
                }
                else if (strstr(buffer,"AVGF"))
                {
                    fgets(buffer,MAXBUF,fp1);
                    t = strtok(buffer," \t\r\n");
                    AvgF[0] = atoi(t);
                    t = strtok(NULL," \t\r\n");
                    AvgF[1] = atoi(t);
                }
                else if (strstr(buffer,"PRVEC"))
                {
                    fgets(buffer,MAXBUF,fp1);
                    t = strtok(buffer," \t\r\n");
                    for (i = 0; i < KAges-1; i++)
                    {
                        PRVec[i] = atof(t);
                        t = strtok(NULL," \t\r\n");
                    }
                }
                else if (strstr(buffer,"STOCKEST"))
                {
                    fgets(buffer,MAXBUF,fp1);
                    t = strtok(buffer," \t\r\n");
                    for (i = 0; i < KAges; i++)
                    {
                        StockEst[i] = atof(t);
                        t = strtok(NULL," \t\r\n");
                    }
                }
                else if (strstr(buffer,"OPTIONS"))
                {
                    fgets(buffer,MAXBUF,fp1);
                    t = strtok(buffer," \t\r\n");
                    for (i = 0; i < NOPT; i++)
                    {
                        OptionFlags[i] = atoi(t);
                        t = strtok(NULL," \t\r\n");
                    }
                }
                else if (strstr(buffer,"RANGE"))
                {
                    fgets(buffer,MAXBUF,fp1);
                    t = strtok(buffer," \t\r\n");
                    NRange[0] = atoi(t);
                    t = strtok(NULL," \t\r\n");
                    NRange[1] = atoi(t);
                }
                else if (strstr(buffer,"NATMORT"))
                {
                    for (j = 0; j < MaxYears; j++)
                    {
                        fgets(buffer,MAXBUF,fp1);
                        t = strtok(buffer," \t\r\n");
                        for (i = 0; i < KAges; i++)
                        {
                            NMortMat[i][j] = atof(t);
                            t = strtok(NULL," \t\r\n");
                        }
                    }
                    for (i = 0; i < KAges; i++)
                        NMortMat[i][MaxYears] = NMortMat[i][MaxYears-1];
                }
                else if (strstr(buffer,"BOOTSTRAP"))
                {
                    fgets(buffer,MAXBUF,fp1);
                    t = strtok(buffer," \t\r\n");
                    NBoot = atol(t);
                }
            }

            fclose(fp1);

        }
        void WriteVPAInputFile(char *fn)
        {
            char xname[FILBUF];
            char *c;
            long i, j, k, n;
            long ny1;
            long KAgest;
            long nfage = 1;
            long k1, k2;
            double xx;
            double fplus = 1.0;
            double stkmin = 1.0;
            double stkmax = 1.0E+08;

            strcpy(xname,fn);

            if ((c = strrchr(xname,'.')) != NULL)
                *c = '\0';

            strcat(xname,"_vpa.dat");

            if ((fp1 = fopen(xname,"w")) == NULL)
            {
                fprintf(stderr,"Unable to Open VPA File: %s\n",xname);
                exit(1);
            }

            ny1 = NYears + 1;

            for (i = 0; i < KAges; i++)
                AgeEst[i] = -1;

            KAgest = 0;
            for (i = 0; i < KAges; i++)
            {
                if (StockEst[i] > 0.0)
                {
                    AgeEst[i] = i + 1;
                    KAgest++;
                }
            }

            for (j = 0; j < NIndex; j++)
            {
                k1 = IndexData[j].StartAge;
                k2 = IndexData[j].EndAge;
                n  = IndexData[j].SurveyIndex ;

                for (i = 0; i < ny1; i++)
                {
                    xx = 0.0;

                    for (k = k1; k <= k2; k++)
                    {
                        if (IndexData[j].Type == 0)
                            xx += SurveySamples[k-1][i][n-1];
                        else
                            xx += SurveySamples[k-1][i][n-1] * StockWeight[k-1][i];
                    }

                    IndexValues[i][j] = xx;

                }
            }

            fprintf(fp1,"VPA/ADAPT V 3.0\n");
            fprintf(fp1,"##\n");
            fprintf(fp1,"MODEL ID\n");
            fprintf(fp1,"Test Case\n");
            fprintf(fp1,"PARAM\n");
            fprintf(fp1,"%ld ",NYears);
            fprintf(fp1,"%ld ",KAges);
            fprintf(fp1,"%ld ",NIndex);
            fprintf(fp1,"%ld ",NIndex);
            fprintf(fp1,"%ld ",KAgest);
            fprintf(fp1,"%ld ",NFYear);
            fprintf(fp1,"%ld ",nfage);
            fprintf(fp1,"%ld ",NBoot);
            fprintf(fp1,"%ld\n",ISeed);

            fprintf(fp1,"PARTIAL RECRUIT\n");
            for (i = 0; i < KAges-1; i++)
                fprintf(fp1,"%-6.4f ",PRVec[i]);
            fprintf(fp1,"\n");

            fprintf(fp1,"AGE ESTIMATE\n");
            for (i = 0; i < KAges; i++)
            {
                if (AgeEst[i] > -1)
                    fprintf(fp1,"%-2d ",AgeEst[i]);
            }
            fprintf(fp1,"\n");

            fprintf(fp1,"STOCK ESTIMATE\n");
            for (i = 0; i < KAges; i++)
            {
                if (AgeEst[i] > -1)
                    fprintf(fp1,"%10.0f ",StockEst[i]);
            }
            fprintf(fp1,"\n");

            fprintf(fp1,"F-PLUS\n");
            fprintf(fp1,"%-6.4f*%-ld\n",fplus,NYears);

            fprintf(fp1,"STOCK MIN-MAX\n");
            fprintf(fp1,"%-8.2E %-8.2E\n",stkmin,stkmax);

            fprintf(fp1,"AGES-F\n");
            fprintf(fp1,"%-d-%-d\n",Fold[0],Fold[1]);

            fprintf(fp1,"AGES-SUMMARY\n");
            fprintf(fp1,"%-d-%-d\n",AvgF[0],AvgF[1]);

            fprintf(fp1,"MFSPAWN\n");
            fprintf(fp1,"%-6.4f  %-6.4f\n",TM[NYears-1],TF[NYears-1]);

            fprintf(fp1,"CATCH AT AGE\n");
            for (j = 0; j < NYears; j++)
            {
                for (i = 0; i < NAgePlus-1; i++)
                {
                    xx = CatchSamples[i][j];

                    if (xx < CatchFill)
                        xx = CatchFill;

                    fprintf(fp1,"%-12.3f ",xx);
                }

                xx = 0.0;
                for (i = NAgePlus-1; i < MaxAge; i++)
                    xx += CatchSamples[i][j];

                if (xx < CatchFill)
                    xx = CatchFill;

                fprintf(fp1,"%-12.3f\n",xx);

            }

            fprintf(fp1,"WEIGHT AT AGE\n");
            for (i = 0; i < NYears; i++)
            {
                for (j = 0; j < KAges; j++)
                    fprintf(fp1,"%-12.3f ",CatchWeight[j][i]);
                fprintf(fp1,"\n");
            }

            fprintf(fp1,"BIOMASS\n");
            for (i = 0; i < NYears+1; i++)
            {
                for (j = 0; j < KAges; j++)
                    fprintf(fp1,"%-12.3f ",StockWeight[j][i]);
                fprintf(fp1,"\n");
            }

            fprintf(fp1,"SSB\n");
            for (i = 0; i < NYears; i++)
            {
                for (j = 0; j < KAges; j++)
                    fprintf(fp1,"%-12.3f ",SpawnWeight[j][i]);
                fprintf(fp1,"\n");
            }

            fprintf(fp1,"MATURITY\n");
            for (i = 0; i < NYears; i++)
            {
                for (j = 0; j < KAges; j++)
                    fprintf(fp1,"%-12.3f ",Mature[j][i]);
                fprintf(fp1,"\n");
            }

            fprintf(fp1,"M MATRIX\n");
            for (i = 0; i < NYears; i++)
            {
                for (j = 0; j < KAges; j++)
                    fprintf(fp1,"%-12.3f ",NMortMat[j][i]);
                fprintf(fp1,"\n");
            }

            fprintf(fp1,"SURVEY INDEX\n");
            for (i = 0; i < NIndex; i++)
                fprintf(fp1,"INDEX_%2.2ld        ",i+1);
            fprintf(fp1,"\n");

            for (i = 0; i < NIndex; i++)
                if (IndexData[i].StartAge == IndexData[i].EndAge)
                    fprintf(fp1,"%-15d ",IndexData[i].StartAge);
                else if (IndexData[i].StartAge < NAgePlus && IndexData[i].EndAge <= NAgePlus)
                    fprintf(fp1,"%7d:%-7d ",IndexData[i].StartAge,IndexData[i].EndAge);
                else
                    fprintf(fp1,"%-15d ",NAgePlus);
            fprintf(fp1,"\n");

            for (i = 0; i < NIndex; i++)
                fprintf(fp1,"1-Jan           ");
            fprintf(fp1,"\n");

            for (i = 0; i < NIndex; i++)
                if (IndexData[i].Type == 0)
                    fprintf(fp1,"Number          ");
                else
                    fprintf(fp1,"Weight          ");
            fprintf(fp1,"\n");

            for (j = 0; j < NIndex; j++)
            {
                xx = 0.0;

                for (i = 0; i < NYears+1; i++)
                    xx += IndexValues[i][j];

                if (xx < XTOL)
                    IndexValues[0][j] = 0.1;
            }


            for (i = 0; i < NYears+1; i++)
            {
                for (j = 0; j < NIndex; j++)
                    fprintf(fp1,"%-15.7E ",IndexValues[i][j]);

                fprintf(fp1,"\n");
            }

            fprintf(fp1,"CHECKED INDEX\n");
            for (i = 0; i < NIndex; i++)
                fprintf(fp1,"%-d ",i+1);
            fprintf(fp1,"\n");

            fprintf(fp1,"CHECKED RETRO\n");
            for (i = 0; i < NIndex; i++)
                fprintf(fp1,"%-d ",1);
            fprintf(fp1,"\n");

            fprintf(fp1,"OPTIONS\n");
            for (i = 0; i < NOPT; i++)
                fprintf(fp1,"%d  ",OptionFlags[i]);
            fprintf(fp1,"\n");

            if (AgeEst[0] != 1)
            {
                fprintf(fp1,"YEAR_RANGE\n");
                fprintf(fp1,"%d-%d\n",NRange[0],NRange[1]);
            }

            fprintf(fp1,"XPARM\n");
            fprintf(fp1,"80.0  1.0\n" );

            fclose(fp1);
        }

        void ReadAgeProTemplateFile(char *fn)
        {
            char xname[FILBUF];
            char buffer[MAXBUF];
            char *c;
            char *tok;
            long i, j;

            strcpy(xname,fn);

            if ((c = strrchr(xname,'.')) != NULL)
                *c = '\0';

            strcat(xname,".tpx");

            if ((fp6 = fopen(xname,"r")) == NULL)
            {
                fprintf(stderr,"Unable to Open AgePro Template File: %s\n",xname);
                exit(1);
            }

            fgets(AgeProPath,MAXBUF,fp6);

            if ((c = strchr(AgeProPath,'\n')) != NULL)
                *c = '\0';

            fgets(buffer,MAXBUF,fp6);
            tok = strtok(buffer," \t\r\n");
            NAgeProSim = atoi(tok);

            fgets(buffer,MAXBUF,fp6);
            tok = strtok(buffer," \t\r\n");
            cvAgeProErr = atof(tok);

            tok = strtok(NULL," \t\r\n");
            biasAgeProErr = atof(tok);

            fgets(buffer,MAXBUF,fp6);
            tok = strtok(buffer," \t\r\n");
            NMortEst = atof(tok);

            for (j = 0; j < AGEPROCV; j++)
            {
                fgets(buffer,MAXBUF,fp6);
                tok = strtok(buffer," \t\r\n");
                for (i = 0; i < NAgePlus; i++)
                {
                    AgeProCV[i][j] = atof(tok);
                    tok = strtok(NULL," \t\r\n");
                }
            }

            fclose(fp6);
        }

        void WriteAgeProInputFile(char *fn)
        {
            char xname[FILBUF];
            char *c;
            long i, n;
            double xx, tl;
            double sd,zscore;

            strcpy(xname,fn);

            if ((c = strrchr(xname,'.')) != NULL)
                *c = '\0';

            strcat(xname,"_agepro.inp");

            if ((fp1 = fopen(xname,"w")) == NULL)
            {
                fprintf(stderr,"Unable to Open AgePro Input File: %s\n",xname);
                exit(1);
            }

            n = NFYear + NYears - 1;

            fprintf(fp1,"AGEPRO VERSION 4.0\n");
            fprintf(fp1,"[CASEID]\n");
            fprintf(fp1,"MSE Version 4.1\n");
            fprintf(fp1,"[GENERAL]\n");
            fprintf(fp1,"%ld  %ld  1  %ld  %ld  1  1  %ld %ld\n",n+1,n+2,NAgePlus,NAgeProSim,DiscFlag,ISeed);
            fprintf(fp1,"[BOOTSTRAP]\n");
            fprintf(fp1,"%ld  1\n",NBoot);

            strcpy(xname,fn);

            if ((c = strrchr(xname,'.')) != NULL)
                *c = '\0';

            if (ModelFlag == 1)
                strcat(xname,"_vpa.bsn");
            else if (ModelFlag == 2)
                strcat(xname,"_asap.bsn");

            fprintf(fp1,"%s\n",xname);

            fprintf(fp1,"[STOCK_WEIGHT]\n");
            fprintf(fp1,"0  0\n");

            for (i = 0; i < NAgePlus; i++)
                fprintf(fp1,"%10.6f ",StockWeight[i][NYears-1]);
            fprintf(fp1,"\n");
            for (i = 0; i < NAgePlus; i++)
                fprintf(fp1,"%10.6f ",AgeProCV[i][0]);
            fprintf(fp1,"\n");

            fprintf(fp1,"[CATCH_WEIGHT]\n");
            fprintf(fp1,"0  0\n");

            for (i = 0; i < NAgePlus; i++)
                fprintf(fp1,"%10.6f ",CatchWeight[i][NYears-1]);
            fprintf(fp1,"\n");
            for (i = 0; i < NAgePlus; i++)
                fprintf(fp1,"%10.6f ",AgeProCV[i][1]);
            fprintf(fp1,"\n");

            fprintf(fp1,"[MEAN_WEIGHT]\n");
            fprintf(fp1,"0  0\n");

            for (i = 0; i < NAgePlus; i++)
                fprintf(fp1,"%10.6f ",CatchWeight[i][NYears-1]);
            fprintf(fp1,"\n");
            for (i = 0; i < NAgePlus; i++)
                fprintf(fp1,"%10.6f ",AgeProCV[i][2]);
            fprintf(fp1,"\n");

            fprintf(fp1,"[SSB_WEIGHT]\n");
            fprintf(fp1,"0  0\n");

            for (i = 0; i < NAgePlus; i++)
                fprintf(fp1,"%10.6f ",SpawnWeight[i][NYears-1]);
            fprintf(fp1,"\n");
            for (i = 0; i < NAgePlus; i++)
                fprintf(fp1,"%10.6f ",AgeProCV[i][3]);
            fprintf(fp1,"\n");

            if (DiscFlag)
            {
                fprintf(fp1,"[DISC_WEIGHT]\n");
                fprintf(fp1,"0  0\n");

                for (i = 0; i < NAgePlus; i++)
                    fprintf(fp1,"%10.6f ",CatchWeight[i][NYears-1]);
                fprintf(fp1,"\n");
                for (i = 0; i < NAgePlus; i++)
                    fprintf(fp1,"%10.6f ",AgeProCV[i][4]);
                fprintf(fp1,"\n");
            }

            fprintf(fp1,"[NATMORT]\n");
            fprintf(fp1,"0  0\n");

            for (i = 0; i < NAgePlus; i++)
                fprintf(fp1,"%10.6f ",NMortEst);
            fprintf(fp1,"\n");
            for (i = 0; i < NAgePlus; i++)
                fprintf(fp1,"%10.6f ",AgeProCV[i][5]);
            fprintf(fp1,"\n");

            fprintf(fp1,"[MATURITY]\n");
            fprintf(fp1,"0  0\n");

            for (i = 0; i < NAgePlus; i++)
                fprintf(fp1,"%10.6f ",Mature[i][NYears-1]);
            fprintf(fp1,"\n");
            for (i = 0; i < NAgePlus; i++)
                fprintf(fp1,"%10.6f ",AgeProCV[i][6]);
            fprintf(fp1,"\n");

            fprintf(fp1,"[BIOLOGICAL]\n");
            fprintf(fp1,"0\n");

            for (i = 0; i < 2; i++)
                fprintf(fp1,"%10.6f ",TF[i]);
            fprintf(fp1,"\n");
            for (i = 0; i < 2; i++)
                fprintf(fp1,"%10.6f ",TM[i]);
            fprintf(fp1,"\n");

            fprintf(fp1,"[FISHERY]\n");
            fprintf(fp1,"0  0\n");

            for (i = 0; i < NAgePlus; i++)
                fprintf(fp1,"%10.6f ",FSelec[i][NYears]);
            fprintf(fp1,"\n");
            for (i = 0; i < NAgePlus; i++)
                fprintf(fp1,"%10.6f ",AgeProCV[i][7]);
            fprintf(fp1,"\n");

            if (DiscFlag)
            {
                fprintf(fp1,"[DISCARD]\n");
                fprintf(fp1,"0  0\n");

                for (i = 0; i < NAgePlus; i++)
                    fprintf(fp1,"%10.6f ",DiscFrac[i][NYears-1]);
                fprintf(fp1,"\n");
                for (i = 0; i < NAgePlus; i++)
                    fprintf(fp1,"%10.6f ",AgeProCV[i][8]);
                fprintf(fp1,"\n");
            }

            fprintf(fp1,"[RECRUIT]\n");
            fprintf(fp1,"1  1  %ld\n",MaxYears);
            fprintf(fp1,"14\n");
            fprintf(fp1,"1.0\n");
            fprintf(fp1,"1.0\n");
            fprintf(fp1,"%ld\n",NYears);

            if (ModelFlag == 1)
            {
                for (i = 0; i < NYears; i++)
                    fprintf(fp1,"%18.3f ",StockVPA[0][i]);
                fprintf(fp1,"\n");
            }
            else if (ModelFlag == 2)
            {
                for (i = 0; i < NYears; i++)
                    fprintf(fp1,"%18.3f ",StockAsap[i][0]);
                fprintf(fp1,"\n");
            }

            fprintf(fp1,"[HARVEST]\n");
            fprintf(fp1,"1  0\n");

            sd = sqrt(log(cvAgeProErr * cvAgeProErr + 1.0));
            zscore = gsl_box_muller();
            xx = exp(sd * zscore);
            tl = TargetLand * (1.0 + biasAgeProErr) * xx;

            fprintf(fp1,"%18.3f  %18.6f\n",tl,FTarg);

            fclose(fp1);

        }

        void ScanAgeProResults(char *fn)
        {
            char xname[FILBUF];
            char buffer[MAXBUF];
            char *c;
            char *t;
            short i;
            short flag;

            strcpy(xname,fn);

            if ((c = strrchr(xname,'.')) != NULL)
                *c = '\0';

            strcat(xname,"_agepro.out");

            if ((fp1 = fopen(xname,"r")) == NULL)
            {
                fprintf(stderr,"Unable to Open AgePro Output File: %s\n",xname);
                exit(1);
            }

            flag = 0;
            while(!feof(fp1))
            {
                fgets(buffer,MAXBUF,fp1);
                if (strstr(buffer,"Combined Catch Distribution"))
                {
                    flag = 1;
                    break;
                }
            }

            if (!flag)
            {
                fprintf(stderr,"Missing AgePro Catch\n");
                exit(1);
            }

            for (i = 0; i < 3; i++)
                fgets(buffer,MAXBUF,fp1);

            fgets(buffer,MAXBUF-1,fp1);
            t = strtok(buffer," \t\r\n");
            for (i = 0; i < 5; i++)
                t = strtok(NULL," \t\r\n");

            ProjLand = atof(t);

            ProjLand = ProjLand * 1000.0;

            fclose(fp1);
        }

        void RemoveOutputFiles(char *fn)
        {
            char xname[FILBUF];
            char *c;

            strcpy(xname,fn);

            if ((c = strrchr(xname,'.')) != NULL)
                *c = '\0';

            strcat(xname,"_vpa.out");

            remove(xname);

            if ((c = strrchr(xname,'.')) != NULL)
                *c = '\0';

            strcat(xname,".aux");

            remove(xname);

            *c = '\0';

            strcat(xname,".pp2");

            remove(xname);

            *c = '\0';

            strcat(xname,".bsn");

            remove(xname);

            *c = '\0';

            strcat(xname,".bsf");

            remove(xname);

            *c = '\0';

            strcat(xname,".bsq");

            remove(xname);

            *c = '\0';

            strcat(xname,".bsb");

            remove(xname);

            *c = '\0';

            strcat(xname,".rsd");

            remove(xname);

            *c = '\0';

            strcat(xname,".bxx");

            remove(xname);

            *c = '\0';

            strcat(xname,".cv");

            remove(xname);

            *c = '\0';

            strcat(xname,".log");

            remove(xname);

            strcpy(xname,fn);

            if ((c = strrchr(xname,'.')) != NULL)
                *c = '\0';

            strcat(xname,"_apro.in");

            remove(xname);

            if ((c = strrchr(xname,'.')) != NULL)
                *c = '\0';

            strcat(xname,".out");

            remove(xname);

            /* Aspic */

            strcpy(xname,fn);

            if ((c = strrchr(xname,'.')) != NULL)
                *c = '\0';

            strcat(xname,"_aspic.bot");

            remove(xname);


            if ((c = strrchr(xname,'.')) != NULL)
                *c = '\0';

            strcat(xname,".bio");

            remove(xname);

            if ((c = strrchr(xname,'.')) != NULL)
                *c = '\0';

            strcat(xname,".fit");

            remove(xname);

            if ((c = strrchr(xname,'.')) != NULL)
                *c = '\0';

            strcat(xname,".ctl");

            remove(xname);


            if ((c = strrchr(xname,'.')) != NULL)
                *c = '\0';

            strcat(xname,".prb");

            remove(xname);
            if ((c = strrchr(xname,'.')) != NULL)
                *c = '\0';

            strcat(xname,".prj");

            remove(xname);


            if (ModelFlag == 2)
            {
                strcpy(xname,fn);

                if ((c = strrchr(xname,'\\')) != NULL)
                {
                    *c = '\0';
                    strcat(xname,"\\asap3.dat");
                }
                else
                    strcpy(xname,"asap3.dat");

                remove(xname);

                if ((c = strrchr(xname,'.')) != NULL)
                    *c = '\0';

                strcat(xname,".rep");
                remove(xname);

                if ((c = strrchr(xname,'.')) != NULL)
                    *c = '\0';

                strcat(xname,".bsn");
                remove(xname);


                if ((c = strrchr(xname,'.')) != NULL)
                    *c = '\0';

                strcat(xname,".std");
                remove(xname);

                if ((c = strrchr(xname,'.')) != NULL)
                    *c = '\0';

                strcat(xname,".mcm");
                remove(xname);

                if ((c = strrchr(xname,'\\')) != NULL)
                {
                    *c = '\0';
                    strcat(xname,"\\asap3mcmc.dat");
                }
                else
                    strcpy(xname,"asap3mcmc.dat");

                remove(xname);

                strcpy(xname,fn);

                if ((c = strrchr(xname,'.')) != NULL)
                    *c = '\0';

                strcat(xname,"_asap.dat");

                remove(xname);

                if ((c = strrchr(xname,'.')) != NULL)
                    *c = '\0';

                strcat(xname,".rep");
                remove(xname);

                if ((c = strrchr(xname,'.')) != NULL)
                    *c = '\0';

                strcat(xname,".bsn");
                remove(xname);
            }
        }

        void ScanVPAResults(char *fn)
        {
            char buffer[MAXBUF];
            char xname[FILBUF];
            char *tok, *c;
            short i, j;

            /* Open VPA Auxilliary Output File */

            strcpy(xname,fn);

            if ((c = strrchr(xname,'.')) != NULL)
                *c = '\0';

            strcat(xname,"_vpa.pp2");

            if ((fp1 = fopen(xname,"r")) == NULL)
            {
                fprintf(stderr,"Unable to Open VPA Output File: %s\n",xname);
                exit(1);
            }

            for (i = 0; i < KAges; i++)
            {
                fgets(buffer,MAXBUF-1,fp1);
                tok = strtok(buffer," \t\r\n");
                for (j = 0; j < NYears; j++)
                {
                    StockVPA[i][j] = atof(tok);
                    tok = strtok(NULL," \t\r\n");
                }
            }

            fclose(fp1);

            if ((c = strrchr(xname,'.')) != NULL)
                *c = '\0';

            strcat(xname,".aux");

            if ((fp1 = fopen(xname,"r")) == NULL)
            {
                fprintf(stderr,"Unable to Open VPA Output File: %s\n",xname);
                exit(1);
            }

            fgets(buffer,MAXBUF-1,fp1);
            for (i = 0; i < KAges; i++)
            {
                fgets(buffer,MAXBUF-1,fp1);
                tok = strtok(buffer," \t\r\n");
                for (j = 0; j < NYears; j++)
                {
                    FVPA[i][j] = atof(tok);
                    tok = strtok(NULL," \t\r\n");
                }
            }

            /* Skip JAN-1 Biomass */

            for (i = 0; i < KAges+2; i++)
                fgets(buffer,MAXBUF-1,fp1);

            for (i = 0; i < KAges; i++)
            {
                fgets(buffer,MAXBUF-1,fp1);
                tok = strtok(buffer," \t\r\n");
                for (j = 0; j < NYears; j++)
                {
                    SSBVPA[i][j] = atof(tok);
                    tok = strtok(NULL," \t\r\n");
                }
            }

            fclose(fp1);

        }

        long LaunchVPA(char *fn)
        {
            char xname[FILBUF];
            char cmdline[MAXBUF];
            char *c;
            long k;
            char qt[2];

            qt[0] = '\x22';
            qt[1] = '\0';

            strcpy(xname,fn);

            if ((c = strrchr(xname,'.')) != NULL)
                *c = '\0';

            strcat(xname,"_vpa.dat");


            strcpy(cmdline,ModelPath);
            strcat(cmdline," ");
            strcat(cmdline,qt);
            strcat(cmdline,xname);
            strcat(cmdline,qt);
            LowerCase(cmdline);

            k = CreateConsoleProcess(cmdline);

            return k;
        }

        long LaunchAgePro(char *fn)
        {
            char xname[FILBUF];
            char cmdline[MAXBUF];
            char *c;
            long k;
            char qt[2];

            qt[0] = '\x22';
            qt[1] = '\0';

            strcpy(xname,fn);

            if ((c = strrchr(xname,'.')) != NULL)
                *c = '\0';

            strcat(xname,"_agepro.inp");

            strcpy(cmdline,AgeProPath);
            strcat(cmdline," ");
            strcat(cmdline,qt);
            strcat(cmdline,xname);
            strcat(cmdline,qt);
            LowerCase(cmdline);

            k = CreateConsoleProcess(cmdline);

            return k;
        }

        void ReadAsapTemplateFile(char *fn)
        {
            char buffer[MAXBUF];
            char *c, *tok;
            long i, j, k;

            if ((c = strrchr(fn,'.')) != NULL)
                *c = '\0';

            strcat(fn,".tx3");

            if ((fp1 = fopen(fn,"r")) == NULL)
            {
                fprintf(stderr,"Unable to Open Template File: %s\n",fn);
                exit(1);
            }


            Fleetblock   = (short *) calloc(MaxYears+1,sizeof(short));


            Fleetsamp = (double *) calloc(MaxYears+1,sizeof(double));
            Fleetcv   = (double *) calloc(MaxYears+1,sizeof(double));
            Discsamp  = (double *) calloc(MaxYears+1,sizeof(double));
            Disccv    = (double *) calloc(MaxYears+1,sizeof(double));


            Recruitcv = (double *) calloc(MaxYears+1,sizeof(double));


            InitNAA     = (double *) calloc(MaxAge,sizeof(double));


            StockAsap   = AllocMatrix(MaxYears+1,KAges);

            NMortMat  = AllocMatrix(MaxYears+1,KAges);

            while(!feof(fp1))
            {
                fgets(buffer,MAXBUF,fp1);
                if (strstr(buffer,"PATH"))
                {
                    fgets(buffer,MAXBUF,fp1);
                    if ((c = strchr(buffer,'\n')) != NULL)
                        *c = '\0';
                    strcpy(ModelPath,buffer);
                }
                else if (strstr(buffer,"INDEX"))
                {
                    fgets(buffer,MAXBUF,fp1);
                    tok = strtok(buffer," \t\r\n");
                    NIndex = atol(tok);

                    IndexValues = AllocMatrix(MaxYears,NIndex);

                    IndexData = (struct si *) calloc(NIndex,(sizeof(struct si)));
                    for (i = 0; i < NIndex; i++)
                    {
                        fgets(buffer,MAXBUF,fp1);
                        tok = strtok(buffer," \t\r\n");
                        IndexData[i].StartAge = atoi(tok);
                        tok = strtok(NULL," \t\r\n");
                        IndexData[i].EndAge = atoi(tok);
                        tok = strtok(NULL," \t\r\n");
                        IndexData[i].Type = atoi(tok);
                        tok = strtok(NULL," \t\r\n");
                        IndexData[i].SurveyIndex = atoi(tok);
                    }

                    IndexSelType = (short *) calloc(NIndex,sizeof(short));
                    LambdaIndex = (double *) calloc(NIndex,sizeof(double));
                    Lambdaq     = (double *) calloc(NIndex,sizeof(double));
                    Lambdaqdev  = (double *) calloc(NIndex,sizeof(double));
                    cvq         = (double *) calloc(NIndex,sizeof(double));
                    cvqdev      = (double *) calloc(NIndex,sizeof(double));
                    Initq       = (double *) calloc(NIndex,sizeof(double));

                    k = (NSELMOD + KAges);
                    k = k * NIndex;
                    Indexsel  = (struct seldat *) calloc(k,sizeof(struct seldat));

                }
                else if (strstr(buffer,"NATMORT"))
                {
                    for (i = 0; i < MaxYears; i++)
                    {
                        fgets(buffer,MAXBUF,fp1);
                        tok = strtok(buffer," \t\r\n");
                        for (j = 0; j < KAges; j++)
                        {
                            NMortMat[i][j] = atof(tok);
                            tok = strtok(NULL," \t\r\n");
                        }
                    }
                }
                else if (strstr(buffer,"BLOCK"))
                {
                    fgets(buffer,MAXBUF,fp1);
                    NBlock = atoi(buffer);
                    FleetSelType = (short *) calloc(NBlock,sizeof(short));
                    k = NSELMOD + KAges;
                    k = k * NBlock;
                    Fleetsel  = (struct seldat *) calloc(k,sizeof(struct seldat));

                    fgets(buffer,MAXBUF,fp1);
                    tok = strtok(buffer," \t\r\n");
                    for (j = 0; j < MaxYears; j++)
                    {
                        Fleetblock[j] = atoi(tok);
                        tok = strtok(NULL," \t\r\n");
                    }

                    for (i = 0; i < NBlock; i++)
                    {
                        fgets(buffer,MAXBUF,fp1);
                        FleetSelType[i] = atoi(buffer);
                        for (j = 0; j < NSELMOD + KAges; j++)
                        {
                            k = i * (NSELMOD + KAges) + j;
                            fgets(buffer,MAXBUF,fp1);
                            tok = strtok(buffer," \t\r\n");
                            Fleetsel[k].value = atof(tok);
                            tok = strtok(NULL," \t\r\n");
                            Fleetsel[k].phase = atoi(tok);
                            tok = strtok(NULL," \t\r\n");
                            Fleetsel[k].lambda = atof(tok);
                            tok = strtok(NULL," \t\r\n");
                            Fleetsel[k].cv = atof(tok);
                        }
                    }
                }
                else if (strstr(buffer,"AVGF"))
                {
                    fgets(buffer,MAXBUF,fp1);
                    tok = strtok(buffer," \t\r\n");
                    AvgF[0] = atoi(tok);
                    tok = strtok(NULL," \t\r\n");
                    AvgF[1] = atoi(tok);
                }
                else if (strstr(buffer,"INDX_SELEC"))
                {
                    fgets(buffer,MAXBUF,fp1);
                    tok = strtok(buffer," \t\r\n");
                    IndexDatcv = atof(tok);
                    for (i = 0; i < NIndex; i++)
                    {
                        fgets(buffer,MAXBUF,fp1);
                        IndexSelType[i] = atoi(buffer);
                        for (j = 0; j < NSELMOD + KAges; j++)
                        {
                            k = i * (NSELMOD + KAges) + j;
                            fgets(buffer,MAXBUF,fp1);
                            tok = strtok(buffer," \t\r\n");
                            Indexsel[k].value = atof(tok);
                            tok = strtok(NULL," \t\r\n");
                            Indexsel[k].phase = atoi(tok);
                            tok = strtok(NULL," \t\r\n");
                            Indexsel[k].lambda = atof(tok);
                            tok = strtok(NULL," \t\r\n");
                            Indexsel[k].cv = atof(tok);
                        }
                    }
                }
                else if (strstr(buffer,"CATCH"))
                {

                    for (i = 0; i < MaxYears; i++)
                    {
                        fgets(buffer,MAXBUF,fp1);
                        tok = strtok(buffer," \t\r\n");
                        Fleetsamp[i] = atof(tok);
                        tok = strtok(NULL," \t\r\n");
                        Fleetcv[i] = atof(tok);
                    }
                }
                else if (strstr(buffer,"DISCARD"))
                {
                    for (i = 0; i < MaxYears; i++)
                    {
                        fgets(buffer,MAXBUF,fp1);
                        tok = strtok(buffer," \t\r\n");
                        Discsamp[i] = atof(tok);
                        tok = strtok(NULL," \t\r\n");
                        Disccv[i] = atof(tok);
                    }
                }
                else if (strstr(buffer,"RECRUIT"))
                {
                    for (i = 0; i < MaxYears; i++)
                    {
                        fgets(buffer,MAXBUF,fp1);
                        tok = strtok(buffer," \t\r\n");
                        Recruitcv[i] = atof(tok);
                    }
                }
                else if (strstr(buffer,"PHASE"))
                {
                    for (i = 0; i < NPHASE; i++)
                    {
                        fgets(buffer,MAXBUF,fp1);
                        tok = strtok(buffer," \t\r\n");
                        PhaseAsap[i] = atoi(tok);
                    }
                }
                else if (strstr(buffer,"LAMBDA"))
                {
                    fgets(buffer,MAXBUF,fp1);
                    tok = strtok(buffer," \t\r\n");
                    for (j = 0; j < NIndex; j++)
                    {
                        LambdaIndex[j] = atof(tok);
                        tok = strtok(NULL," \t\r\n");
                    }

                    fgets(buffer,MAXBUF,fp1);
                    tok = strtok(buffer," \t\r\n");
                    LambdaFleet = atof(tok);

                    fgets(buffer,MAXBUF,fp1);
                    tok = strtok(buffer," \t\r\n");
                    LambdaDisc = atof(tok);

                    fgets(buffer,MAXBUF,fp1);
                    tok = strtok(buffer," \t\r\n");
                    LambdaFMult = atof(tok);

                    fgets(buffer,MAXBUF,fp1);
                    tok = strtok(buffer," \t\r\n");
                    cvFMult = atof(tok);

                    fgets(buffer,MAXBUF,fp1);
                    tok = strtok(buffer," \t\r\n");
                    LambdaFDev = atof(tok);

                    fgets(buffer,MAXBUF,fp1);
                    tok = strtok(buffer," \t\r\n");
                    cvFDev = atof(tok);

                    fgets(buffer,MAXBUF,fp1);
                    tok = strtok(buffer," \t\r\n");
                    LambdaN1 = atof(tok);
                    tok = strtok(NULL," \t\r\n");
                    cvN1 = atof(tok);

                    fgets(buffer,MAXBUF,fp1);
                    tok = strtok(buffer," \t\r\n");
                    LambdaRec = atof(tok);

                    fgets(buffer,MAXBUF,fp1);
                    tok = strtok(buffer," \t\r\n");
                    for (j = 0; j < NIndex; j++)
                    {
                        Lambdaq[j] = atof(tok);
                        tok = strtok(NULL," \t\r\n");
                    }

                    fgets(buffer,MAXBUF,fp1);
                    tok = strtok(buffer," \t\r\n");
                    for (j = 0; j < NIndex; j++)
                    {
                        cvq[j] = atof(tok);
                        tok = strtok(NULL," \t\r\n");
                    }

                    fgets(buffer,MAXBUF,fp1);
                    tok = strtok(buffer," \t\r\n");
                    for (j = 0; j < NIndex; j++)
                    {
                        Lambdaqdev[j] = atof(tok);
                        tok = strtok(NULL," \t\r\n");
                    }

                    fgets(buffer,MAXBUF,fp1);
                    tok = strtok(buffer," \t\r\n");
                    for (j = 0; j < NIndex; j++)
                    {
                        cvqdev[j] = atof(tok);
                        tok = strtok(NULL," \t\r\n");
                    }

                    fgets(buffer,MAXBUF,fp1);
                    tok = strtok(buffer," \t\r\n");
                    LambdaSteep = atof(tok);
                    tok = strtok(NULL," \t\r\n");
                    cvSteep = atof(tok);

                    fgets(buffer,MAXBUF,fp1);
                    tok = strtok(buffer," \t\r\n");
                    LambdaVirgin = atof(tok);
                    tok = strtok(NULL," \t\r\n");
                    cvVirgin = atof(tok);
                }
                else if (strstr(buffer,"INITIAL"))
                {
                    fgets(buffer,MAXBUF,fp1);
                    tok = strtok(buffer," \t\r\n");
                    NAAFlag = atoi(tok);

                    fgets(buffer,MAXBUF,fp1);
                    tok = strtok(buffer," \t\r\n");
                    for (j = 0; j < KAges; j++)
                    {
                        InitNAA[j] = atof(tok);
                        tok = strtok(NULL," \t\r\n");
                    }

                    fgets(buffer,MAXBUF,fp1);
                    tok = strtok(buffer," \t\r\n");
                    InitFMult = atof(tok);

                    fgets(buffer,MAXBUF,fp1);
                    tok = strtok(buffer," \t\r\n");
                    for (j = 0; j < NIndex; j++)
                    {
                        Initq[j] = atof(tok);
                        tok = strtok(NULL," \t\r\n");
                    }

                    fgets(buffer,MAXBUF,fp1);
                    tok = strtok(buffer," \t\r\n");
                    SRRFlag = atoi(tok);

                    fgets(buffer,MAXBUF,fp1);
                    tok = strtok(buffer," \t\r\n");
                    InitVirgin = atof(tok);

                    fgets(buffer,MAXBUF,fp1);
                    tok = strtok(buffer," \t\r\n");
                    InitSteep = atof(tok);

                    fgets(buffer,MAXBUF,fp1);
                    tok = strtok(buffer," \t\r\n");
                    InitFMax = atof(tok);
                }
                else if (strstr(buffer,"MCMC"))
                {
                    fgets(buffer,MAXBUF,fp1);
                    tok = strtok(buffer," \t\r\n");
                    NBoot = atol(tok);
                    tok = strtok(NULL," \t\r\n");
                    KThin = atol(tok);

                    fgets(buffer,MAXBUF,fp1);
                    tok = strtok(buffer," \t\r\n");
                    KStart = atol(tok);
                    tok = strtok(NULL," \t\r\n");
                    KEnd = atol(tok);
                }
            }

            fclose(fp1);

            /* Expand Data */

            Fleetcv[MaxYears]    = Fleetcv[MaxYears-1];
            Fleetsamp[MaxYears]  = Fleetsamp[MaxYears-1];
            Disccv[MaxYears]     = Disccv[MaxYears-1];
            Discsamp[MaxYears]   = Discsamp[MaxYears-1];
            Fleetcv[MaxYears]    = Fleetcv[MaxYears-1];
            Fleetsamp[MaxYears]  = Fleetsamp[MaxYears-1];
            Recruitcv[MaxYears]  = Recruitcv[MaxYears-1];
            Fleetblock[MaxYears] = Fleetblock[MaxYears-1];

        }

        void WriteAsapInputFile(char *fn)
        {
            char xname[FILBUF];
            char *c;
            long i, j, k;
            long k1,k2;
            double fd;
            double xx;

            strcpy(xname,fn);

            if ((c = strrchr(xname,'\\')) != NULL)
            {
                *c = '\0';
                strcat(xname,"\\asap3.dat");
            }
            else
                strcpy(xname,"asap3.dat");


            if ((fp1 = fopen(xname,"w")) == NULL)
            {
                fprintf(stderr,"Unable to Open Asap Input File: %s\n",xname);
                exit(1);
            }

            fprintf(fp1,"# ASAP VERSION 3.0\n");
            fprintf(fp1,"#MSE Version 4.1\n");
            fprintf(fp1,"#\n");
            fprintf(fp1,"# Number of Years\n");
            fprintf(fp1,"%d\n",NYears);
            fprintf(fp1,"# First Year\n");
            fprintf(fp1,"%d\n",NFYear);
            fprintf(fp1,"# Number of Ages\n");
            fprintf(fp1,"%d\n",KAges);
            fprintf(fp1,"# Number of Fleets\n");
            fprintf(fp1,"1\n");
            fprintf(fp1,"# Number of Selectivity Blocks\n");
            fprintf(fp1,"%d\n",NBlock);
            fprintf(fp1,"# Number of Indices\n");
            fprintf(fp1,"%d\n",NIndex);

            fprintf(fp1,"# Natural Mortality\n");
            for (i = 0; i < NYears; i++)
            {
                for (j = 0; j < KAges; j++)
                    fprintf(fp1,"%8.4f  ",NMortMat[i][j]);
                fprintf(fp1,"\n");
            }
            fprintf(fp1,"# Fecundity Option\n0\n");
            fprintf(fp1,"# Fraction of year that elapses prior to SSB calculation (0=Jan-1)\n");
            fprintf(fp1,"%8.4f\n",TF[NYears-1]);
            fprintf(fp1,"# Maturity\n");
            for (i = 0; i < NYears; i++)
            {
                for (j = 0; j < KAges; j++)
                    fprintf(fp1,"%-12.4f ",Mature[j][i]);
                fprintf(fp1,"\n");
            }

            fprintf(fp1,"# Number of Weight at Age Matrices\n");
            fprintf(fp1,"3\n");

            fprintf(fp1,"# Weight at Age for Catch Matrix\n");
            for (i = 0; i < NYears; i++)
            {
                for (j = 0; j < KAges; j++)
                    fprintf(fp1,"%-12.4f ",CatchWeight[j][i]);
                fprintf(fp1,"\n");
            }

            fprintf(fp1,"# Weight at Age for Spawning Stock\n");
            for (i = 0; i < NYears; i++)
            {
                for (j = 0; j < KAges; j++)
                    fprintf(fp1,"%-12.4f ",SpawnWeight[j][i]);
                fprintf(fp1,"\n");
            }

            fprintf(fp1,"# Weight at Age for Jan-1 Stock\n");
            for (i = 0; i < NYears; i++)
            {
                for (j = 0; j < KAges; j++)
                    fprintf(fp1,"%-12.4f ",StockWeight[j][i]);
                fprintf(fp1,"\n");
            }

            fprintf(fp1,"# Weights at Age Pointers\n");
            k = 4;

            for (i = 0; i < k; i++)
                fprintf(fp1,"1\n");

            fprintf(fp1,"2\n");
            fprintf(fp1,"3\n");


            fprintf(fp1,"# Selectivity Blocks\n");

            fprintf(fp1,"# Selectivity Blocks for Fleet-%d\n",k+1);
            for (i = 0; i < NYears; i++)
                fprintf(fp1,"%d\n",Fleetblock[i]);


            fprintf(fp1,"#  Selectivity Options for each block 1=by age, 2=logisitic, 3=double logistic\n");
            for (i = 0; i < NBlock; i++)
                fprintf(fp1,"%d  ",FleetSelType[i]);
            fprintf(fp1,"\n");

            fprintf(fp1,"# Selectivity Input by Block\n");
            for (i = 0; i < NBlock; i++)
            {
                fprintf(fp1,"# Block-%d\n",i+1);
                for (j = 0; j < NSELMOD + KAges; j++)
                {
                    k = i * (NSELMOD + KAges) + j;
                    fprintf(fp1,"%10.4f  ",Fleetsel[k].value);
                    fprintf(fp1,"%4d  ",Fleetsel[k].phase);
                    fprintf(fp1,"%10.4f  ",Fleetsel[k].lambda);
                    fprintf(fp1,"%10.4f\n",Fleetsel[k].cv);
                }
            }
            fprintf(fp1,"#  Fleet Start Age\n");
            fprintf(fp1,"1  ");
            fprintf(fp1,"\n");
            fprintf(fp1,"#  Fleet End Age\n");
            fprintf(fp1,"%d  ",KAges);
            fprintf(fp1,"\n");
            fprintf(fp1,"# Ages for Average F\n");
            fprintf(fp1,"%d   %d\n",AvgF[0],AvgF[1]);
            fprintf(fp1,"# Average F Option\n");
            fprintf(fp1,"1\n");
            fprintf(fp1,"# Use Likelihood\n");
            fprintf(fp1,"1\n");

            fprintf(fp1,"# Release Mortality by Fleet\n");
            fprintf(fp1,"1.0  ");
            fprintf(fp1,"\n");

            fprintf(fp1,"# Fleet-1 Catch at Age\n");
            for (j = 0; j < NYears; j++)
            {

                xx = 0.0;
                for (i = 0; i < KAges; i++)
                {
                    fprintf(fp1,"%-16.6E ",CatchSamples[i][j]);
                    xx += CatchWeight[i][j] * CatchSamples[i][j];
                }
                fprintf(fp1,"%16.6E\n",xx);
            }

            fprintf(fp1,"# Fleet-1 Discards at Age\n");
            for (j = 0; j < NYears; j++)
            {

                xx = 0.0;
                for (i = 0; i < KAges; i++)
                {
                    fprintf(fp1,"%-16.6E ",DiscFrac[i][j] * CatchSamples[i][j]);
                    xx += DiscFrac[i][j] * CatchWeight[i][j] * CatchSamples[i][j];
                }
                fprintf(fp1,"%16.6E\n",xx);
            }

            fprintf(fp1,"# Fleet-1 Release Proportion at Age\n");
            for (j = 0; j < NYears; j++)
            {
                for (i = 0; i < KAges; i++)
                {
                    fd = 0.0;
                    fprintf(fp1,"%-10.6f ",fd);
                }
                fprintf(fp1,"\n");
            }

            fprintf(fp1,"# Aggregate Index Units\n");
            for (i = 0; i < NIndex; i++)
            {
                if (IndexData[i].Type == 0)
                    fprintf(fp1,"2  ");
                else
                    fprintf(fp1,"1  ");
            }
            fprintf(fp1,"\n");

            fprintf(fp1,"# Age Proportion Units\n");
            for (i = 0; i < NIndex; i++)
            {
                if (IndexData[i].Type == 0)
                    fprintf(fp1,"2  ");
                else
                    fprintf(fp1,"1  ");
            }
            fprintf(fp1,"\n");

            fprintf(fp1,"# Weight at Age Matrix\n");
            for (i = 0; i < NIndex; i++)
                fprintf(fp1,"1  ");
            fprintf(fp1,"\n");

            fprintf(fp1,"# Index Month\n");
            for (i = 0; i < NIndex; i++)
                fprintf(fp1,"1  ");
            fprintf(fp1,"\n");

            fprintf(fp1,"# Index Selectivity Link to Fleet\n");
            for (i = 0; i < NIndex; i++)
                fprintf(fp1,"-1  ");
            fprintf(fp1,"\n");

            fprintf(fp1,"# Index Selectivity Option\n");
            for (i = 0; i < NIndex; i++)
                fprintf(fp1,"%d  ",IndexSelType[i]);
            fprintf(fp1,"\n");

            fprintf(fp1,"# Index Start Age\n");
            for (i = 0; i < NIndex; i++)
                fprintf(fp1,"%d  ",IndexData[i].StartAge);
            fprintf(fp1,"\n");

            fprintf(fp1,"# Index End Age\n");
            for (i = 0; i < NIndex; i++)
                if (IndexData[i].EndAge <= KAges)
                    fprintf(fp1,"%d  ",IndexData[i].EndAge);
                else
                    fprintf(fp1,"%d  ",KAges);

            fprintf(fp1,"\n");

            fprintf(fp1,"# Estimate Proportion\n");
            for (i = 0; i < NIndex; i++)
                fprintf(fp1,"0  ");
            fprintf(fp1,"\n");

            fprintf(fp1,"# Use Index\n");
            for (i = 0; i < NIndex; i++)
                fprintf(fp1,"1  ");
            fprintf(fp1,"\n");

            for (i = 0; i < NIndex; i++)
            {
                fprintf(fp1,"# INDEX-%d  Selectivity\n",i+1);
                for (j = 0; j < NSELMOD + KAges; j++)
                {
                    k = i * (NSELMOD + KAges) + j;
                    fprintf(fp1,"%10.4f  ",Indexsel[k].value);
                    fprintf(fp1,"%4d  ",Indexsel[k].phase);
                    fprintf(fp1,"%10.4f  ",Indexsel[k].lambda);
                    fprintf(fp1,"%10.4f\n",Indexsel[k].cv);
                }
            }

            fprintf(fp1,"# Index Data\n");
            for (k = 0; k < NIndex; k++)
            {
                fprintf(fp1,"# INDEX-%d\n",k+1);
                for (j = 0; j < NYears; j++)
                {

                    xx = 0.0;

                    k1 = IndexData[k].StartAge -1;
                    k2 = IndexData[k].EndAge - 1;

                    for (i = k1; i <= k2; i++)
                        xx += SurveySamples[i][j][k];


                    if (xx < XTOL)
                        xx = -999.0;

                    fprintf(fp1,"%4ld   %18.6E  %8.4f   ",NFYear+j,xx,IndexDatcv);
                    for (i = 0; i < KAges; i++)
                        fprintf(fp1,"0.0000   ");
                    fprintf(fp1,"0\n");
                }
            }

            fprintf(fp1,"# Phase Control\n");
            for (i = 0; i < NPHASE; i++)
            {
                fprintf(fp1,"#\n");
                fprintf(fp1,"%d\n",PhaseAsap[i]);
            }

            fprintf(fp1,"# Recruit CV\n");
            for (i = 0; i < NYears; i++)
                fprintf(fp1,"%10.4f\n",Recruitcv[i]);


            fprintf(fp1,"# Lambda Index\n");
            for (i = 0; i < NIndex; i++)
                fprintf(fp1,"%10.2f  ",LambdaIndex[i]);
            fprintf(fp1,"\n");

            fprintf(fp1,"# Lambda Total Catch\n");
            fprintf(fp1,"%10.2f  ",LambdaFleet);
            fprintf(fp1,"\n");

            fprintf(fp1,"# Lambda Total Discards\n");
            fprintf(fp1,"%10.2f  ",LambdaDisc);
            fprintf(fp1,"\n");

            fprintf(fp1,"# Catch cv\n");
            for (i = 0; i < NYears; i++)
            {
                fprintf(fp1,"%10.4f  ",Fleetcv[i]);
                fprintf(fp1,"\n");
            }

            fprintf(fp1,"# Discard cv\n");
            for (i = 0; i < NYears; i++)
            {
                fprintf(fp1,"%10.4f  ",Disccv[i]);
                fprintf(fp1,"\n");
            }


            fprintf(fp1,"# Catch Input Effective Saple Size\n");
            for (i = 0; i < NYears; i++)
            {
                fprintf(fp1,"%10.0f  ",Fleetsamp[i]);
                fprintf(fp1,"\n");
            }

            fprintf(fp1,"# Discard Input Effective Sample Size\n");
            for (i = 0; i < NYears; i++)
            {
                fprintf(fp1,"%10.0f  ",Discsamp[i]);
                fprintf(fp1,"\n");
            }

            fprintf(fp1,"# Lambda F Mult\n");
            fprintf(fp1,"%10.2f  ",LambdaFMult);
            fprintf(fp1,"\n");

            fprintf(fp1,"# CV F Mult\n");
            fprintf(fp1,"%10.2f  ",cvFMult);
            fprintf(fp1,"\n");

            fprintf(fp1,"# Lambda F Dev\n");
            fprintf(fp1,"%10.2f  ",LambdaFDev);
            fprintf(fp1,"\n");


            fprintf(fp1,"# CV F Dev\n");
                fprintf(fp1,"%10.2f  ",cvFDev);
            fprintf(fp1,"\n");

            fprintf(fp1,"# Lambda N 1st Year\n");
            fprintf(fp1,"%10.2f\n",LambdaN1);
            fprintf(fp1,"# CV N 1st Year\n");
            fprintf(fp1,"%10.4f\n",cvN1);

            fprintf(fp1,"# Lambda Rec Dev\n");
            fprintf(fp1,"%10.2f\n",LambdaRec);

            fprintf(fp1,"# Lambda Catchability\n");
            for (i = 0; i < NIndex; i++)
                fprintf(fp1,"%10.2f  ",Lambdaq[i]);
            fprintf(fp1,"\n");

            fprintf(fp1,"# CV Catchability\n");
            for (i = 0; i < NIndex; i++)
                fprintf(fp1,"%10.4f  ",cvq[i]);
            fprintf(fp1,"\n");

            fprintf(fp1,"# Lambda Catchability Devs\n");
            for (i = 0; i < NIndex; i++)
                fprintf(fp1,"%10.2f  ",Lambdaqdev[i]);
            fprintf(fp1,"\n");

            fprintf(fp1,"# CV Catchability Devs\n");
            for (i = 0; i < NIndex; i++)
                fprintf(fp1,"%10.4f  ",cvqdev[i]);
            fprintf(fp1,"\n");

            fprintf(fp1,"# Lambda Dev Initial Steepness\n");
            fprintf(fp1,"%10.2f\n",LambdaSteep);
            fprintf(fp1,"# CV Dev Initial Steepness\n");
            fprintf(fp1,"%10.4f\n",cvSteep);

            fprintf(fp1,"# Lambda Dev from Unexploited Biomass\n");
            fprintf(fp1,"%10.2f\n",LambdaVirgin);
            fprintf(fp1,"# CV Dev from Unexploited Biomass\n");
            fprintf(fp1,"%10.4f\n",cvVirgin);

            fprintf(fp1,"# Initial Guesses\n");

            fprintf(fp1,"# Numbers at Age 1st Year Option\n");
            fprintf(fp1,"%d\n",NAAFlag);

            fprintf(fp1,"# Init Numbers at Age 1st Year\n");
            for (i = 0; i < KAges; i++)
                fprintf(fp1,"%12.1f  ",InitNAA[i]);
            fprintf(fp1,"\n");

            fprintf(fp1,"# Init FMult\n");
            fprintf(fp1,"%10.5f  ",InitFMult);
            fprintf(fp1,"\n");

            fprintf(fp1,"# Init Catchability\n");
            for (i = 0; i < NIndex; i++)
                fprintf(fp1,"%10.5f  ",Initq[i]);
            fprintf(fp1,"\n");

            fprintf(fp1,"# Stock Recruitment Option\n");
            fprintf(fp1,"%d\n",SRRFlag);


            fprintf(fp1,"# Init Unexploited Biomass\n");
            fprintf(fp1,"%10.5f\n",InitVirgin);

            fprintf(fp1,"# Init Steepness\n");
            fprintf(fp1,"%10.5f\n",InitSteep);

            fprintf(fp1,"# Max F\n");
            fprintf(fp1,"%10.5f\n",InitFMax);

            fprintf(fp1,"# Ignore Guesses\n");
            fprintf(fp1,"0\n");

            fprintf(fp1,"# Projection Spec\n");
            fprintf(fp1,"# Do Projection\n");
            fprintf(fp1,"0\n");

            fprintf(fp1,"# Fleet Directed Flag\n");
            fprintf(fp1,"1  ");
            fprintf(fp1,"\n");

            fprintf(fp1,"# Final Year\n");
            fprintf(fp1,"%d\n",NFYear+NYears);

            fprintf(fp1,"#\n");
            fprintf(fp1,"%d  -1  3  -99  1\n",NFYear+NYears);

            fprintf(fp1,"# MCMC Spec\n");
            fprintf(fp1,"# Do MCMC\n");
            fprintf(fp1,"1\n");

            fprintf(fp1,"#MCMC Year Option\n");
            fprintf(fp1,"1\n");

            fprintf(fp1,"#MCMC Iterations\n");
            fprintf(fp1,"%ld\n",NBoot);
            fprintf(fp1,"#MCMC Thinning Rate\n");
            fprintf(fp1,"%ld\n",KThin);
            fprintf(fp1,"#MCMC Random Seed\n");
            fprintf(fp1,"%ld\n",ISeed);

            fprintf(fp1,"#AgePro R Option\n");
            fprintf(fp1,"0\n");
            fprintf(fp1,"#AgePro Start Year\n");
            fprintf(fp1,"%ld\n",KStart);
            fprintf(fp1,"#AgePro End Year\n");
            fprintf(fp1,"%ld\n",KEnd);

            fprintf(fp1,"# Export R Flag\n");
            fprintf(fp1,"0\n");

            fprintf(fp1,"#Test Value\n");
            fprintf(fp1,"-23456\n");
            fprintf(fp1,"# ---- FINIS ---\n");

            fprintf(fp1,"# Fleet Names\n");
            fprintf(fp1,"#$FLEET-1\n");

            fprintf(fp1,"# Survey Names\n");
            for (i = 0; i < NIndex; i++)
                fprintf(fp1,"#$INDEX-%d\n",i+1);

            fclose(fp1);
        }

        long LaunchAsap(char *fn)
        {
            char cmdline[MAXBUF];
            char xbuff[MCLEN];
            char xname[FILBUF];
            char zname[FILBUF];
            char *c;
            long k, n;
            char qt[2];

            qt[0] = '\x22';
            qt[1] = '\0';

            n = KThin * NBoot;

            strcpy(cmdline,WorkPath);
            strcat(cmdline,"\\");
            strcat(cmdline,"asap3.exe");

            sprintf(xbuff," -mcmc  %ld",n);
            strcat(cmdline,xbuff);

            sprintf(xbuff," -mcsave  %ld",KThin);
            strcat(cmdline,xbuff);

            sprintf(xbuff," -mcseed  %ld",ISeed);
            strcat(cmdline,xbuff);

            LowerCase(cmdline);

            k = CreateConsoleProcess(cmdline);
            if (k == 0)
                return k;

            strcpy(cmdline,WorkPath);
            strcat(cmdline,"\\");
            strcat(cmdline,"asap3.exe  -mceval \n");

            k = CreateConsoleProcess(cmdline);
            if (k == 0)
                return k;


            strcpy(xname,fn);

            if ((c = strrchr(xname,'.')) != NULL)
                *c = '\0';

            strcat(xname,"_asap.dat");

            strcpy(zname,WorkPath);
            strcat(zname,"\\asap3.dat");

            rename(zname,xname);

            if ((c = strrchr(xname,'.')) != NULL)
                *c = '\0';

            strcat(xname,".rep");

            if ((c = strrchr(zname,'.')) != NULL)
                *c = '\0';

            strcat(zname,".rep");

            rename(zname,xname);


            if ((c = strrchr(xname,'.')) != NULL)
                *c = '\0';

            strcat(xname,".bsn");

            if ((c = strrchr(zname,'.')) != NULL)
                *c = '\0';

            strcat(zname,".bsn");

            rename(zname,xname);

            return k;
        }

        void ScanAsapResults(char *fn)
        {
            char buffer[MAXBUF];
            char xname[FILBUF];
            char *c;
            char *tok;
            long i, j;
            short flag;

            strcpy(xname,fn);
            if ((c = strrchr(xname,'.')) != NULL)
                *c = '\0';

            strcat(xname,"_asap.rep");

            if ((fp1 = fopen(xname,"r")) == NULL)
            {
                fprintf(stderr,"Unable to Open ASAP Report File: %s\n",xname);
                exit(1);
            }

            flag = 0;
            while (!feof(fp1))
            {
                fgets(buffer,MAXBUF,fp1);
                if (strstr(buffer,"Population Numbers at the Start of the Year"))
                {
                    flag = 1;
                    break;
                }
            }

            if (!flag)
            {
                fprintf(stderr,"Missing Data From ASAP Report File\n");
                exit(1);
            }


            for (i = 0; i < NYears; i++)
            {
                fgets(buffer,MAXBUF,fp1);
                tok = strtok(buffer," \t\r\n");
                for (j = 0; j < KAges; j++)
                    StockAsap[i][j] = atof(tok);

                tok = strtok(NULL," \t\r\n");
            }


            flag = 0;
            while (!feof(fp1))
            {
                fgets(buffer,MAXBUF,fp1);
                if (strstr(buffer,"Biomass Time Series"))
                {
                    flag = 1;
                    break;
                }

            }

            if (!flag)
            {
                fprintf(stderr,"Missing Data From ASAP Report File\n");
                exit(1);
            }

            for (i = 0; i < NYears+1; i++)
                fgets(buffer,MAXBUF,fp1);

            tok = strtok(buffer," \t\r\n");
            tok = strtok(NULL," \t\r\n");
            tok = strtok(NULL," \t\r\n");

            SSBAsap = atof(tok);

            SSBAsap = SSBAsap / 1000.0;

            fclose(fp1);

        }

        void ReadAspicTemplateFile(char *fn)
        {
            char buffer[MAXBUF];
            char *c;
            char *tok;
            short i;

            if ((c = strrchr(fn,'.')) != NULL)
                *c = '\0';

            strcat(fn,".tx4");

            if ((fp1 = fopen(fn,"r")) == NULL)
            {
                fprintf(stderr,"Unable to Open Template File: %s\n",fn);
                exit(1);
            }


            AspicF = (double *) calloc(MaxYears+1,sizeof(double));
            AspicB = (double *) calloc(MaxYears+1,sizeof(double));
            AspicY = (double *) calloc(MaxYears+1,sizeof(double));


            fgets(buffer,MAXBUF,fp1);

            if ((c = strchr(buffer,'\n')) != NULL)
                *c = '\0';

            strcpy(AspicPath,buffer);

            strcpy(AspicpPath,buffer);

            if ((c = strchr(AspicpPath,'.')) != NULL)
                *c = '\0';

            strcat(AspicpPath,"p.exe");


            for (i = 0; i < NSTRING; i++)
            {
                fgets(buffer,MAXBUF,fp1);
                PString[i] = _strdup(buffer);
            }

            fgets(buffer,MAXBUF,fp1);
            tok = strtok(buffer," \t\r\n");
            cvAspicErr = atof(tok);
            tok = strtok(NULL," \t\r\n");
            biasAspicErr = atof(tok);

            fclose(fp1);

        }

        void WriteAspicInputFile(char *fn)
        {
            char xname[FILBUF];
            char *c;
            char quote = '\x22';
            short i, j, k;
            double xt, xx;

            strcpy(xname,fn);

            if ((c = strrchr(xname,'.')) != NULL)
                *c = '\0';

            strcat(xname,"_aspic.inp");

            if ((fp1 = fopen(xname,"w")) == NULL)
            {
                fprintf(stderr,"Unable to Open ASPIC Input File: %s\n",xname);
                exit(1);
            }

            for (i = 0; i < NSTRING; i++)
                fputs(PString[i],fp1);

            fprintf(fp1,"%ld\n",NYears);

            fprintf(fp1,"%cCPUE  & Yield%c\n",quote,quote);

            fprintf(fp1,"CC\n");

            for (i = 0; i < NYears; i++)
            {
                xt = 0.0;
                for (j = 0; j < MaxAge; j++)
                    xt += CatchSamples[j][i] * CatchWeight[j][i];
                xx = 0.0;
                for (j = 0; j < MaxAge; j++)
                    xx += SurveySamples[j][i][0] * StockWeight[j][i];

                fprintf(fp1,"%d    %18.8E  %18.8E\n",NFYear+i,xx,xt);
            }
            for (k = 1; k < NSurvey; k++)
            {
                fprintf(fp1,"%cSurvey #%d Biomass Index%c\n",quote,k+1,quote);
                fprintf(fp1,"I0\n");
                for (i = 0; i < NYears; i++)
                {

                    xx = 0.0;
                    for (j = 0; j < MaxAge; j++)
                        xx += SurveySamples[j][i][k] * StockWeight[j][i];

                    fprintf(fp1,"%d   %18.8E\n",NFYear+i,xx);
                }
            }

            fclose(fp1);
        }

        void CaptureAspicResults(char *fn)
        {
            char buffer[MAXBUF];
            char xname[FILBUF];
            char *tok, *c;
            short i;

            /* Open ASPIC Output File */

            strcpy(xname,fn);

            if ((c = strrchr(xname,'.')) != NULL)
                *c = '\0';

            strcat(xname,"_aspic.bot");

            if ((fp1 = fopen(xname,"r")) == NULL)
            {
                fprintf(stderr,"Unable to Open ASPIC Output File: %s\n",xname);
                fprintf(stderr,"Iteration Dropped\n");
                return;
            }

            while (!feof(fp1))
            {

                fgets(buffer,MAXBUF-1,fp1);

                if (strstr(buffer,"ESTIMATED POPULATION TRAJECTORY"))
                {
                    for (i = 0; i < 6; i++)
                        fgets(buffer,MAXBUF-1,fp1);

                    for (i = 0; i < NYears; i++)
                    {
                        fgets(buffer,MAXBUF-1,fp1);
                        tok = strtok(buffer," \t\r\n");
                        tok = strtok(NULL," \t\r\n");
                        tok = strtok(NULL," \t\r\n");
                        AspicF[i] = atof(tok);
                        tok = strtok(NULL," \t\r\n");
                        AspicB[i] = atof(tok);
                        tok = strtok(NULL," \t\r\n");
                        tok = strtok(NULL," \t\r\n");
                        tok = strtok(NULL," \t\r\n");
                        AspicY[i] = atof(tok);
                    }
                }

            }

            fclose(fp1);

        }

        void WriteAspicPFile(char *fn)
        {
            char xname[FILBUF];
            char *c;
            char quote = '\x22';
            double sd, zscore, tl;
            double xx;

            strcpy(xname,fn);

            if ((c = strrchr(xname,'.')) != NULL)
                *c = '\0';

            strcat(xname,"_aspic.ctl");

            if ((fp1 = fopen(xname,"w")) == NULL)
            {
                fprintf(stderr,"Unable to Open ASPIC Projection Input File: %s\n",xname);
                exit(1);
            }


            if ((c = strrchr(xname,'.')) != NULL)
                *c = '\0';

            strcat(xname,".bio");

            fprintf(fp1,"%cMSE Version 4.1%c\n",quote,quote);

            fprintf(fp1,"%-s\n",xname);

            fprintf(fp1,"XX\n");

            fprintf(fp1,"PC  1\n");

            fprintf(fp1,"0\n");

            fprintf(fp1,"2\n");

            sd = sqrt(log(cvAspicErr * cvAspicErr + 1.0));
            zscore = gsl_box_muller();
            xx = exp(sd * zscore);
            tl = TargetLand * (1.0 + biasAspicErr) * xx;

            fprintf(fp1,"%10.3f   Y\n",tl);

            fprintf(fp1,"%10.4f   F\n",FTarg/FullF[NYears-1]);

            fclose(fp1);
        }

        long LaunchAspicProj(char *fn)
        {
            char xname[FILBUF];
            char cmdline[MAXBUF];
            char *c;
            long k;
            char qt[2];

            qt[0] = '\x22';
            qt[1] = '\0';

            strcpy(xname,fn);

            if ((c = strrchr(xname,'.')) != NULL)
                *c = '\0';

            strcat(xname,"_aspic.ctl");


            strcpy(cmdline,AspicpPath);
            strcat(cmdline," ");
            strcat(cmdline,qt);
            strcat(cmdline,xname);
            strcat(cmdline,qt);
            LowerCase(cmdline);

            k = CreateConsoleProcess(cmdline);

            return k;
        }

        void ScanAspicProj(char *fn)
        {
            char xname[FILBUF];
            char buffer[MAXBUF];
            char *c;
            char *tok;
            short i, flag;

            strcpy(xname,fn);

            if ((c = strrchr(xname,'.')) != NULL)
                *c = '\0';

            strcat(xname,"_aspic.prj");

            if ((fp1 = fopen(xname,"r")) == NULL)
            {
                fprintf(stderr,"Unable to Open ASPIC Projection Input File: %s\n",xname);
                exit(1);
            }

            flag = 0;
            while(!feof(fp1))
            {
                fgets(buffer,MAXBUF,fp1);
                if (strstr(buffer,"TABLE OF PROJECTED YIELDS"))
                {
                    flag = 1;
                    break;
                }
            }

            if (!flag)
            {
                fprintf(stderr,"Missing Data From Aspic Projection Report\n");
                exit(1);
            }

            for (i = 0; i < 3; i++)
                fgets(buffer,MAXBUF,fp1);

            tok = strtok(buffer, " \t\r\n");
            tok = strtok(NULL, " \t\r\n");

            ProjLand = atof(tok);

            ProjLand = ProjLand / 1000.0;

            fclose(fp1);

        }

        long LaunchAspic(char *fn)
        {
            char xname[FILBUF];
            char cmdline[MAXBUF];
            char *c;
            long k;
            char qt[2];

            qt[0] = '\x22';
            qt[1] = '\0';

            strcpy(xname,fn);

            if ((c = strrchr(xname,'.')) != NULL)
                *c = '\0';

            strcat(xname,"_aspic.inp");


            strcpy(cmdline,AspicPath);
            strcat(cmdline," ");
            strcat(cmdline,qt);
            strcat(cmdline,xname);
            strcat(cmdline,qt);
            LowerCase(cmdline);

            k = CreateConsoleProcess(cmdline);

            return k;
        }

        void ExpandData()
        {
            long i,k;

            for (i = 0; i < MaxAge; i++)
            {
                StockWeight[i][MaxYears]   = StockWeight[i][MaxYears-1];
                StockWeight[i][MaxYears+1] = StockWeight[i][MaxYears-1];
                SpawnWeight[i][MaxYears]   = SpawnWeight[i][MaxYears-1];
                CatchWeight[i][MaxYears]   = CatchWeight[i][MaxYears-1];
                DiscFrac[i][MaxYears]      = DiscFrac[i][MaxYears-1];
                FSelec[i][MaxYears]        = FSelec[i][MaxYears-1];
                NatMort[i][MaxYears]       = NatMort[i][MaxYears-1];
                Mature[i][MaxYears]        = Mature[i][MaxYears-1];
            }
            TF[MaxYears] = TF[MaxYears-1];
            TM[MaxYears] = TM[MaxYears-1];

            for (k = 0; k < NSurvey; k++)
            {
                SurveyCV[MaxYears][k]   = SurveyCV[MaxYears-1][k];
                SurveyDraw[MaxYears][k] = SurveyDraw[MaxYears-1][k];
                SurveyQ[MaxYears][k]    = SurveyQ[MaxYears-1][k];

                for (i = 0; i < MaxAge; i++)
                    SurveySelec[i][MaxYears][k] = SurveySelec[i][MaxYears-1][k];
            }
        }

    };

}

#endif /* MSE_HPP */

