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
 */

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
#define   XTOL    1.0E-15
#define   XXLW         6
#define   XXNM         2
#define   XXDF         2

#define   XXSD     10000

#define   XLBOUND  0.0001
#define   XUBOUND  5.0

#define   TAGSIZE  16

#define   MXITER  100

#define   FORMFEED  '\x0C'

#define   NOPT  10

#define   NSTRING  21

#define   NSELMOD  6
#define   NPHASE   8

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

        FILE *fp1, *fp2, *fp3, *fp4, *fp5, *fp6;

        short assessment_frequency;

        short NYears;
        short NAges;
        short KAges;
        short NBins;
        short NFYear;
        short NXYear;
        short NFAge = 1;
        short NAgePlus;
        short MaxAge;
        short MinLen;
        short MaxLen;
        short NMarket;
        short NSurvey;
        short NIndex;
        short DiscFlag;
        short SSBFlag;
        short RecruitFlag;
        short NAgeErr;


        short BaseYears;
        short MXYear;
        short MaxYears;

        short KFish;

        short ModelFlag;

        short *SampFlag;
        short *MisAgeFlag;

        short *AgeErrStart;
        short *AgeErrEnd;


        long ISeed;
        long KSeed;
        long NIter;

        long BaseSeed;

        long NSample;

        long *KSamples;

        double InitPop;
        double BasePop;
        double PrevPop;

        double LInfMale;
        double KParmMale;
        double TZeroMale;
        double GrowthStdDevMale;

        double LInfFemale;
        double KParmFemale;
        double TZeroFemale;
        double GrowthStdDevFemale;

        double FMAlpha;
        double FMBeta;

        double FFAlpha;
        double FFBeta;

        double TF;
        double TM;

        double RAlpha;
        double RBeta;
        double RKPar;

        double CVRecruit;
        double CVLand;
        double PopCv;

        double SProdDelta;

        double *LenStdDevMale;
        double *LenStdDevFemale;

        double *SampFrac;

        double *InitMale;
        double *InitFemale;
        double *InitStock;

        double *GrowthIntMale;
        double *GrowthIntFemale;

        double *GrowthAvgMale;
        double *GrowthAvgFemale;

        double *InitRecruits;
        double *Recruits;
        double *RecSexFemale;
        double *FracMatMale;
        double *FracMatFemale;

        double *FFull;

        double *RecDev;

        double *SSB;

        double *StockBiomass;
        double *CatchBiomass;
        double *SurplusProd;


        double *TotalLand;
        double *TotalDisc;

        double *LandDev;

        double *PopDev;

        double *DFrac;


        double *SampleLength;
        double *SelectLength;

        double *SurveyStock;
        double *SurveySamples;
        double *SurveyDist;

        double *WorkVec;

        double **AgeProbMale;
        double **AgeProbFemale;
        double **MaleBins;
        double **FemaleBins;
        double **StockBins;

        double **StockAge;
        double **StockLength;
        double **StockWeights;
        double **StockAgeMax;

        double **MaleAge;
        double **MaleLength;
        double **MaleWeights;

        double **FemaleAge;
        double **FemaleLength;
        double **FemaleWeights;

        double **MaleGrowth;
        double **FemaleGrowth;

        double **LenWtCoeff;
        double **DiscFrac;

        double **NMort;
        double **NMSelect;
        double **FSelect;

        double **MaleCatch;
        double **MaleDisc;
        double **MaleSurvive;
        double **MaleSpStock;

        double **FemaleCatch;
        double **FemaleDisc;
        double **FemaleSurvive;
        double **FemaleSpStock;

        double **CatchBins;
        double **DiscBins;
        double **SpawnStock;
        double **CatchAgeKey;
        double **StockAgeKey;

        double **CatchLength;
        double **DiscLength;

        double **MaleCatchLength;
        double **FemaleCatchLength;

        double **CatchAge;
        double **DiscAge;
        double **SpStockAge;
        double **SSBAge;

        double **MaleCatchAge;
        double **FemaleCatchAge;

        double **MaleSpawnAge;
        double **FemaleSpawnAge;

        double **CatchWeights;
        double **MaleCatchWeights;
        double **FemaleCatchWeights;

        double **MaleSpawnWeights;
        double **FemaleSpawnWeights;
        double **SpStockWeights;

        double **Mature;

        double **Survivorx;

        double **FCalc;
        double **MCalc;

        double **MarketProb;
        double **MarketCatch;
        double **MarketLand;
        double **MarketDist;
        double **MarketSampleRate;
        double **MarketSampleSize;
        double **MarketTarg;

        double **SampleBins;
        double **SampleAges;

        double **MisAge;

        double **AgeErr;


        double **AgeLengthKey;
        double **AgeKeyCalc;

        double **ExpandLength;
        double **ExpandCatch;
        double **ExpandWeights;

        double **SurveySelect;
        double **SurveyQ;
        double **SurveyLength;
        double **SurveyAge;
        double **SurveyBins;
        double **SelectAge;

        double **IndexValues;

        double **WorkData;


        char **MarketTags;

        struct si {
            char Name[TAGSIZE];
            short StartAge;
            short EndAge;
            short Type;
            short SurveyIndex;
        };

        struct si *IndexData;


        /* MSR */

        short MRule;
        short ProjFlag;

        double Fmsy;
        double Bmsy;
        double Blim;
        double TargetF;
        double TargetSSB;
        double TargetLand;
        double ProjLand;
        double ProjF;
        double DeltaLand;
        double PrevLand;

        double *RuleDev;
        double CvRule;
        double BiasRule;
        double ProjSave;

        double *SSBDev;
        double CvSSB;
        double BiasSSB;

        double *FProjDev;
        double CvFProj;
        double BiasFProj;



        /* VPA */

        char VPAPath[FILBUF];

        short Fold[2];

        double *PRVec;
        double *StockEst;
        double NMortEst;

        double *CatchFactors;

        double **StockVPA;
        double **FVPA;
        double **SSBVPA;

        double **VPASurveyFactors;

        short CMStart;
        short CMEnd;

        double CMxlb;
        double CMxub;

        short *AgeEst;

        short OptionFlags[NOPT];

        short NRange[2];

        short SurveyFacFlag;
        short CMultFlag;

        short NBoot;

        /* AgePro */

        char AgeProPath[FILBUF];

        short NAgeProSim;

        double AgeProCV;
        double AgeProBias;

        double *AgeProDev;

        /*  ASPIC */

        char AspicPath[FILBUF];
        char AspicpPath[FILBUF];

        double *AspicF;
        double *AspicB;
        double *AspicY;

        char *PString[NSTRING];

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
            char fname[FILBUF];
            char *fn;
            char *c;
            long iopt = 5;
            long kx, nx;
            long Iter;
            short k, n;
            double zero = 0.0;
            double zx, sd, xt, err;
            double pland, FPrev;
            double SSBEstim, FEstim;
            double FHigh, FLow;


            if (argc < 2) {
                fprintf(stderr, "Usage: mseproj  filename\n");
                exit(1);
            }

            fn = *++argv;

            strcpy(fname, fn);


            if ((fp1 = fopen(fname, "r")) == NULL) {
                fprintf(stderr, "Unable to Open Input File: %s\n", fname);
                exit(1);
            }

            if ((c = strrchr(fname, '.')) != NULL)
                *c = '\0';

            strcat(fname, ".rep");

            if ((fp2 = fopen(fname, "w")) == NULL) {
                fprintf(stderr, "Unable to Open Report File: %s\n", fname);
                exit(1);
            }


            if ((c = strrchr(fname, '.')) != NULL)
                *c = '\0';

            strcat(fname, ".xx1");

            if ((fp3 = fopen(fname, "w")) == NULL) {
                fprintf(stderr, "Unable to Open Report File: %s\n", fname);
                exit(1);
            }


            if ((c = strrchr(fname, '.')) != NULL)
                *c = '\0';

            strcat(fname, ".xx2");

            if ((fp4 = fopen(fname, "w")) == NULL) {
                fprintf(stderr, "Unable to Open Report File: %s\n", fname);
                exit(1);
            }


            if ((c = strrchr(fname, '.')) != NULL)
                *c = '\0';

            strcat(fname, ".xx3");

            if ((fp5 = fopen(fname, "w")) == NULL) {
                fprintf(stderr, "Unable to Open Report File: %s\n", fname);
                exit(1);
            }



            if ((c = strrchr(fname, '.')) != NULL)
                *c = '\0';


            /* Read Input Data File */

            ReadDataFile();
            std::cout << "assessment frequency = " << this->assessment_frequency << "\n";
            std::cout << "done.\n";
            //            exit(0);
            BaseSeed = ISeed;

            BaseYears = NYears;

            BasePop = InitPop;

            PrevPop = InitPop;

            if (ModelFlag == 1) {
                ReadVPATemplateFile(fname);

                SetupAgepro(fname);

                ReadAgeproTemplate(fname);
            } else if (ModelFlag == 2)
                ReadAspicTemplateFile(fname);




            fprintf(fp2, "Management Strategy Evaluation Version 3.4\n\n");

            /* Initialize Random Number Generator */

            RNOPT(&iopt);

            RNSET(&ISeed);

            DRNNOR(&NIter, PopDev);

            RNGET(&KSeed);

            /*  Iterate on Initial Population */

            for (Iter = 0; Iter < NIter; Iter++) {

                printf("Iteration # %ld\n", Iter + 1);
                fprintf(fp2, "Iteration # %ld\n\n", Iter + 1);

                NYears = BaseYears;

                zx = PopDev[Iter];
                sd = sqrt(log(PopCv * PopCv + 1.0));

                InitPop = BasePop * exp(zx * sd);

                /* Create Initial Population */

                CreateInitialDistribution();

                /* Create Growth Projection Matrix for Each Sex */

                CreateGrowthMatrix();

                /* Calculate Fraction Mature at Length */

                CalcMaturityFraction();

                /* Get Random Deviates  */

                nx = (long) MaxYears;

                DRNNOR(&nx, RecDev);
                DRNNOR(&nx, RuleDev);
                DRNNOR(&nx, AgeProDev);
                DRNNOR(&nx, SSBDev);
                DRNNOR(&nx, FProjDev);

                /* Calculate True State of Nature Across Time Series */

                for (k = 0; k < NYears; k++) {

                    if (!(k % this->assessment_frequency) == 0) {
                                                continue;
                    }
                    /* Stock Age Key */

                    CalcStockAgeKey(k);

                    /* Calculate Stock Weights */

                    CalcStockWeights(k);

                    /* Apply Mortality */

                    CalculateCatch(k);

                    /* Calculate Catch Weights */

                    CalcCatchWeights(k);

                    /* Calculate Spawning Stock Weights & SSB */

                    CalcSpStockWeights(k);

                    /* Calculate Maturity at Age */

                    CalcMaturity(k);

                    /* Apply Growth to Survivors & Advance Stock Age */

                    ApplyGrowthMatrix();

                    /* Insert Recruits */

                    CalcRecruits(k);

                    /* Save Catch & Stock Results for Year k */

                    SaveAnnualData(k);

                    /* Calculate Landings & Discards */

                    CalcLandings(k);

                    CalcMarketSamples(k);

                }

                /* Generate Samples for Base Years */

                for (k = 0; k < NYears; k++) {
                    /* Get Samples by Market Category */
                    if (!(k % this->assessment_frequency) == 0) {
                                                continue;
                    }
                    GetMarketSamples(k);

                    /* Get Age Length Key for Expanded Catch */

                    if (NAgeErr > 0)
                        CopyMisAgeMatrix(k);

                    if (SampFlag[0] == 0)
                        OriginalCatchAgeKey(k);
                    else
                        SampleCatchAgeKey(k);

                    /* Expand Samples by Market */

                    ExpandSamples(k);

                    CalcExpandWeights(k);
                }



                /* Stock Age Key */

                CalcStockAgeKey(NYears);

                /* Calculate Stock Weights in Terminal Year + 1 */

                CalcStockWeights(NYears);


                /* Get Survey Samples */

                n = NYears + 1;
                ZeroMatrix(IndexValues, n, NIndex);

                for (n = 0; n < NSurvey; n++) {
                    for (k = 0; k < NYears + 1; k++) {
                        if (!(k % this->assessment_frequency) == 0) {
                                                continue;
                    }
                        GetSurveySamples(k, n);
                        ApplySurveyAgeKey(k, n);
                    }

                    CreateIndexData(n);

                }


                if (ModelFlag == 0) {
                    fprintf(fp2, "Year      ");
                    fprintf(fp2, "Total Catch   ");
                    fprintf(fp2, "   SSB Estim     ");
                    fprintf(fp2, " Target F       ");
                    fprintf(fp2, " F Full (N+1) \n\n");
                } else {
                    fprintf(fp2, "Year      ");
                    fprintf(fp2, "Total Catch   ");
                    fprintf(fp2, "   SSB (VPA)         ");
                    fprintf(fp2, " Target F     ");
                    fprintf(fp2, "AgePro Catch     ");
                    fprintf(fp2, "Proj Catch        ");
                    fprintf(fp2, "Full F (N+1)\n\n");
                }


                ProjFlag = 0;

                ProjSave = 0.0;


                while (NYears < MaxYears + 1) {

                    kx = NYears * NMarket;

                    DRNNOR(&kx, LandDev);

                    if (ModelFlag == 0) {
                        TargetLand = TotalLand[NYears - 1] + TotalDisc[NYears - 1];

                        /* Apply Noise to SSB */

                        zx = SSBDev[NYears - 1];
                        sd = sqrt(log(CvSSB * CvSSB + 1.0));
                        err = exp(zx * sd);

                        SSBEstim = SSB[NYears - 1] * (1.0 + BiasSSB) * err;

                        SSBEstim = SSBEstim / 1000.;

                        /* Apply Management Rule */

                        if (SSBEstim <= Blim)
                            FEstim = 0.0;
                        else if (SSBEstim < Bmsy)
                            FEstim = SSBEstim * Fmsy / Bmsy;
                        else
                            FEstim = Fmsy;

                        /* Apply Noise to F */

                        zx = FProjDev[NYears - 1];
                        sd = sqrt(log(CvFProj * CvFProj + 1.0));
                        err = exp(zx * sd);


                        if (ProjFlag) {
                            if (MRule) {
                                pland = TargetLand * (1.0 + DeltaLand);
                                FHigh = CalculateProjectedF(pland);

                                pland = TargetLand * (1.0 - DeltaLand);
                                FLow = CalculateProjectedF(pland);

                                FFull[NYears] = FPrev * (1.0 + BiasFProj) * err;

                                if (FFull[NYears] > FHigh)
                                    FFull[NYears] = FHigh;
                                else if (FFull[NYears] < FLow)
                                    FFull[NYears] = FLow;
                            } else
                                FFull[NYears] = FPrev * (1.0 + BiasFProj) * err;
                        } else {
                            zx = RuleDev[NYears - 1];
                            sd = sqrt(log(CvRule * CvRule + 1.0));
                            err = exp(zx * sd);

                            pland = TargetLand * (1.0 + BiasRule) * err;

                            FFull[NYears] = CalculateProjectedF(pland);
                        }

                        fprintf(fp2, "%d %15.3f ", NFYear + NYears - 1, TargetLand);
                        fprintf(fp2, "%15.3f %15.4f %15.4f\n", SSBEstim, FEstim, FFull[NYears]);

                        fprintf(fp3, "%d %15.3f ", NFYear + NYears - 1, TargetLand);
                        fprintf(fp3, "%15.3f %15.4f %15.4f\n", SSBEstim, FEstim, FFull[NYears]);

                        printf("Year Completed %d\n", NFYear + NYears - 1);

                        ProjFlag = 1;

                        FPrev = FEstim;

                    } else if (ModelFlag == 1) {
                        /* Launch VPA */

                        strcpy(fname, fn);

                        WriteVPAInputFile(fname);

                        RemoveOutputFiles(fname);

                        kx = LaunchVPA(fname);

                        if (kx == 0) {
                            fprintf(stderr, "Error Launching VPA\n");
                            exit(1);
                        }

                        ScanVPAResults(fname);

                        /* Apply Management Rules */

                        TargetF = GetTargetF(&TargetSSB);

                        TargetLand = TotalLand[NYears - 1] + TotalDisc[NYears - 1];

                        /* Launch AgePro */

                        WriteAgeProInputFile(fname);

                        kx = LaunchAgePro();

                        if (kx == 0) {
                            fprintf(stderr, "Error Launching AgePro\n");
                            exit(1);
                        }

                        ScanAgeProReport(fname);


                        if (ProjFlag) {
                            pland = ProjSave;

                            if (MRule) {
                                xt = (pland - PrevLand) / PrevLand;
                                if (xt > DeltaLand)
                                    pland = PrevLand * (1.0 + DeltaLand);
                                else if (fabs(xt) > DeltaLand)
                                    pland = PrevLand * (1.0 - DeltaLand);
                            }

                        } else
                            pland = TargetLand;

                        /* Apply Implementation Error */

                        zx = RuleDev[NYears - 1];
                        sd = sqrt(log(CvRule * CvRule + 1.0));
                        err = exp(zx * sd);

                        pland = pland * (1.0 + BiasRule) * err;


                        ProjF = CalculateProjectedF(pland);


                        fprintf(fp2, "%d %15.3f ", NFYear + NYears - 1, TargetLand);
                        fprintf(fp2, "%15.3f %15.4f ", TargetSSB, TargetF);
                        fprintf(fp2, "%15.3f %15.3f %15.4f\n", ProjLand, pland, ProjF);

                        fprintf(fp3, "%d %15.3f ", NFYear + NYears - 1, TargetLand);
                        fprintf(fp3, "%15.3f %15.4f ", TargetSSB, TargetF);
                        fprintf(fp3, "%15.3f %15.3f %15.4f\n", ProjLand, pland, ProjF);

                        printf("Year Completed %d\n", NFYear + NYears - 1);

                        FFull[NYears] = ProjF;

                        ProjFlag = 1;

                        ProjSave = ProjLand;

                        PrevLand = pland;
                    } else if (ModelFlag == 2) {
                        /* Launch Aspic */

                        strcpy(fname, fn);

                        WriteAspicInputFile(fname);

                        RemoveOutputFiles(fname);

                        kx = LaunchAspic(fname);

                        if (kx == 0) {
                            fprintf(stderr, "Error Launching Aspic\n");
                            exit(1);
                        }

                        CaptureAspicResults(fname);

                        /* Apply Management Rules */

                        TargetF = GetAspicTargetF(&TargetSSB);

                        TargetLand = TotalLand[NYears - 1] + TotalDisc[NYears - 1];


                        /* Launch  Aspic Projection */

                        WriteAspicPFile(fname);

                        kx = LaunchAspicProj(fname);

                        if (kx == 0) {
                            fprintf(stderr, "Error Launching Aspic Projection\n");
                            exit(1);
                        }

                        ScanAspicProj(fname);

                        if (ProjFlag) {
                            pland = ProjSave;

                            if (MRule) {
                                xt = (pland - PrevLand) / PrevLand;
                                if (xt > DeltaLand)
                                    pland = PrevLand * (1.0 + DeltaLand);
                                else if (fabs(xt) > DeltaLand)
                                    pland = PrevLand * (1.0 - DeltaLand);
                            }

                        } else
                            pland = TargetLand;

                        /* Apply Implementation Error */

                        zx = RuleDev[NYears - 1];
                        sd = sqrt(log(CvRule * CvRule + 1.0));
                        err = exp(zx * sd);

                        pland = pland * (1.0 + BiasRule) * err;


                        ProjF = CalculateProjectedF(pland);


                        fprintf(fp2, "%d %15.3f ", NFYear + NYears - 1, TargetLand);
                        fprintf(fp2, "%15.3f %15.4f ", TargetSSB, TargetF);
                        fprintf(fp2, "%15.3f %15.3f %15.4f\n", ProjLand, pland, ProjF);

                        fprintf(fp3, "%d %15.3f ", NFYear + NYears - 1, TargetLand);
                        fprintf(fp3, "%15.3f %15.4f ", TargetSSB, TargetF);
                        fprintf(fp3, "%15.3f %15.3f %15.4f\n", ProjLand, pland, ProjF);

                        printf("Year Completed %d\n", NFYear + NYears - 1);

                        FFull[NYears] = ProjF;

                        ProjFlag = 1;

                        ProjSave = ProjLand;

                        PrevLand = pland;

                    }




                    if (NYears == MaxYears)
                        break;

                    /* Apply Mortality */

                    CalculateCatch(NYears);

                    /* Calculate Catch Weights */

                    CalcCatchWeights(NYears);

                    /* Calculate Spawning Stock Weights & SSB */

                    CalcSpStockWeights(NYears);

                    /* Calculate Maturity at Age */

                    CalcMaturity(NYears);

                    /* Apply Growth to Survivors & Advance Stock Age */

                    ApplyGrowthMatrix();

                    /* Insert Recruits */

                    CalcRecruits(NYears);

                    /* Save Catch & Stock Results for Year k */

                    SaveAnnualData(NYears);

                    /* Calculate Landings & Discards */

                    CalcLandings(NYears);

                    CalcMarketSamples(NYears);


                    /* Use Alternate Random Sequence for Sampling*/

                    RNGET(&ISeed);
                    RNSET(&KSeed);

                    /* Get Samples by Market Category */

                    GetMarketSamples(NYears);

                    /* Get Age Length Key for Expanded Catch */

                    if (NAgeErr > 0)
                        CopyMisAgeMatrix(NYears);

                    if (SampFlag[0] == 0)
                        OriginalCatchAgeKey(NYears);
                    else
                        SampleCatchAgeKey(NYears);

                    /* Expand Samples by Market */

                    ExpandSamples(NYears);

                    CalcExpandWeights(NYears);



                    NYears++;


                    /* Stock Age Key */

                    CalcStockAgeKey(NYears);

                    /* Calculate Stock Weights in Terminal Year + 1 */

                    CalcStockWeights(NYears);



                    /* Get Survey Samples */


                    for (n = 0; n < NSurvey; n++) {

                        GetSurveySamples(NYears, n);
                        ApplySurveyAgeKey(NYears, n);

                        AddIndexData(n);

                    }


                    RNGET(&KSeed);
                    RNSET(&ISeed);

                }

                fprintf(fp2, "\n");
                fprintf(fp3, "\n");


                /* Basic Reports */

                BasicReports();

                PrevPop = InitPop;

            }

            fclose(fp2);
            fclose(fp3);
            fclose(fp4);
            fclose(fp5);

        }

        void ReadDataFile() {
            char buffer[MAXBUF];
            char *c;
            char *tok;
            short i, j, k, n;



            fgets(buffer, MAXBUF - 1, fp1);
            if (!strstr(buffer, "MSE VERSION 3.1")) {
                fprintf(stderr, "Invalid Input Data File\n");
                exit(1);
            }

            DiscFlag = 0;

            SProdDelta = 1.0;

            ModelFlag = 0;

            while (!feof(fp1)) {
                fgets(buffer, MAXBUF - 1, fp1);

                if (strstr(buffer, "GENERAL DATA")) {
                    fgets(buffer, MAXBUF - 1, fp1);
                    tok = strtok(buffer, " \t\r\n");
                    NFYear = atoi(tok);
                    tok = strtok(NULL, " \t\r\n");
                    NXYear = atoi(tok);
                    tok = strtok(NULL, " \t\r\n");
                    MXYear = atoi(tok);
                    tok = strtok(NULL, " \t\r\n");
                    NAgePlus = atoi(tok);
                    tok = strtok(NULL, " \t\r\n");
                    MaxAge = atoi(tok);
                    tok = strtok(NULL, " \t\r\n");
                    MinLen = atoi(tok);
                    tok = strtok(NULL, " \t\r\n");
                    MaxLen = atoi(tok);
                    tok = strtok(NULL, " \t\r\n");
                    NMarket = atoi(tok);
                    tok = strtok(NULL, " \t\r\n");
                    NSurvey = atoi(tok);
                    tok = strtok(NULL, " \t\r\n");
                    NIndex = atoi(tok);
                    tok = strtok(NULL, " \t\r\n");
                    NAgeErr = atoi(tok);
                    tok = strtok(NULL, " \t\r\n");
                    ISeed = atol(tok);
                    tok = strtok(NULL, " \t\r\n");
                    NIter = atol(tok);
                    tok = strtok(NULL, " \t\r\n");
                    ModelFlag = atoi(tok);

                    AllocData();
                } else if (strstr(buffer, "ASSESSMENT FREQUENCY")) {
                    fgets(buffer, MAXBUF - 1, fp1);
                    tok = strtok(buffer, " \t\r\n");
                    this->assessment_frequency = atof(tok);
                } else if (strstr(buffer, "INITPOP")) {
                    fgets(buffer, MAXBUF - 1, fp1);
                    tok = strtok(buffer, " \t\r\n");
                    PopCv = atof(tok);

                    InitPop = 0.0;
                    for (i = 0; i < NAges; i++) {
                        fgets(buffer, MAXBUF - 1, fp1);
                        tok = strtok(buffer, " \t\r\n");
                        InitMale[i] = atof(tok);
                        tok = strtok(NULL, " \t\r\n");
                        InitFemale[i] = atof(tok);
                        tok = strtok(NULL, " \t\r\n");
                        LenStdDevMale[i] = atof(tok);
                        tok = strtok(NULL, " \t\r\n");
                        LenStdDevFemale[i] = atof(tok);

                        InitStock[i] = InitMale[i] + InitFemale[i];
                        InitPop += InitStock[i];
                    }


                } else if (strstr(buffer, "GROWTH")) {
                    /* Male */

                    fgets(buffer, MAXBUF - 1, fp1);
                    tok = strtok(buffer, " \t\r\n");
                    LInfMale = atof(tok);
                    tok = strtok(NULL, " \t\r\n");
                    KParmMale = atof(tok);
                    tok = strtok(NULL, " \t\r\n");
                    TZeroMale = atof(tok);

                    fgets(buffer, MAXBUF - 1, fp1);
                    tok = strtok(buffer, " \t\r\n");
                    GrowthStdDevMale = atof(tok);

                    /* Female */

                    fgets(buffer, MAXBUF - 1, fp1);
                    tok = strtok(buffer, " \t\r\n");
                    LInfFemale = atof(tok);
                    tok = strtok(NULL, " \t\r\n");
                    KParmFemale = atof(tok);
                    tok = strtok(NULL, " \t\r\n");
                    TZeroFemale = atof(tok);

                    fgets(buffer, MAXBUF - 1, fp1);
                    tok = strtok(buffer, " \t\r\n");
                    GrowthStdDevFemale = atof(tok);



                } else if (strstr(buffer, "MATURE")) {
                    /* SSB Flag */

                    fgets(buffer, MAXBUF - 1, fp1);
                    tok = strtok(buffer, " \t\r\n");
                    SSBFlag = atoi(tok);

                    fgets(buffer, MAXBUF - 1, fp1);
                    tok = strtok(buffer, " \t\r\n");
                    FMAlpha = atof(tok);
                    tok = strtok(NULL, " \t\r\n");
                    FMBeta = atof(tok);

                    fgets(buffer, MAXBUF - 1, fp1);
                    tok = strtok(buffer, " \t\r\n");
                    FFAlpha = atof(tok);
                    tok = strtok(NULL, " \t\r\n");
                    FFBeta = atof(tok);

                } else if (strstr(buffer, "TFM")) {
                    fgets(buffer, MAXBUF - 1, fp1);
                    tok = strtok(buffer, " \t\r\n");
                    TF = atof(tok);
                    tok = strtok(NULL, " \t\r\n");
                    TM = atof(tok);
                } else if (strstr(buffer, "LEN-WT")) {
                    for (i = 0; i < MaxYears + 1; i++) {
                        fgets(buffer, MAXBUF - 1, fp1);
                        tok = strtok(buffer, " \t\r\n");
                        for (j = 0; j < XXLW; j++) {
                            LenWtCoeff[i][j] = atof(tok);
                            tok = strtok(NULL, " \t\r\n");
                        }
                    }
                } else if (strstr(buffer, "DISCARD")) {
                    DiscFlag = 1;
                    for (i = 0; i < MaxYears; i++) {
                        fgets(buffer, MAXBUF - 1, fp1);
                        tok = strtok(buffer, " \t\r\n");
                        for (j = 0; j < XXDF; j++) {
                            DiscFrac[i][j] = atof(tok);
                            tok = strtok(NULL, " \t\r\n");
                        }
                    }
                } else if (strstr(buffer, "FFULL")) {
                    for (i = 0; i < NYears; i++) {
                        fgets(buffer, MAXBUF - 1, fp1);
                        tok = strtok(buffer, " \t\r\n");
                        FFull[i] = atof(tok);
                    }
                } else if (strstr(buffer, "FSELECT")) {
                    for (i = 0; i < NBins; i++) {
                        fgets(buffer, MAXBUF - 1, fp1);
                        tok = strtok(buffer, " \t\r\n");
                        for (j = 0; j < MaxYears; j++) {
                            FSelect[i][j] = atof(tok);
                            tok = strtok(NULL, " \t\r\n");
                        }
                    }
                } else if (strstr(buffer, "NSELECT")) {
                    for (i = 0; i < NAges; i++) {
                        fgets(buffer, MAXBUF - 1, fp1);
                        tok = strtok(buffer, " \t\r\n");
                        for (j = 0; j < XXNM; j++) {
                            NMSelect[i][j] = atof(tok);
                            tok = strtok(NULL, " \t\r\n");
                        }
                    }
                } else if (strstr(buffer, "NMORT")) {
                    for (i = 0; i < MaxYears; i++) {
                        fgets(buffer, MAXBUF - 1, fp1);
                        tok = strtok(buffer, " \t\r\n");
                        for (j = 0; j < XXNM; j++) {
                            NMort[i][j] = atof(tok);
                            tok = strtok(NULL, " \t\r\n");
                        }
                    }
                } else if (strstr(buffer, "RECRUITS")) {
                    fgets(buffer, MAXBUF - 1, fp1);
                    tok = strtok(buffer, " \t\r\n");
                    RecruitFlag = atoi(tok);

                    fgets(buffer, MAXBUF - 1, fp1);
                    tok = strtok(buffer, " \t\r\n");
                    CVRecruit = atof(tok);

                    for (i = 0; i < MaxYears; i++) {
                        fgets(buffer, MAXBUF - 1, fp1);
                        tok = strtok(buffer, " \t\r\n");
                        RecSexFemale[i] = atof(tok);
                    }

                    if (RecruitFlag == 0) {
                        for (i = 0; i < MaxYears; i++) {
                            fgets(buffer, MAXBUF - 1, fp1);
                            tok = strtok(buffer, " \t\r\n");
                            InitRecruits[i] = atof(tok);
                        }
                    } else {
                        fgets(buffer, MAXBUF - 1, fp1);
                        tok = strtok(buffer, " \t\r\n");
                        RAlpha = atof(tok);
                        tok = strtok(NULL, " \t\r\n");
                        RBeta = atof(tok);
                        tok = strtok(NULL, " \t\r\n");
                        RKPar = atof(tok);
                    }
                } else if (strstr(buffer, "MARKET")) {
                    fgets(buffer, MAXBUF - 1, fp1);
                    tok = strtok(buffer, " \t\r\n");
                    CVLand = atof(tok);
                    for (i = 0; i < NMarket; i++) {
                        fgets(buffer, MAXBUF - 1, fp1);
                        if ((c = strchr(buffer, '\n')) != NULL)
                            *c = '\0';
                        c = buffer;
                        *(c + 15) = '\0';
                        MarketTags[i] = _strdup(buffer);
                    }

                    for (i = 0; i < NBins; i++) {
                        fgets(buffer, MAXBUF - 1, fp1);
                        tok = strtok(buffer, " \t\r\n");
                        for (j = 0; j < NMarket; j++) {
                            MarketProb[i][j] = atof(tok);
                            tok = strtok(NULL, " \t\r\n");
                        }
                    }

                } else if (strstr(buffer, "SAMPLE")) {
                    for (i = 0; i < NSurvey + 1; i++) {
                        fgets(buffer, MAXBUF - 1, fp1);
                        tok = strtok(buffer, " \t\r\n");
                        MisAgeFlag[i] = atoi(tok);
                        tok = strtok(NULL, " \t\r\n");
                        SampFlag[i] = atoi(tok);
                        tok = strtok(NULL, " \t\r\n");
                        SampFrac[i] = atof(tok);
                    }


                    for (j = 0; j < MaxYears; j++) {
                        fgets(buffer, MAXBUF - 1, fp1);
                        tok = strtok(buffer, " \t\r\n");
                        for (i = 0; i < NMarket; i++) {
                            MarketSampleRate[i][j] = atof(tok);
                            tok = strtok(NULL, " \t\r\n");
                        }
                    }
                } else if (strstr(buffer, "MISAGE")) {
                    for (i = 0; i < NAgeErr; i++) {
                        fgets(buffer, MAXBUF - 1, fp1);
                        tok = strtok(buffer, " \t\r\n");
                        AgeErrStart[i] = atoi(tok);
                        tok = strtok(NULL, " \t\r\n");
                        AgeErrEnd[i] = atoi(tok);
                    }

                    for (k = 0; k < NAgeErr; k++) {
                        for (i = 0; i < NAges; i++) {
                            n = k * NAges + i;
                            fgets(buffer, MAXBUF - 1, fp1);
                            tok = strtok(buffer, " \t\r\n");
                            for (j = 0; j < NAges; j++) {
                                AgeErr[n][j] = atof(tok);
                                tok = strtok(NULL, " \t\r\n");
                            }
                        }
                    }

                } else if (strstr(buffer, "SSELECT")) {
                    for (k = 0; k < NSurvey; k++) {
                        for (i = 0; i < NBins; i++) {
                            n = k * NBins + i;
                            fgets(buffer, MAXBUF - 1, fp1);
                            tok = strtok(buffer, " \t\r\n");
                            for (j = 0; j < MaxYears + 1; j++) {
                                SurveySelect[n][j] = atof(tok);
                                tok = strtok(NULL, " \t\r\n");
                            }
                        }
                    }
                } else if (strstr(buffer, "SURVEY")) {
                    for (i = 0; i < NSurvey; i++) {
                        fgets(buffer, MAXBUF - 1, fp1);
                        tok = strtok(buffer, " \t\r\n");
                        for (j = 0; j < MaxYears + 1; j++) {
                            SurveyQ[i][j] = atof(tok);
                            tok = strtok(NULL, " \t\r\n");
                        }
                    }
                } else if (strstr(buffer, "INDEX")) {
                    for (i = 0; i < NIndex; i++) {
                        fgets(buffer, MAXBUF - 1, fp1);
                        tok = strtok(buffer, " \t\r\n");
                        strcpy(IndexData[i].Name, tok);

                        tok = strtok(NULL, " \t\r\n");
                        IndexData[i].StartAge = atoi(tok);

                        tok = strtok(NULL, " \t\r\n");
                        IndexData[i].EndAge = atoi(tok);

                        tok = strtok(NULL, " \t\r\n");
                        IndexData[i].Type = atoi(tok);

                        tok = strtok(NULL, " \t\r\n");
                        IndexData[i].SurveyIndex = atoi(tok);

                    }
                } else if (strstr(buffer, "SPROD")) {
                    fgets(buffer, MAXBUF - 1, fp1);
                    tok = strtok(buffer, " \t\r\n");
                    SProdDelta = atof(tok);
                } else if (strstr(buffer, "MSR")) {

                    fgets(buffer, MAXBUF - 1, fp1);
                    tok = strtok(buffer, " \t\r\n");
                    MRule = atoi(tok);

                    fgets(buffer, MAXBUF - 1, fp1);
                    tok = strtok(buffer, " \t\r\n");
                    Fmsy = atof(tok);
                    tok = strtok(NULL, " \t\r\n");
                    Bmsy = atof(tok);
                    tok = strtok(NULL, " \t\r\n");
                    Blim = atof(tok);
                    tok = strtok(NULL, " \t\r\n");
                    DeltaLand = atof(tok);

                    fgets(buffer, MAXBUF - 1, fp1);
                    tok = strtok(buffer, " \t\r\n");
                    CvRule = atof(tok);
                    tok = strtok(NULL, " \t\r\n");
                    BiasRule = atof(tok);

                    fgets(buffer, MAXBUF - 1, fp1);
                    tok = strtok(buffer, " \t\r\n");
                    CvSSB = atof(tok);
                    tok = strtok(NULL, " \t\r\n");
                    BiasSSB = atof(tok);

                    fgets(buffer, MAXBUF - 1, fp1);
                    tok = strtok(buffer, " \t\r\n");
                    CvFProj = atof(tok);
                    tok = strtok(NULL, " \t\r\n");
                    BiasFProj = atof(tok);
                }
            }
        }

        void AllocData() {

            short ka, kc, kc1, kx, ky1, kmb, ks1, ks2, ns1, nx, nxa;
            long kb;


            NYears = NXYear - NFYear + 1;
            NAges = MaxAge;
            KAges = NAgePlus;
            NBins = MaxLen - MinLen + 1;
            MaxYears = MXYear - NFYear + 1;

            ky1 = MaxYears + 1;
            kx = NBins * NAges;
            ka = NAges + 1;
            kb = NBins * NAges * ky1;
            kc = NBins * MaxYears;
            kc1 = NBins * ky1;
            kmb = NBins * NMarket;
            ks1 = NSurvey + 1;
            ks2 = NSurvey * NBins;
            ns1 = NSurvey * ky1;
            nxa = NAgeErr * NAges;
            nx = (short) NIter;

            InitMale = (double *) calloc(NAges, sizeof (double));
            InitFemale = (double *) calloc(NAges, sizeof (double));
            InitStock = (double *) calloc(NAges, sizeof (double));
            GrowthAvgMale = (double *) calloc(NAges, sizeof (double));
            GrowthAvgFemale = (double *) calloc(NAges, sizeof (double));
            GrowthIntMale = (double *) calloc(NAges, sizeof (double));
            GrowthIntFemale = (double *) calloc(NAges, sizeof (double));
            GrowthIntFemale = (double *) calloc(NAges, sizeof (double));

            LenStdDevMale = (double *) calloc(NAges, sizeof (double));
            LenStdDevFemale = (double *) calloc(NAges, sizeof (double));

            SampFrac = (double *) calloc(ks1, sizeof (double));
            Recruits = (double *) calloc(ky1, sizeof (double));
            InitRecruits = (double *) calloc(MaxYears, sizeof (double));
            RecSexFemale = (double *) calloc(MaxYears, sizeof (double));
            FracMatMale = (double *) calloc(NBins, sizeof (double));
            FracMatFemale = (double *) calloc(NBins, sizeof (double));
            RecDev = (double *) calloc(MaxYears, sizeof (double));
            FFull = (double *) calloc(MaxYears, sizeof (double));
            SSB = (double *) calloc(MaxYears, sizeof (double));
            StockBiomass = (double *) calloc(ky1, sizeof (double));
            CatchBiomass = (double *) calloc(MaxYears, sizeof (double));
            SurplusProd = (double *) calloc(MaxYears, sizeof (double));

            TotalLand = (double *) calloc(MaxYears, sizeof (double));
            TotalDisc = (double *) calloc(MaxYears, sizeof (double));

            WorkVec = (double *) calloc(NIter, sizeof (double));

            MarketTags = (char **) calloc(NMarket, sizeof (char *));

            LandDev = (double *) calloc(NMarket*MaxYears, sizeof (double));

            SampleLength = (double *) calloc(NBins, sizeof (double));

            SelectLength = (double *) calloc(NBins, sizeof (double));

            SurveyStock = (double *) calloc(NBins, sizeof (double));
            SurveySamples = (double *) calloc(NBins, sizeof (double));
            SurveyDist = (double *) calloc(NBins, sizeof (double));

            MisAgeFlag = (short *) calloc(ks1, sizeof (short));
            SampFlag = (short *) calloc(ks1, sizeof (short));

            AgeErrStart = (short *) calloc(NAgeErr, sizeof (short));
            AgeErrEnd = (short *) calloc(NAgeErr, sizeof (short));



            KSamples = (long *) calloc(ns1, sizeof (long));

            IndexData = (struct si *) calloc(NIndex, sizeof (struct si));

            PopDev = (double *) calloc(NIter, sizeof (double));
            RuleDev = (double *) calloc(MaxYears, sizeof (double));
            SSBDev = (double *) calloc(MaxYears, sizeof (double));
            FProjDev = (double *) calloc(MaxYears, sizeof (double));

            DFrac = (double *) calloc(KAges, sizeof (double));


            AgeProbMale = AllocMatrix(NBins, NAges);
            AgeProbFemale = AllocMatrix(NBins, NAges);

            MaleBins = AllocMatrix(NBins, NAges);
            FemaleBins = AllocMatrix(NBins, NAges);
            StockBins = AllocMatrix(NBins, NAges);

            StockLength = AllocMatrix(NBins, ky1);
            StockAge = AllocMatrix(KAges, ky1);
            StockAgeMax = AllocMatrix(NAges, ky1);
            StockWeights = AllocMatrix(KAges, ky1);

            MaleLength = AllocMatrix(NBins, ky1);
            MaleAge = AllocMatrix(KAges, ky1);
            MaleWeights = AllocMatrix(KAges, ky1);

            FemaleLength = AllocMatrix(NBins, ky1);
            FemaleAge = AllocMatrix(KAges, ky1);
            FemaleWeights = AllocMatrix(KAges, ky1);

            MaleGrowth = AllocMatrix(kx, NBins);
            FemaleGrowth = AllocMatrix(kx, NBins);

            LenWtCoeff = AllocMatrix(ky1, XXLW);
            DiscFrac = AllocMatrix(MaxYears, XXDF);

            NMort = AllocMatrix(MaxYears, XXNM);
            NMSelect = AllocMatrix(NAges, XXNM);
            FSelect = AllocMatrix(NBins, MaxYears);

            MaleCatch = AllocMatrix(NBins, NAges);
            MaleDisc = AllocMatrix(NBins, NAges);
            MaleSurvive = AllocMatrix(NBins, NAges);
            MaleSpStock = AllocMatrix(NBins, NAges);

            FemaleCatch = AllocMatrix(NBins, NAges);
            FemaleDisc = AllocMatrix(NBins, NAges);
            FemaleSurvive = AllocMatrix(NBins, NAges);
            FemaleSpStock = AllocMatrix(NBins, NAges);

            CatchBins = AllocMatrix(NBins, NAges);
            DiscBins = AllocMatrix(NBins, NAges);
            SpawnStock = AllocMatrix(NBins, NAges);
            CatchAgeKey = AllocMatrix(kc, NAges);

            CatchLength = AllocMatrix(NBins, MaxYears);
            DiscLength = AllocMatrix(NBins, MaxYears);

            MaleCatchLength = AllocMatrix(NBins, MaxYears);
            FemaleCatchLength = AllocMatrix(NBins, MaxYears);

            CatchAge = AllocMatrix(KAges, MaxYears);
            DiscAge = AllocMatrix(KAges, MaxYears);
            MaleCatchAge = AllocMatrix(KAges, MaxYears);
            FemaleCatchAge = AllocMatrix(KAges, MaxYears);

            SpStockAge = AllocMatrix(KAges, MaxYears);
            MaleSpawnAge = AllocMatrix(KAges, MaxYears);
            FemaleSpawnAge = AllocMatrix(KAges, MaxYears);
            SSBAge = AllocMatrix(KAges, MaxYears);

            CatchWeights = AllocMatrix(KAges, MaxYears);
            MaleCatchWeights = AllocMatrix(KAges, MaxYears);
            FemaleCatchWeights = AllocMatrix(KAges, MaxYears);

            SpStockWeights = AllocMatrix(KAges, MaxYears);
            MaleSpawnWeights = AllocMatrix(KAges, MaxYears);
            FemaleSpawnWeights = AllocMatrix(KAges, MaxYears);

            Mature = AllocMatrix(KAges, MaxYears);

            Survivorx = AllocMatrix(NBins, NAges);

            FCalc = AllocMatrix(KAges, MaxYears);
            MCalc = AllocMatrix(KAges, MaxYears);

            MarketProb = AllocMatrix(NBins, NMarket);
            MarketCatch = AllocMatrix(kmb, MaxYears);
            MarketDist = AllocMatrix(kmb, MaxYears);
            MarketLand = AllocMatrix(NMarket, MaxYears);
            MarketTarg = AllocMatrix(NMarket, MaxYears);
            MarketSampleRate = AllocMatrix(NMarket, MaxYears);
            MarketSampleSize = AllocMatrix(NMarket, MaxYears);
            MisAge = AllocMatrix(NAges, NAges);
            SampleBins = AllocMatrix(NBins, NMarket);
            SampleAges = AllocMatrix(NBins, NAges);
            AgeLengthKey = AllocMatrix(NBins, NAges);
            AgeKeyCalc = AllocMatrix(NBins, NAges);
            AgeErr = AllocMatrix(nxa, NAges);


            ExpandLength = AllocMatrix(NBins, MaxYears);
            ExpandCatch = AllocMatrix(KAges, MaxYears);
            ExpandWeights = AllocMatrix(KAges, MaxYears);

            SurveySelect = AllocMatrix(ks2, ky1);
            SurveyQ = AllocMatrix(NSurvey, ky1);
            StockAgeKey = AllocMatrix(kc1, NAges);
            SurveyLength = AllocMatrix(ks2, ky1);
            SurveyAge = AllocMatrix(ns1, KAges);
            SurveyBins = AllocMatrix(NBins, NAges);
            SelectAge = AllocMatrix(ns1, KAges);
            IndexValues = AllocMatrix(ky1, NIndex);

            WorkData = AllocMatrix(nx, ky1);

            /* VPA */

            PRVec = (double *) calloc(KAges, sizeof (double));
            StockEst = (double *) calloc(KAges, sizeof (double));

            CatchFactors = (double *) calloc(MaxYears, sizeof (double));

            AgeEst = (short *) calloc(KAges, sizeof (short));

            StockVPA = AllocMatrix(KAges, ky1);
            FVPA = AllocMatrix(KAges, MaxYears);
            SSBVPA = AllocMatrix(KAges, MaxYears);

            VPASurveyFactors = AllocMatrix(ky1, NIndex);

            /* AgePro */

            AgeProDev = (double *) calloc(MaxYears, sizeof (double));
        }

        void CreateInitialDistribution() {
            short i, j;
            double xm, tt, xxt, xxf;
            double x1, x2, p1, p2;

            /* Scale Input Numbers at Age to New Total Population */

            for (j = 0; j < NAges; j++) {
                InitMale[j] = InitMale[j] * InitPop / PrevPop;
                InitFemale[j] = InitFemale[j] * InitPop / PrevPop;
            }



            /* Male Distribution */

            ZeroMatrix(AgeProbMale, NBins, NAges);

            for (j = 0; j < NAges; j++) {
                tt = (double) (NFAge + j);
                xm = LInfMale * (1.0 - exp(-KParmMale * (tt - TZeroMale)));
                GrowthAvgMale[j] = xm;
                xxt = 0.0;
                for (i = 0; i < NBins; i++) {
                    x1 = (double) (MinLen + i);
                    x2 = x1 + 1.0;
                    x1 = (x1 - xm) / LenStdDevMale[j];
                    x2 = (x2 - xm) / LenStdDevMale[j];
                    p2 = DNORDF(&x2);
                    p1 = DNORDF(&x1);
                    AgeProbMale[i][j] = p2 - p1;
                    xxt += AgeProbMale[i][j];
                }

                for (i = 0; i < NBins; i++)
                    AgeProbMale[i][j] = AgeProbMale[i][j] / xxt;
            }
            for (j = 0; j < NAges - 1; j++)
                GrowthIntMale[j] = GrowthAvgMale[j + 1] - GrowthAvgMale[j];



            /* Male Stock Distribution */

            for (j = 0; j < NAges; j++) {
                for (i = 0; i < NBins; i++)
                    MaleBins[i][j] = InitMale[j] * AgeProbMale[i][j];
            }


            /* Female Distribution */

            ZeroMatrix(AgeProbFemale, NBins, NAges);

            for (j = 0; j < NAges; j++) {
                tt = (double) (NFAge + j);
                xm = LInfFemale * (1.0 - exp(-KParmFemale * (tt - TZeroFemale)));
                GrowthAvgFemale[j] = xm;
                xxt = 0.0;
                for (i = 0; i < NBins; i++) {
                    x1 = (double) (MinLen + i);
                    x2 = x1 + 1.0;
                    x1 = (x1 - xm) / LenStdDevFemale[j];
                    x2 = (x2 - xm) / LenStdDevFemale[j];
                    p2 = DNORDF(&x2);
                    p1 = DNORDF(&x1);
                    AgeProbFemale[i][j] = p2 - p1;
                    xxt += AgeProbFemale[i][j];
                }

                for (i = 0; i < NBins; i++)
                    AgeProbFemale[i][j] = AgeProbFemale[i][j] / xxt;
            }
            for (j = 0; j < NAges - 1; j++)
                GrowthIntFemale[j] = GrowthAvgFemale[j + 1] - GrowthAvgFemale[j];



            /* Female Stock Distribution */

            for (j = 0; j < NAges; j++) {
                for (i = 0; i < NBins; i++)
                    FemaleBins[i][j] = InitFemale[j] * AgeProbFemale[i][j];
            }

            /* Combined Stock Distribution */


            xm = 0.0;
            for (j = 0; j < NAges; j++) {
                xxt = InitMale[j] + InitFemale[j];

                tt = xxt / InitPop;

                xm += tt;

                if (xxt > 0.0)
                    xxf = InitFemale[j] / xxt;
                else
                    xxf = 0.0;

                for (i = 0; i < NBins; i++)
                    StockBins[i][j] = MaleBins[i][j] + FemaleBins[i][j];


            }




            /* Save Male Stock Numbers */

            for (i = 0; i < NBins; i++) {
                MaleLength[i][0] = 0.0;
                for (j = 0; j < NAges; j++)
                    MaleLength[i][0] += MaleBins[i][j];
            }

            for (j = 0; j < KAges - 1; j++) {
                MaleAge[j][0] = 0.0;
                for (i = 0; i < NBins; i++)
                    MaleAge[j][0] += MaleBins[i][j];
            }


            MaleAge[KAges - 1][0] = 0.0;
            for (j = KAges - 1; j < NAges; j++)
                for (i = 0; i < NBins; i++)
                    MaleAge[KAges - 1][0] += MaleBins[i][j];



            /* Save Female Stock Numbers */

            for (i = 0; i < NBins; i++) {
                FemaleLength[i][0] = 0.0;
                for (j = 0; j < NAges; j++)
                    FemaleLength[i][0] += FemaleBins[i][j];
            }

            for (j = 0; j < KAges - 1; j++) {
                FemaleAge[j][0] = 0.0;
                for (i = 0; i < NBins; i++)
                    FemaleAge[j][0] += FemaleBins[i][j];
            }


            FemaleAge[KAges - 1][0] = 0.0;
            for (j = KAges - 1; j < NAges; j++)
                for (i = 0; i < NBins; i++)
                    FemaleAge[KAges - 1][0] += FemaleBins[i][j];


            /* Save Combined Stock Numbers */

            for (i = 0; i < NBins; i++) {
                StockLength[i][0] = 0.0;
                for (j = 0; j < NAges; j++)
                    StockLength[i][0] += StockBins[i][j];
            }

            for (j = 0; j < KAges - 1; j++) {
                StockAge[j][0] = 0.0;
                for (i = 0; i < NBins; i++)
                    StockAge[j][0] += StockBins[i][j];
            }


            StockAge[KAges - 1][0] = 0.0;
            for (j = KAges - 1; j < NAges; j++)
                for (i = 0; i < NBins; i++)
                    StockAge[KAges - 1][0] += StockBins[i][j];


            for (j = 0; j < NAges; j++) {
                StockAgeMax[j][0] = 0.0;
                for (i = 0; i < NBins; i++)
                    StockAgeMax[j][0] += StockBins[i][j];
            }


            Recruits[0] = 0.0;

            for (i = 0; i < NBins; i++)
                Recruits[0] += StockBins[i][0];


        }

        void CreateGrowthMatrix() {
            short i, j, k;
            double x1, x2, xx1, xx2, z1, z2, p, px;

            /* Male Growth Projection Matrix */

            for (k = 0; k < NAges; k++) {
                for (i = 0; i < NBins; i++)
                    for (j = 0; j < NBins; j++)
                        MaleGrowth[k * NBins + i][j] = 0.0;
            }

            for (k = 0; k < NAges; k++) {
                for (j = 0; j < NBins; j++) {
                    xx1 = MinLen + j;
                    xx2 = xx1 + 1.0;
                    z1 = (xx1 + xx2) / 2.;

                    z2 = z1 + GrowthIntMale[k];

                    px = 0.0;
                    for (i = j; i < NBins; i++) {
                        xx1 = MinLen + i;
                        xx2 = xx1 + 1.0;
                        x1 = (xx1 - z2) / GrowthStdDevMale;
                        x2 = (xx2 - z2) / GrowthStdDevMale;
                        p = DNORDF(&x2) - DNORDF(&x1);
                        MaleGrowth[k * NBins + i][j] = p;
                        px += p;
                    }
                    if (px < 1.00000000)
                        for (i = j; i < NBins; i++)
                            MaleGrowth[k * NBins + i][j] = MaleGrowth[k * NBins + i][j] / px;
                }
            }

            /* Female Growth Projection Matrix */

            for (k = 0; k < NAges; k++) {
                for (i = 0; i < NBins; i++)
                    for (j = 0; j < NBins; j++)
                        FemaleGrowth[k * NBins + i][j] = 0.0;
            }

            for (k = 0; k < NAges; k++) {
                for (j = 0; j < NBins; j++) {
                    xx1 = MinLen + j;
                    xx2 = xx1 + 1.0;
                    z1 = (xx1 + xx2) / 2.;

                    z2 = z1 + GrowthIntFemale[k];

                    px = 0.0;
                    for (i = j; i < NBins; i++) {
                        xx1 = MinLen + i;
                        xx2 = xx1 + 1.0;
                        x1 = (xx1 - z2) / GrowthStdDevFemale;
                        x2 = (xx2 - z2) / GrowthStdDevFemale;
                        p = DNORDF(&x2) - DNORDF(&x1);
                        FemaleGrowth[k * NBins + i][j] = p;
                        px += p;
                    }
                    if (px < 1.00000000)

                        for (i = j; i < NBins; i++)
                            FemaleGrowth[k * NBins + i][j] = FemaleGrowth[k * NBins + i][j] / px;
                }
            }


        }

        void CalcMaturityFraction() {
            short i;
            double xl;
            double a, b;
            double xp;

            a = FMAlpha;
            b = FMBeta;

            for (i = 0; i < NBins; i++) {
                xl = (double) (MinLen + i);
                xp = exp(-b * (xl - a));
                FracMatMale[i] = 1.0 / (1.0 + xp);
            }

            a = FFAlpha;
            b = FFBeta;

            for (i = 0; i < NBins; i++) {

                xl = (double) (MinLen + i);
                xp = exp(-b * (xl - a));
                FracMatFemale[i] = 1.0 / (1.0 + xp);
            }


        }

        void CalcStockWeights(short k) {
            if (!(k % this->assessment_frequency) == 0) {
                                                return;
                    }
            short i, j;
            double xl, xt, xw, xx;
            double alpha, beta;


            /* Males */

            alpha = LenWtCoeff[k][0];
            beta = LenWtCoeff[k][1];


            for (j = 0; j < KAges - 1; j++) {
                xt = 0.0;
                xx = 0.0;
                for (i = 0; i < NBins; i++) {
                    xl = (double) (MinLen + i);
                    xw = exp(alpha + beta * log(xl));
                    xt += MaleBins[i][j];
                    xx += xw * MaleBins[i][j];
                }
                if (xt > 0.0)
                    MaleWeights[j][k] = xx / xt;
                else
                    MaleWeights[j][k] = 0.0;
            }

            xx = 0.0;
            xt = 0.0;
            for (j = KAges - 1; j < NAges; j++) {
                for (i = 0; i < NBins; i++) {
                    xl = (double) (MinLen + i);
                    xw = exp(alpha + beta * log(xl));
                    xt += MaleBins[i][j];
                    xx += xw * MaleBins[i][j];
                }
            }
            if (xt > 0.0)
                MaleWeights[KAges - 1][k] = xx / xt;
            else
                MaleWeights[KAges - 1][k] = 0.0;

            /* Females */

            alpha = LenWtCoeff[k][2];
            beta = LenWtCoeff[k][3];


            for (j = 0; j < KAges - 1; j++) {
                xt = 0.0;
                xx = 0.0;
                for (i = 0; i < NBins; i++) {
                    xl = (double) (MinLen + i);
                    xw = exp(alpha + beta * log(xl));
                    xt += FemaleBins[i][j];
                    xx += xw * FemaleBins[i][j];
                }
                if (xt > 0.0)
                    FemaleWeights[j][k] = xx / xt;
                else
                    FemaleWeights[j][k] = 0.0;
            }

            xx = 0.0;
            xt = 0.0;
            for (j = KAges - 1; j < NAges; j++) {
                for (i = 0; i < NBins; i++) {
                    xl = (double) (MinLen + i);
                    xw = exp(alpha + beta * log(xl));
                    xt += FemaleBins[i][j];
                    xx += xw * FemaleBins[i][j];
                }
            }
            if (xt > 0.0)
                FemaleWeights[KAges - 1][k] = xx / xt;
            else
                FemaleWeights[KAges - 1][k] = 0.0;

            /* Combined */

            StockBiomass[k] = 0.0;
            for (j = 0; j < KAges; j++) {
                if (StockAge[j][k] > 0.0)
                    StockWeights[j][k] = (MaleAge[j][k] * MaleWeights[j][k] + FemaleAge[j][k] * FemaleWeights[j][k]) / StockAge[j][k];

                else
                    StockWeights[j][k] = 0.0;

                StockBiomass[k] += MaleAge[j][k] * MaleWeights[j][k] + FemaleAge[j][k] * FemaleWeights[j][k];
            }


        }

        double DiscardFraction(double xl, double a, double b) {
            double xp, df;


            if (!DiscFlag)
                return 0.0;

            if (a < XTOL && b < XTOL)
                return 0.0;


            xp = exp(-b * (xl - a));
            df = 1.0 / (1.0 + xp);

            return df;

        }

        void CalculateCatch(short k) {
            short i, j, n;
            double f, m, z, fd, zx;
            double xc, xl, xs, xd, xt;
            double ad, bd;

            ad = DiscFrac[k][0];
            bd = DiscFrac[k][1];

            /* Apply Mortality to Males */


            for (i = 0; i < NBins; i++) {
                xl = (double) (MinLen + i);

                fd = DiscardFraction(xl, ad, bd);
                xd = 1.0 - fd;

                f = FFull[k] * FSelect[i][k];


                for (j = 0; j < NAges; j++) {

                    m = NMort[k][0] * NMSelect[j][0];

                    z = f + m;

                    xc = f * (1.0 - exp(-z)) / z;


                    MaleCatch[i][j] = MaleBins[i][j] * xc * xd;
                    MaleDisc[i][j] = MaleBins[i][j] * xc * fd;

                    xs = exp(-z);

                    MaleSurvive[i][j] = MaleBins[i][j] * xs;

                    if (SSBFlag == 0) {
                        zx = f * TF + m * TM;

                        MaleSpStock[i][j] = MaleBins[i][j] * FracMatMale[i] * exp(-zx);
                    } else
                        MaleSpStock[i][j] = 0.0;
                }
            }



            /* Apply Mortality to Females */

            for (i = 0; i < NBins; i++) {
                xl = (double) (MinLen + i);

                fd = DiscardFraction(xl, ad, bd);
                xd = 1.0 - fd;

                f = FFull[k] * FSelect[i][k];


                for (j = 0; j < NAges; j++) {
                    m = NMort[k][1] * NMSelect[j][1];

                    z = f + m;

                    xc = f * (1.0 - exp(-z)) / z;

                    FemaleCatch[i][j] = FemaleBins[i][j] * xc * xd;
                    FemaleDisc[i][j] = FemaleBins[i][j] * xc * fd;

                    xs = exp(-z);

                    FemaleSurvive[i][j] = FemaleBins[i][j] * xs;

                    zx = f * TF + m * TM;

                    FemaleSpStock[i][j] = FemaleBins[i][j] * FracMatFemale[i] * exp(-zx);
                }
            }

            /* Combine Males & Females */

            for (i = 0; i < NBins; i++) {
                for (j = 0; j < NAges; j++) {
                    CatchBins[i][j] = MaleCatch[i][j] + FemaleCatch[i][j];
                    DiscBins[i][j] = MaleDisc[i][j] + FemaleDisc[i][j];
                    SpawnStock[i][j] = MaleSpStock[i][j] + FemaleSpStock[i][j];
                }
            }

            /* Catch Age Key */

            for (i = 0; i < NBins; i++) {
                n = k * NBins + i;

                xt = 0.0;

                for (j = 0; j < NAges; j++)
                    xt += CatchBins[i][j];

                for (j = 0; j < NAges; j++) {
                    if (xt > 0.0)
                        CatchAgeKey[n][j] = CatchBins[i][j] / xt;
                    else
                        CatchAgeKey[n][j] = 0.0;
                }
            }

            for (i = 0; i < NBins; i++) {
                CatchLength[i][k] = 0.0;
                for (j = 0; j < NAges; j++)
                    CatchLength[i][k] += CatchBins[i][j];
            }

            for (i = 0; i < NBins; i++) {
                MaleCatchLength[i][k] = 0.0;
                for (j = 0; j < NAges; j++)
                    MaleCatchLength[i][k] += MaleCatch[i][j];
            }

            for (i = 0; i < NBins; i++) {
                FemaleCatchLength[i][k] = 0.0;
                for (j = 0; j < NAges; j++)
                    FemaleCatchLength[i][k] += FemaleCatch[i][j];
            }

            for (i = 0; i < NBins; i++) {
                DiscLength[i][k] = 0.0;
                for (j = 0; j < NAges; j++)
                    DiscLength[i][k] += DiscBins[i][j];
            }


            /* Catch & Discards at Age */

            for (j = 0; j < KAges - 1; j++) {
                CatchAge[j][k] = 0.0;
                for (i = 0; i < NBins; i++)
                    CatchAge[j][k] += CatchBins[i][j];
            }


            CatchAge[KAges - 1][k] = 0.0;
            for (j = KAges - 1; j < NAges; j++)
                for (i = 0; i < NBins; i++)
                    CatchAge[KAges - 1][k] += CatchBins[i][j];


            for (j = 0; j < KAges - 1; j++) {
                DiscAge[j][k] = 0.0;
                for (i = 0; i < NBins; i++)
                    DiscAge[j][k] += DiscBins[i][j];
            }


            DiscAge[KAges - 1][k] = 0.0;
            for (j = KAges - 1; j < NAges; j++)
                for (i = 0; i < NBins; i++)
                    DiscAge[KAges - 1][k] += DiscBins[i][j];

            /* Male Catch at Age */

            for (j = 0; j < KAges - 1; j++) {
                MaleCatchAge[j][k] = 0.0;
                for (i = 0; i < NBins; i++)
                    MaleCatchAge[j][k] += MaleCatch[i][j];
            }


            MaleCatchAge[KAges - 1][k] = 0.0;
            for (j = KAges - 1; j < NAges; j++)
                for (i = 0; i < NBins; i++)
                    MaleCatchAge[KAges - 1][k] += MaleCatch[i][j];


            /* Female Catch at Age */

            for (j = 0; j < KAges - 1; j++) {
                FemaleCatchAge[j][k] = 0.0;
                for (i = 0; i < NBins; i++)
                    FemaleCatchAge[j][k] += FemaleCatch[i][j];
            }


            FemaleCatchAge[KAges - 1][k] = 0.0;
            for (j = KAges - 1; j < NAges; j++)
                for (i = 0; i < NBins; i++)
                    FemaleCatchAge[KAges - 1][k] += FemaleCatch[i][j];

            /* Spawning Stock at Age */

            for (j = 0; j < KAges - 1; j++) {
                SpStockAge[j][k] = 0.0;
                for (i = 0; i < NBins; i++)
                    SpStockAge[j][k] += SpawnStock[i][j];
            }


            SpStockAge[KAges - 1][k] = 0.0;
            for (j = KAges - 1; j < NAges; j++)
                for (i = 0; i < NBins; i++)
                    SpStockAge[KAges - 1][k] += SpawnStock[i][j];



            /* Male Spawning Stock at Age */

            for (j = 0; j < KAges - 1; j++) {
                MaleSpawnAge[j][k] = 0.0;
                for (i = 0; i < NBins; i++)
                    MaleSpawnAge[j][k] += MaleSpStock[i][j];
            }


            MaleSpawnAge[KAges - 1][k] = 0.0;
            for (j = KAges - 1; j < NAges; j++)
                for (i = 0; i < NBins; i++)
                    MaleSpawnAge[KAges - 1][k] += MaleSpStock[i][j];


            /* Female Spawning Stock at Age */

            for (j = 0; j < KAges - 1; j++) {
                FemaleSpawnAge[j][k] = 0.0;
                for (i = 0; i < NBins; i++)
                    FemaleSpawnAge[j][k] += FemaleSpStock[i][j];
            }


            FemaleSpawnAge[KAges - 1][k] = 0.0;
            for (j = KAges - 1; j < NAges; j++)

                for (i = 0; i < NBins; i++)
                    FemaleSpawnAge[KAges - 1][k] += FemaleSpStock[i][j];

        }

        void CalcCatchWeights(short k) {
            short i, j;
            double xl, xt, xw, xx;
            double alpha, beta;


            /* Males */

            alpha = LenWtCoeff[k][0];
            beta = LenWtCoeff[k][1];


            for (j = 0; j < KAges - 1; j++) {
                xt = 0.0;
                xx = 0.0;
                for (i = 0; i < NBins; i++) {
                    xl = (double) (MinLen + i);
                    xw = exp(alpha + beta * log(xl));
                    xt += MaleCatch[i][j];
                    xx += xw * MaleCatch[i][j];
                }
                if (xt > 0.0)
                    MaleCatchWeights[j][k] = xx / xt;
                else
                    MaleCatchWeights[j][k] = 0.0;
            }

            xx = 0.0;
            xt = 0.0;
            for (j = KAges - 1; j < NAges; j++) {
                for (i = 0; i < NBins; i++) {
                    xl = (double) (MinLen + i);
                    xw = exp(alpha + beta * log(xl));
                    xt += MaleCatch[i][j];
                    xx += xw * MaleCatch[i][j];
                }
            }
            if (xt > 0.0)
                MaleCatchWeights[KAges - 1][k] = xx / xt;
            else
                MaleCatchWeights[KAges - 1][k] = 0.0;

            /* Females */

            alpha = LenWtCoeff[k][2];
            beta = LenWtCoeff[k][3];


            for (j = 0; j < KAges - 1; j++) {
                xt = 0.0;
                xx = 0.0;
                for (i = 0; i < NBins; i++) {
                    xl = (double) (MinLen + i);
                    xw = exp(alpha + beta * log(xl));
                    xt += FemaleCatch[i][j];
                    xx += xw * FemaleCatch[i][j];
                }
                if (xt > 0.0)
                    FemaleCatchWeights[j][k] = xx / xt;
                else
                    FemaleCatchWeights[j][k] = 0.0;
            }

            xx = 0.0;
            xt = 0.0;
            for (j = KAges - 1; j < NAges; j++) {
                for (i = 0; i < NBins; i++) {
                    xl = (double) (MinLen + i);
                    xw = exp(alpha + beta * log(xl));
                    xt += FemaleCatch[i][j];
                    xx += xw * FemaleCatch[i][j];
                }
            }
            if (xt > 0.0)
                FemaleCatchWeights[KAges - 1][k] = xx / xt;
            else
                FemaleCatchWeights[KAges - 1][k] = 0.0;

            /* Combined */


            CatchBiomass[k] = 0.0;
            for (j = 0; j < KAges; j++) {
                if (CatchAge[j][k] > 0.0)
                    CatchWeights[j][k] = (MaleCatchAge[j][k] * MaleCatchWeights[j][k] + FemaleCatchAge[j][k] * FemaleCatchWeights[j][k]) / CatchAge[j][k];

                else
                    CatchWeights[j][k] = 0.0;

                CatchBiomass[k] += MaleCatchAge[j][k] * MaleCatchWeights[j][k] + FemaleCatchAge[j][k] * FemaleCatchWeights[j][k];
            }


        }

        void CalcSpStockWeights(short k) {
            short i, j;
            double xl, xt, xw, xx;
            double alpha, beta;


            /* Males */

            alpha = LenWtCoeff[k][0];
            beta = LenWtCoeff[k][1];


            for (j = 0; j < KAges - 1; j++) {
                xt = 0.0;
                xx = 0.0;
                for (i = 0; i < NBins; i++) {
                    xl = (double) (MinLen + i);
                    xw = exp(alpha + beta * log(xl));
                    xt += MaleSpStock[i][j];
                    xx += xw * MaleSpStock[i][j];
                }
                if (xt > 0.0)
                    MaleSpawnWeights[j][k] = xx / xt;
                else
                    MaleSpawnWeights[j][k] = 0.0;
            }

            xx = 0.0;
            xt = 0.0;
            for (j = KAges - 1; j < NAges; j++) {
                for (i = 0; i < NBins; i++) {
                    xl = (double) (MinLen + i);
                    xw = exp(alpha + beta * log(xl));
                    xt += MaleSpStock[i][j];
                    xx += xw * MaleSpStock[i][j];
                }
            }
            if (xt > 0.0)
                MaleSpawnWeights[KAges - 1][k] = xx / xt;
            else
                MaleSpawnWeights[KAges - 1][k] = 0.0;

            /* Females */

            alpha = LenWtCoeff[k][2];
            beta = LenWtCoeff[k][3];


            for (j = 0; j < KAges - 1; j++) {
                xt = 0.0;
                xx = 0.0;
                for (i = 0; i < NBins; i++) {
                    xl = (double) (MinLen + i);
                    xw = exp(alpha + beta * log(xl));
                    xt += FemaleSpStock[i][j];
                    xx += xw * FemaleSpStock[i][j];
                }
                if (xt > 0.0)
                    FemaleSpawnWeights[j][k] = xx / xt;
                else
                    FemaleSpawnWeights[j][k] = 0.0;
            }

            xx = 0.0;
            xt = 0.0;
            for (j = KAges - 1; j < NAges; j++) {
                for (i = 0; i < NBins; i++) {
                    xl = (double) (MinLen + i);
                    xw = exp(alpha + beta * log(xl));
                    xt += FemaleSpStock[i][j];
                    xx += xw * FemaleSpStock[i][j];
                }
            }
            if (xt > 0.0)
                FemaleSpawnWeights[KAges - 1][k] = xx / xt;
            else
                FemaleSpawnWeights[KAges - 1][k] = 0.0;

            /* Combined */

            SSB[k] = 0.0;
            for (j = 0; j < KAges; j++) {
                if (CatchAge[j][k] > 0.0)
                    SpStockWeights[j][k] = (MaleSpawnAge[j][k] * MaleSpawnWeights[j][k] + FemaleSpawnAge[j][k] * FemaleSpawnWeights[j][k]) / SpStockAge[j][k];

                else
                    SpStockWeights[j][k] = 0.0;

                SSBAge[j][k] = MaleSpawnAge[j][k] * MaleSpawnWeights[j][k] + FemaleSpawnAge[j][k] * FemaleSpawnWeights[j][k];

                SSB[k] += MaleSpawnAge[j][k] * MaleSpawnWeights[j][k] + FemaleSpawnAge[j][k] * FemaleSpawnWeights[j][k];
            }


        }

        void CalcMaturity(short k) {
            short i, j;
            double xn, xx;

            if (SSBFlag == 0) {

                for (j = 0; j < KAges - 1; j++) {
                    xn = 0.0;
                    xx = 0.0;
                    for (i = 0; i < NBins; i++) {
                        xn += MaleBins[i][j];
                        xx += MaleBins[i][j] * FracMatMale[i];

                        xn += FemaleBins[i][j];
                        xx += FemaleBins[i][j] * FracMatFemale[i];
                    }

                    Mature[j][k] = xx / xn;
                }


                xn = 0.0;
                xx = 0.0;
                for (j = KAges - 1; j < NAges; j++) {
                    for (i = 0; i < NBins; i++) {
                        xn += MaleBins[i][j];
                        xx += MaleBins[i][j] * FracMatMale[i];

                        xn += FemaleBins[i][j];
                        xx += FemaleBins[i][j] * FracMatFemale[i];
                    }

                    Mature[KAges - 1][k] = xx / xn;
                }
            } else {
                for (j = 0; j < KAges - 1; j++) {
                    xn = 0.0;
                    xx = 0.0;
                    for (i = 0; i < NBins; i++) {
                        xn += FemaleBins[i][j];
                        xx += FemaleBins[i][j] * FracMatFemale[i];
                    }

                    Mature[j][k] = xx / xn;
                }


                xn = 0.0;
                xx = 0.0;
                for (j = KAges - 1; j < NAges; j++) {
                    for (i = 0; i < NBins; i++) {

                        xn += FemaleBins[i][j];
                        xx += FemaleBins[i][j] * FracMatFemale[i];
                    }

                    Mature[KAges - 1][k] = xx / xn;
                }
            }
        }

        void ApplyGrowthMatrix() {
            short i, j, k;

            /* Males */

            for (k = 0; k < NAges; k++) {
                for (i = 0; i < NBins; i++) {
                    Survivorx[i][k] = 0.0;
                    for (j = 0; j < NBins; j++)
                        Survivorx[i][k] += MaleGrowth[k * NBins + i][j] * MaleSurvive[j][k];
                }
            }


            /* Advance Stock Matrix */

            ZeroMatrix(MaleBins, NBins, NAges);

            for (j = 0; j < NAges - 1; j++) {
                for (i = 0; i < NBins; i++)
                    MaleBins[i][j + 1] = Survivorx[i][j];
            }

            /* Females */


            for (k = 0; k < NAges; k++) {
                for (i = 0; i < NBins; i++) {
                    Survivorx[i][k] = 0.0;
                    for (j = 0; j < NBins; j++)
                        Survivorx[i][k] += FemaleGrowth[k * NBins + i][j] * FemaleSurvive[j][k];
                }
            }

            /* Advance Stock Matrix */

            ZeroMatrix(FemaleBins, NBins, NAges);

            for (j = 0; j < NAges - 1; j++) {
                for (i = 0; i < NBins; i++)
                    FemaleBins[i][j + 1] = Survivorx[i][j];
            }


            /* Combined */


            for (j = 0; j < NAges; j++) {

                for (i = 0; i < NBins; i++)
                    StockBins[i][j] = MaleBins[i][j] + FemaleBins[i][j];
            }

        }

        void CalcRecruits(short k) {
            short i;
            double xs, xrec;
            double zx, sd;
            double xx, xxx;

            xs = SSB[k];
            zx = RecDev[k];
            sd = sqrt(log(CVRecruit * CVRecruit + 1.0));

            switch (RecruitFlag) {
                case 0:
                    xrec = InitRecruits[k];
                    break;
                case 1:
                    xx = xs / 1.0E+06;
                    xrec = RAlpha * xx / (RBeta + xx);
                    xrec = xrec * 1.0E+06;
                    break;
                case 2:
                    xx = xs / 1.0E+06;
                    xrec = xx * exp(RAlpha + RBeta * xx);
                    xrec = xrec * 1.0E+06;
                    break;
                case 3:
                    xx = xs / 1.0E+06;
                    xxx = xx / RKPar;
                    xxx = pow(xxx, RBeta);
                    xrec = RAlpha * xx / (1.0 + xxx);
                    xrec = xrec * 1.0E+06;
                    break;
                default:
                    xrec = InitRecruits[k];
            }

            Recruits[k + 1] = xrec * exp(sd * zx);

            /* Insert Recruits */

            for (i = 0; i < NBins; i++) {

                MaleBins[i][0] = Recruits[k + 1] * (1.0 - RecSexFemale[k]) * AgeProbMale[i][0];
                FemaleBins[i][0] = Recruits[k + 1] * RecSexFemale[k] * AgeProbFemale[i][0];
                StockBins[i][0] = MaleBins[i][0] + FemaleBins[i][0];
            }
        }

        void SaveAnnualData(short k) {
            short i, j;

            /*  Stock Numbers */

            for (i = 0; i < NBins; i++) {
                StockLength[i][k + 1] = 0.0;
                for (j = 0; j < NAges; j++)
                    StockLength[i][k + 1] += StockBins[i][j];
            }

            for (j = 0; j < KAges - 1; j++) {
                StockAge[j][k + 1] = 0.0;
                for (i = 0; i < NBins; i++)
                    StockAge[j][k + 1] += StockBins[i][j];
            }


            StockAge[KAges - 1][k + 1] = 0.0;
            for (j = KAges - 1; j < NAges; j++)
                for (i = 0; i < NBins; i++)
                    StockAge[KAges - 1][k + 1] += StockBins[i][j];


            for (j = 0; j < NAges; j++) {
                StockAgeMax[j][k + 1] = 0.0;
                for (i = 0; i < NBins; i++)
                    StockAgeMax[j][k + 1] += StockBins[i][j];
            }


            /* Males */


            for (j = 0; j < KAges - 1; j++) {
                MaleAge[j][k + 1] = 0.0;
                for (i = 0; i < NBins; i++)
                    MaleAge[j][k + 1] += MaleBins[i][j];
            }


            MaleAge[KAges - 1][k + 1] = 0.0;
            for (j = KAges - 1; j < NAges; j++)
                for (i = 0; i < NBins; i++)
                    MaleAge[KAges - 1][k + 1] += MaleBins[i][j];


            for (i = 0; i < NBins; i++) {
                MaleLength[i][k + 1] = 0.0;
                for (j = 0; j < NAges; j++)
                    MaleLength[i][k + 1] += MaleBins[i][j];
            }



            /* Females */

            for (j = 0; j < KAges - 1; j++) {
                FemaleAge[j][k + 1] = 0.0;
                for (i = 0; i < NBins; i++)
                    FemaleAge[j][k + 1] += FemaleBins[i][j];
            }


            FemaleAge[KAges - 1][k + 1] = 0.0;
            for (j = KAges - 1; j < NAges; j++)
                for (i = 0; i < NBins; i++)
                    FemaleAge[KAges - 1][k + 1] += FemaleBins[i][j];


            for (i = 0; i < NBins; i++) {
                FemaleLength[i][k + 1] = 0.0;

                for (j = 0; j < NAges; j++)
                    FemaleLength[i][k + 1] += FemaleBins[i][j];
            }


        }

        void CalcStockAgeKey(short k) {
            short i, j, n;
            double xt;
            if (!(k % this->assessment_frequency) == 0) {
                                                return;
                    }

            /* Stock Age Length Key */

            for (i = 0; i < NBins; i++) {
                xt = 0.0;
                n = k * NBins + i;
                for (j = 0; j < NAges; j++) {
                    xt += StockBins[i][j];
                }

                for (j = 0; j < NAges; j++) {
                    if (xt > 0.0)
                        StockAgeKey[n][j] = StockBins[i][j] / xt;

                    else
                        StockAgeKey[n][j] = 0.0;
                }
            }
        }

        double SolveCatch(double c, double n, double m) {
            short iter;
            double xx, err1, err2, slope;
            double f, z, cx;


            /* Check for Missing Data */

            if (c < XTOL) {
                f = XLBOUND;
                return f;
            }

            if (n < XTOL) {
                f = XUBOUND;
                return (f);
            }


            /* Initial Guess */

            f = c / n;

            /* Newton Raphson Solver */

            iter = 0;
            while (iter < MXITER) {
                z = m + f;
                cx = n * (1.0 - exp(-z)) * f / z;

                err1 = (cx - c) / c;

                if (fabs(err1) < XTOL)
                    return f;

                if (iter == 0) {
                    xx = f;
                    err2 = err1;
                    if (c > cx)
                        f = f * 1.2;
                    else
                        f = f * 0.8;
                } else {
                    slope = (err2 - err1) / (xx - f);
                    xx = f;
                    err2 = err1;
                    f = (slope * xx - err2) / slope;
                }

                if (f > XUBOUND)
                    f = XUBOUND;

                if (f < XLBOUND)
                    f = XLBOUND;

                iter++;
            }

            return f;
        }

        void BasicReports() {
            short i, j, k;
            double zero = 0.0;
            double xt;


            fprintf(fp2, "Initial Population %18.0f\n\n", InitPop);

            fprintf(fp4, "%18.0f\n", InitPop);

            fprintf(fp2, "Stock at Age\n\n");

            fprintf(fp2, "Year     ");
            for (j = 0; j < NYears + 1; j++)
                fprintf(fp2, "%10d ", NFYear + j);
            fprintf(fp2, "\n\n");

            for (i = 0; i < KAges; i++) {
                fprintf(fp2, "Age %3d  ", NFAge + i);
                for (j = 0; j < NYears + 1; j++)
                    fprintf(fp2, "%10.0f ", StockAge[i][j]);
                fprintf(fp2, "\n");
            }
            fprintf(fp2, "\n");

            fprintf(fp2, "Total    ");
            for (j = 0; j < NYears + 1; j++) {
                xt = 0.0;
                for (i = 0; i < KAges; i++)
                    xt += StockAge[i][j];

                fprintf(fp2, "%10.0f ", xt);
                fprintf(fp4, "%10.0f ", xt);
            }

            fprintf(fp2, "\n\n");
            fprintf(fp4, "\n");

            fprintf(fp2, "Catch at Age\n\n");

            fprintf(fp2, "Year     ");
            for (j = 0; j < NYears; j++)
                fprintf(fp2, "%10d ", NFYear + j);
            fprintf(fp2, "\n\n");

            for (i = 0; i < KAges; i++) {
                fprintf(fp2, "Age %3d  ", NFAge + i);
                for (j = 0; j < NYears; j++)
                    fprintf(fp2, "%10.0f ", CatchAge[i][j]);
                fprintf(fp2, "\n");
            }
            fprintf(fp2, "\n");

            fprintf(fp2, "Total    ");
            for (j = 0; j < NYears; j++) {
                xt = 0.0;
                for (i = 0; i < KAges; i++)
                    xt += CatchAge[i][j];

                fprintf(fp2, "%10.0f ", xt);
                fprintf(fp4, "%10.0f ", xt);
            }

            fprintf(fp2, "\n\n");
            fprintf(fp4, "\n");

            fprintf(fp2, "Stock Weight at Age\n\n");

            fprintf(fp2, "Year     ");
            for (j = 0; j < NYears + 1; j++)
                fprintf(fp2, "%10d ", NFYear + j);
            fprintf(fp2, "\n\n");

            for (i = 0; i < KAges; i++) {
                fprintf(fp2, "Age %3d  ", NFAge + i);
                for (j = 0; j < NYears + 1; j++)
                    fprintf(fp2, "%10.4f ", StockWeights[i][j]);
                fprintf(fp2, "\n");
            }
            fprintf(fp2, "\n\n");



            fprintf(fp2, "Catch Weight at Age\n\n");

            fprintf(fp2, "Year     ");
            for (j = 0; j < NYears; j++)
                fprintf(fp2, "%10d ", NFYear + j);
            fprintf(fp2, "\n\n");

            for (i = 0; i < KAges; i++) {
                fprintf(fp2, "Age %3d  ", NFAge + i);
                for (j = 0; j < NYears; j++)
                    fprintf(fp2, "%10.4f ", CatchWeights[i][j]);
                fprintf(fp2, "\n");
            }
            fprintf(fp2, "\n\n");


            fprintf(fp2, "Spawning Stock Weight at Age\n\n");

            fprintf(fp2, "Year     ");
            for (j = 0; j < NYears; j++)
                fprintf(fp2, "%10d ", NFYear + j);
            fprintf(fp2, "\n\n");

            for (i = 0; i < KAges; i++) {
                fprintf(fp2, "Age %3d  ", NFAge + i);
                for (j = 0; j < NYears; j++)
                    fprintf(fp2, "%10.4f ", SpStockWeights[i][j]);
                fprintf(fp2, "\n");
            }
            fprintf(fp2, "\n\n");


            fprintf(fp2, "Maturity at Age\n\n");

            fprintf(fp2, "Year     ");
            for (j = 0; j < NYears; j++)
                fprintf(fp2, "%10d ", NFYear + j);
            fprintf(fp2, "\n\n");

            for (i = 0; i < KAges; i++) {
                fprintf(fp2, "Age %3d  ", NFAge + i);
                for (j = 0; j < NYears; j++)
                    fprintf(fp2, "%10.4f ", Mature[i][j]);
                fprintf(fp2, "\n");
            }
            fprintf(fp2, "\n\n");

            /* Recruitment & SSB */

            fprintf(fp2, "Recruitment & SSB\n\n");

            fprintf(fp2, "Year             ");
            for (j = 0; j < NYears + 1; j++)
                fprintf(fp2, "%10d ", NFYear + j);
            fprintf(fp2, "\n\n");

            fprintf(fp2, "Spawning_Biomass ");

            for (j = 0; j < NYears; j++) {
                fprintf(fp2, "%10.0f ", SSB[j]);
                fprintf(fp4, "%10.0f ", SSB[j] / 1000.);
            }
            fprintf(fp2, "\n");
            fprintf(fp4, "\n");

            fprintf(fp2, "Recruits_Calc    ");
            for (j = 0; j < NYears + 1; j++) {
                fprintf(fp2, "%10.0f ", Recruits[j]);
                fprintf(fp4, "%10.0f ", Recruits[j] / 1000.);
            }
            fprintf(fp2, "\n\n");
            fprintf(fp4, "\n");

            fprintf(fp2, "Fully Recruited F");
            for (j = 0; j < NYears + 1; j++) {
                fprintf(fp2, "%10.4f ", FFull[j]);
                fprintf(fp4, "%10.4f ", FFull[j]);
            }

            fprintf(fp2, "\n\n");
            fprintf(fp4, "\n");

            /* Landings & Discards (MT) */

            fprintf(fp2, "Total Landings & Discards(MT)\n\n");


            fprintf(fp2, "Year             ");
            for (k = 0; k < NYears; k++)
                fprintf(fp2, "%10d ", NFYear + k);
            fprintf(fp2, "\n\n");
            fprintf(fp2, "Landings         ");
            for (k = 0; k < NYears; k++) {
                fprintf(fp2, "%10.3f ", TotalLand[k]);
                fprintf(fp4, "%10.3f ", TotalLand[k]);
            }
            fprintf(fp2, "\n");
            fprintf(fp4, "\n");

            fprintf(fp2, "Discards         ");
            for (k = 0; k < NYears; k++) {

                fprintf(fp2, "%10.3f ", TotalDisc[k]);
                fprintf(fp4, "%10.3f ", TotalDisc[k]);
            }
            fprintf(fp2, "\n\n");
            fprintf(fp4, "\n\n");

        }

        double CalcAvgCatchWeight(short i, short k) {
            double alpha1, beta1, xw1;
            double alpha2, beta2, xw2;
            double xl, xw;

            alpha1 = LenWtCoeff[k][0];
            beta1 = LenWtCoeff[k][1];

            alpha2 = LenWtCoeff[k][2];
            beta2 = LenWtCoeff[k][3];

            xl = (double) (MinLen + i);

            xw1 = exp(alpha1 + beta1 * log(xl));
            xw2 = exp(alpha2 + beta2 * log(xl));

            if (CatchLength[i][k] > 0.0)
                xw = (MaleCatchLength[i][k] * xw1 + FemaleCatchLength[i][k] * xw2) / CatchLength[i][k];
            else
                xw = (xw1 + xw2) / 2.0;

            return xw;
        }

        void CalcLandings(short k) {
            short i, j, n;
            double xw, xt, xn, xxt;

            /* Calculate Total Landings (MT) */

            TotalLand[k] = 0.0;

            for (j = 0; j < NMarket; j++) {
                xt = 0.0;
                xxt = 0.0;

                for (i = 0; i < NBins; i++) {
                    n = j * NBins + i;

                    xw = CalcAvgCatchWeight(i, k);

                    MarketCatch[n][k] = CatchLength[i][k] * MarketProb[i][j];

                    xt += MarketCatch[n][k] * xw / 1000.0;

                    xxt += MarketCatch[n][k];

                }

                MarketLand[j][k] = xt;
                TotalLand[k] += xt;

                /* Market Cumulative Distribution */

                xn = 0.0;

                for (i = 0; i < NBins; i++) {
                    n = j * NBins + i;

                    xt = MarketCatch[n][k] / xxt;

                    xn += xt;

                    MarketDist[n][k] = xn;
                }

            }

            /* Calculate Total Discards (MT) */

            TotalDisc[k] = 0.0;

            for (i = 0; i < NBins; i++) {

                xw = CalcAvgCatchWeight(i, k);

                TotalDisc[k] += DiscLength[i][k] * xw / 1000.0;
            }

        }

        void GetMarketSamples(short k) {
            short i, j;
            long ix, ns;
            double xr;


            ZeroMatrix(SampleBins, NBins, NMarket);

            NSample = 0;
            for (j = 0; j < NMarket; j++) {
                ns = (long) MarketSampleSize[j][k];
                for (ix = 0; ix < ns; ix++) {
                    xr = DRNUNF();
                    for (i = 0; i < NBins; i++) {
                        if (xr < MarketDist[i + j * NBins][k]) {
                            ++SampleBins[i][j];
                            break;
                        }
                    }

                }
                NSample += (long) MarketSampleSize[j][k];
            }

            for (i = 0; i < NBins; i++) {
                SampleLength[i] = 0.0;

                for (j = 0; j < NMarket; j++)
                    SampleLength[i] += SampleBins[i][j];
            }
        }

        void OriginalCatchAgeKey(short k) {
            short i, j, n;


            ZeroMatrix(AgeLengthKey, NBins, NAges);

            for (i = 0; i < NBins; i++) {
                n = i + k * NBins;
                for (j = 0; j < NAges; j++)
                    AgeLengthKey[i][j] = CatchAgeKey[n][j];
            }

            if (MisAgeFlag[0])
                MultMatrix(NBins, NAges, NAges, AgeLengthKey, MisAge, AgeKeyCalc);
            else {
                for (i = 0; i < NBins; i++)

                    for (j = 0; j < NAges; j++)
                        AgeKeyCalc[i][j] = AgeLengthKey[i][j];
            }


        }

        void SampleCatchAgeKey(short k) {
            short i, j, n;
            long ns, kk;
            double xt, xr, xa, xn;

            ZeroMatrix(AgeLengthKey, NBins, NAges);

            for (i = 0; i < NBins; i++) {
                xn = SampleLength[i] * SampFrac[0];
                ns = (long) (xn + 0.5);
                if (ns == 0 && SampleLength[i] > 0.0)
                    ns = 1;
                for (kk = 0; kk < ns; kk++) {
                    xr = DRNUNF();
                    xa = 0;
                    for (j = 0; j < NAges; j++) {
                        n = i + k * NBins;
                        xa += CatchAgeKey[n][j];
                        if (xr < xa) {
                            ++AgeLengthKey[i][j];
                            break;
                        }
                    }
                }
            }

            for (i = 0; i < NBins; i++) {
                xt = 0.;
                for (j = 0; j < NAges; j++)
                    xt += AgeLengthKey[i][j];


                for (j = 0; j < NAges; j++) {
                    if (xt > 0.)
                        AgeLengthKey[i][j] = AgeLengthKey[i][j] / xt;
                    else
                        AgeLengthKey[i][j] = 0.;
                }
            }

            if (MisAgeFlag[0])
                MultMatrix(NBins, NAges, NAges, AgeLengthKey, MisAge, AgeKeyCalc);
            else {
                for (i = 0; i < NBins; i++)

                    for (j = 0; j < NAges; j++)
                        AgeKeyCalc[i][j] = AgeLengthKey[i][j];
            }

        }

        void ExpandSamples(short k) {
            short i, j, n;
            double xl, xw, xt, xp;
            double alpha, beta;
            double sd, zx, err;


            alpha = LenWtCoeff[k][4];
            beta = LenWtCoeff[k][5];


            for (j = 0; j < NMarket; j++) {
                xt = 0.0;

                for (i = 0; i < NBins; i++) {
                    xl = (double) (MinLen + i);
                    xw = exp(alpha + beta * log(xl));

                    xt += SampleBins[i][j] * xw;
                }

                n = k * NMarket + j;
                sd = sqrt(log(CVLand * CVLand + 1.0));
                zx = LandDev[n];
                err = exp(zx * sd);
                MarketTarg[j][k] = MarketLand[j][k] * err;

                if (xt > 0.0)
                    xp = MarketTarg[j][k] / xt * 1000.;
                else
                    xp = 0.0;

                for (i = 0; i < NBins; i++)
                    SampleBins[i][j] = SampleBins[i][j] * xp;

            }



            for (i = 0; i < NBins; i++) {
                ExpandLength[i][k] = 0.0;
                for (j = 0; j < NMarket; j++)
                    ExpandLength[i][k] += SampleBins[i][j];
            }



            for (i = 0; i < NBins; i++)
                for (j = 0; j < NAges; j++)
                    SampleAges[i][j] = ExpandLength[i][k] * AgeKeyCalc[i][j];

            for (j = 0; j < KAges - 1; j++) {
                ExpandCatch[j][k] = 0.0;
                for (i = 0; i < NBins; i++)
                    ExpandCatch[j][k] += SampleAges[i][j];
            }

            ExpandCatch[KAges - 1][k] = 0.0;
            for (j = KAges - 1; j < NAges; j++) {

                for (i = 0; i < NBins; i++)
                    ExpandCatch[KAges - 1][k] += SampleAges[i][j];
            }

        }

        void CalcExpandWeights(short k) {
            short i, j;
            double xl, xt, xw, xx;
            double alpha, beta;

            alpha = LenWtCoeff[k][4];
            beta = LenWtCoeff[k][5];


            for (j = 0; j < KAges - 1; j++) {
                xt = 0.0;
                xx = 0.0;
                for (i = 0; i < NBins; i++) {
                    xl = (double) (MinLen + i);
                    xw = exp(alpha + beta * log(xl));
                    xt += SampleAges[i][j];
                    xx += xw * SampleAges[i][j];
                }
                if (xt > 0.0)
                    ExpandWeights[j][k] = xx / xt;
                else
                    ExpandWeights[j][k] = CatchWeights[j][k];
            }

            xx = 0.0;
            xt = 0.0;
            for (j = KAges - 1; j < NAges; j++) {
                for (i = 0; i < NBins; i++) {
                    xl = (double) (MinLen + i);
                    xw = exp(alpha + beta * log(xl));
                    xt += SampleAges[i][j];
                    xx += xw * SampleAges[i][j];
                }
            }
            if (xt > 0.0)
                ExpandWeights[KAges - 1][k] = xx / xt;

            else
                ExpandWeights[KAges - 1][k] = CatchWeights[KAges - 1][k];
        }

        void GetSurveySamples(short k, short ns) {
            short i, n;
            double xt, xn, xr;
            long nx, kx;

            xn = 0.0;
            for (i = 0; i < NBins; i++) {
                n = ns * NBins + i;

                SurveyStock[i] = StockLength[i][k] * SurveySelect[n][k];

                xn += SurveyStock[i];
            }

            xt = ceil(SurveyQ[ns][k] * xn);

            n = ns * (MaxYears + 1) + k;
            KSamples[n] = (long) (xt);

            for (i = 0; i < NBins; i++)
                SurveyStock[i] = SurveyStock[i] * xt / xn;

            xn = 0.0;
            for (i = 0; i < NBins; i++) {
                xn += SurveyStock[i] / xt;
                SurveyDist[i] = xn;

                SurveySamples[i] = 0.0;
            }

            nx = KSamples[n];
            for (kx = 0; kx < nx; kx++) {
                xr = DRNUNF();
                for (i = 0; i < NBins; i++) {
                    xn = SurveyDist[i];
                    if (xr < xn) {
                        SurveySamples[i]++;
                        break;
                    }
                }
            }

            for (i = 0; i < NBins; i++) {

                n = ns * NBins + i;
                SurveyLength[n][k] = SurveySamples[i];
            }


        }

        void ApplySurveyAgeKey(short k, short ns) {
            short i, j, n;
            double xt, xn, xr, xa;
            long nx, kx;


            ZeroMatrix(AgeLengthKey, NBins, NAges);

            if (SampFlag[ns + 1] == 0) {

                for (i = 0; i < NBins; i++) {
                    n = i + k * NBins;
                    for (j = 0; j < NAges; j++)
                        AgeLengthKey[i][j] = StockAgeKey[n][j];
                }

                if (MisAgeFlag[ns + 1])
                    MultMatrix(NBins, NAges, NAges, AgeLengthKey, MisAge, AgeKeyCalc);
                else {
                    for (i = 0; i < NBins; i++)
                        for (j = 0; j < NAges; j++)
                            AgeKeyCalc[i][j] = AgeLengthKey[i][j];
                }
            } else {
                for (i = 0; i < NBins; i++) {
                    n = ns * NBins + i;
                    xn = SurveyLength[n][k] * SampFrac[ns + 1];
                    nx = (long) (xn + 0.5);
                    if (nx == 0 && SurveyLength[n][k] > 0.0)
                        nx = 1;
                    for (kx = 0; kx < nx; kx++) {
                        xr = DRNUNF();
                        xa = 0;
                        for (j = 0; j < NAges; j++) {
                            n = i + k * NBins;
                            xa += StockAgeKey[n][j];
                            if (xr < xa) {
                                ++AgeLengthKey[i][j];
                                break;
                            }
                        }
                    }
                }

                for (i = 0; i < NBins; i++) {
                    xt = 0.;
                    for (j = 0; j < NAges; j++)
                        xt += AgeLengthKey[i][j];


                    for (j = 0; j < NAges; j++) {
                        if (xt > 0.)
                            AgeLengthKey[i][j] = AgeLengthKey[i][j] / xt;
                        else
                            AgeLengthKey[i][j] = 0.;
                    }
                }

                if (MisAgeFlag[ns + 1])
                    MultMatrix(NBins, NAges, NAges, AgeLengthKey, MisAge, AgeKeyCalc);
                else {
                    for (i = 0; i < NBins; i++)
                        for (j = 0; j < NAges; j++)
                            AgeKeyCalc[i][j] = AgeLengthKey[i][j];
                }
            }

            for (i = 0; i < NBins; i++) {
                n = ns * NBins + i;
                for (j = 0; j < NAges; j++) {
                    SurveyBins[i][j] = SurveyLength[n][k] * AgeKeyCalc[i][j];
                }
            }

            n = ns * (MaxYears + 1) + k;
            for (j = 0; j < KAges - 1; j++) {
                xt = 0.0;
                for (i = 0; i < NBins; i++)
                    xt += SurveyBins[i][j];

                SurveyAge[n][j] = xt;
            }

            SurveyAge[n][KAges - 1] = 0.0;
            for (j = KAges - 1; j < NAges; j++) {
                xt = 0.0;
                n = ns * (MaxYears + 1) + k;

                for (i = 0; i < NBins; i++)
                    xt += SurveyBins[i][j];

                SurveyAge[n][KAges - 1] += xt;
            }

        }

        void CreateIndexData(short kx) {
            short j, k, n, nx;
            short j1, j2;

            /* Complete Index Data */

            for (n = 0; n < NIndex; n++) {
                if (kx == IndexData[n].SurveyIndex - 1) {
                    for (k = 0; k < NYears + 1; k++) {
                        nx = kx * (MaxYears + 1) + k;
                        IndexValues[k][n] = 0.0;
                        j1 = IndexData[n].StartAge - NFAge;
                        j2 = IndexData[n].EndAge - NFAge + 1;

                        for (j = j1; j < j2; j++)
                            if (IndexData[n].Type == 0)
                                IndexValues[k][n] += SurveyAge[nx][j];
                            else if (IndexData[n].Type == 1 && j < KAges)
                                IndexValues[k][n] += SurveyAge[nx][j] * StockWeights[j][k];

                            else
                                IndexValues[k][n] += SurveyAge[nx][j] * StockWeights[KAges - 1][k];
                    }
                }
            }

        }

        void AddIndexData(short kx) {
            short j, n, nx;
            short j1, j2;

            /* Complete Index Data */

            for (n = 0; n < NIndex; n++) {
                if (kx == IndexData[n].SurveyIndex - 1) {
                    nx = kx * (MaxYears + 1) + NYears;
                    IndexValues[NYears][n] = 0.0;
                    j1 = IndexData[n].StartAge - NFAge;
                    j2 = IndexData[n].EndAge - NFAge + 1;

                    for (j = j1; j < j2; j++)
                        if (IndexData[n].Type == 0)
                            IndexValues[NYears][n] += SurveyAge[nx][j];
                        else if (IndexData[n].Type == 1 && j < KAges)
                            IndexValues[NYears][n] += SurveyAge[nx][j] * StockWeights[j][NYears];

                        else
                            IndexValues[NYears][n] += SurveyAge[nx][j] * StockWeights[KAges - 1][NYears];
                }
            }

        }

        void ReadVPATemplateFile(char *fn) {
            char buffer[MAXBUF];
            char *tok;
            char *c;
            short i, j;

            if ((c = strrchr(fn, '.')) != NULL)
                *c = '\0';

            strcat(fn, ".tp1");

            if ((fp1 = fopen(fn, "r")) == NULL) {
                fprintf(stderr, "Unable to Open Template File: %s\n", fn);
                exit(1);
            }


            SurveyFacFlag = 0;
            CMultFlag = 0;



            while (!feof(fp1)) {
                fgets(buffer, MAXBUF, fp1);
                if (strstr(buffer, "PATH")) {
                    fgets(buffer, MAXBUF, fp1);
                    if ((c = strchr(buffer, '\n')) != NULL)
                        *c = '\0';
                    strcpy(VPAPath, buffer);
                } else if (strstr(buffer, "FOLD")) {
                    fgets(buffer, MAXBUF, fp1);
                    tok = strtok(buffer, " \t\r\n");
                    Fold[0] = atoi(tok);
                    tok = strtok(NULL, " \t\r\n");
                    Fold[1] = atoi(tok);
                } else if (strstr(buffer, "PRVEC")) {
                    fgets(buffer, MAXBUF, fp1);
                    tok = strtok(buffer, " \t\r\n");
                    for (i = 0; i < KAges - 1; i++) {
                        PRVec[i] = atof(tok);
                        tok = strtok(NULL, " \t\r\n");
                    }
                } else if (strstr(buffer, "STOCKEST")) {
                    fgets(buffer, MAXBUF, fp1);
                    tok = strtok(buffer, " \t\r\n");
                    for (i = 0; i < KAges; i++) {
                        StockEst[i] = atof(tok);
                        tok = strtok(NULL, " \t\r\n");
                    }
                } else if (strstr(buffer, "OPTIONS")) {
                    fgets(buffer, MAXBUF, fp1);
                    tok = strtok(buffer, " \t\r\n");
                    NBoot = atoi(tok);

                    fgets(buffer, MAXBUF, fp1);
                    tok = strtok(buffer, " \t\r\n");
                    for (i = 0; i < NOPT; i++) {
                        OptionFlags[i] = atoi(tok);
                        tok = strtok(NULL, " \t\r\n");
                    }
                } else if (strstr(buffer, "RANGE")) {
                    fgets(buffer, MAXBUF, fp1);
                    tok = strtok(buffer, " \t\r\n");
                    NRange[0] = atoi(tok);
                    tok = strtok(NULL, " \t\r\n");
                    NRange[1] = atoi(tok);
                } else if (strstr(buffer, "NMORT")) {
                    fgets(buffer, MAXBUF, fp1);
                    tok = strtok(buffer, " \t\r\n");
                    NMortEst = atof(tok);
                } else if (strstr(buffer, "CATCH_ADJUST")) {
                    for (i = 0; i < MaxYears; i++) {
                        fgets(buffer, MAXBUF, fp1);
                        tok = strtok(buffer, " \t\r\n");
                        CatchFactors[i] = atof(tok);
                    }
                } else if (strstr(buffer, "SURVEY FACTORS")) {
                    SurveyFacFlag = 1;
                    for (i = 0; i < MaxYears + 1; i++) {
                        fgets(buffer, MAXBUF, fp1);
                        tok = strtok(buffer, " \t\r\n");
                        for (j = 0; j < NIndex; j++) {
                            VPASurveyFactors[i][j] = atof(tok);
                            tok = strtok(NULL, " \t\r\n");
                        }
                    }
                } else if (strstr(buffer, "CMULT")) {

                    CMultFlag = 1;
                    fgets(buffer, MAXBUF, fp1);
                    tok = strtok(buffer, " \t\r\n");
                    CMStart = atoi(tok);
                    tok = strtok(NULL, " \t\r\n");
                    CMEnd = atoi(tok);

                    fgets(buffer, MAXBUF, fp1);
                    tok = strtok(buffer, " \t\r\n");
                    CMxlb = atof(tok);
                    tok = strtok(NULL, " \t\r\n");
                    CMxub = atof(tok);
                }
            }


            fclose(fp1);

        }

        void WriteVPAInputFile(char *fn) {
            char xname[FILBUF];
            char *c;
            short i, j;
            short KAgeEst;
            double fd;
            double fplus = 1.0;
            double stkmin = 1.0;
            double stkmax = 1.0E+08;

            strcpy(xname, fn);

            if ((c = strrchr(xname, '.')) != NULL)
                *c = '\0';

            strcat(xname, "_VPA.dat");

            if ((fp1 = fopen(xname, "w")) == NULL) {
                fprintf(stderr, "Unable to Open VPA Input File: %s\n", xname);
                exit(1);
            }

            for (i = 0; i < KAges; i++)
                AgeEst[i] = -1;

            KAgeEst = 0;
            for (i = 0; i < KAges; i++) {
                if (StockEst[i] > 0.0) {
                    AgeEst[i] = NFAge + i;
                    KAgeEst++;
                }
            }

            fprintf(fp1, "VPA/ADAPT V 3.0\n");
            fprintf(fp1, "##\n");
            fprintf(fp1, "MODEL ID\n");
            fprintf(fp1, "Test Case\n");
            fprintf(fp1, "PARAM\n");
            fprintf(fp1, "%d ", NYears);
            fprintf(fp1, "%d ", KAges);
            fprintf(fp1, "%d ", NIndex);
            fprintf(fp1, "%d ", NIndex);
            fprintf(fp1, "%d ", KAgeEst);
            fprintf(fp1, "%d ", NFYear);
            fprintf(fp1, "%d ", NFAge);
            fprintf(fp1, "%d ", NBoot);
            fprintf(fp1, "%ld\n", BaseSeed);

            fprintf(fp1, "PARTIAL RECRUIT\n");
            for (i = 0; i < KAges - 1; i++)
                fprintf(fp1, "%-6.4f ", PRVec[i]);
            fprintf(fp1, "\n");

            fprintf(fp1, "AGE ESTIMATE\n");
            for (i = 0; i < KAges; i++) {
                if (AgeEst[i] > -1)
                    fprintf(fp1, "%-2d ", AgeEst[i]);
            }
            fprintf(fp1, "\n");

            fprintf(fp1, "STOCK ESTIMATE\n");
            for (i = 0; i < KAges; i++) {
                if (AgeEst[i] > -1)
                    fprintf(fp1, "%10.0f ", StockEst[i]);
            }
            fprintf(fp1, "\n");

            fprintf(fp1, "F-PLUS\n");
            fprintf(fp1, "%-6.4f*%-d\n", fplus, NYears);

            fprintf(fp1, "STOCK MIN-MAX\n");
            fprintf(fp1, "%-8.2E %-8.2E\n", stkmin, stkmax);

            fprintf(fp1, "AGES-F\n");
            fprintf(fp1, "%-d-%-d\n", Fold[0], Fold[1]);

            fprintf(fp1, "AGES-SUMMARY\n");
            fprintf(fp1, "%-d-%-d\n", Fold[0], Fold[1]);

            fprintf(fp1, "MFSPAWN\n");
            fprintf(fp1, "%-6.4f  %-6.4f\n", TM, TF);

            fprintf(fp1, "CATCH AT AGE\n");
            for (i = 0; i < NYears; i++) {
                if (DiscFlag) {
                    for (j = 0; j < KAges; j++) {
                        fd = DiscAge[j][i] / (DiscAge[j][i] + CatchAge[j][i]);
                        fprintf(fp1, "%-18.8E ", ExpandCatch[j][i] * CatchFactors[i] / (1.0 - fd));
                    }
                    fprintf(fp1, "\n");
                } else {
                    for (j = 0; j < KAges; j++)
                        fprintf(fp1, "%-12.1f ", ExpandCatch[j][i] * CatchFactors[i]);
                    fprintf(fp1, "\n");
                }
            }

            for (j = 0; j < KAges; j++)
                DFrac[j] = DiscAge[j][NYears - 1] / (DiscAge[j][NYears - 1] + CatchAge[j][NYears - 1]);

            fprintf(fp1, "WEIGHT AT AGE\n");
            for (i = 0; i < NYears; i++) {
                for (j = 0; j < KAges; j++)
                    fprintf(fp1, "%-12.3f ", ExpandWeights[j][i]);
                fprintf(fp1, "\n");
            }

            fprintf(fp1, "BIOMASS\n");
            for (i = 0; i < NYears + 1; i++) {
                for (j = 0; j < KAges; j++)
                    fprintf(fp1, "%-12.3f ", StockWeights[j][i]);
                fprintf(fp1, "\n");
            }

            fprintf(fp1, "SSB\n");
            for (i = 0; i < NYears; i++) {
                for (j = 0; j < KAges; j++)
                    fprintf(fp1, "%-12.3f ", SpStockWeights[j][i]);
                fprintf(fp1, "\n");
            }

            fprintf(fp1, "MATURITY\n");
            for (i = 0; i < NYears; i++) {
                for (j = 0; j < KAges; j++)
                    fprintf(fp1, "%-12.3f ", Mature[j][i]);
                fprintf(fp1, "\n");
            }

            fprintf(fp1, "M MATRIX\n");
            for (i = 0; i < NYears; i++) {
                for (j = 0; j < KAges; j++)
                    fprintf(fp1, "%-12.3f ", NMortEst);
                fprintf(fp1, "\n");
            }

            if (SurveyFacFlag) {
                fprintf(fp1, "SURVEY FACTORS\n");
                for (i = 0; i < NYears + 1; i++) {
                    for (j = 0; j < NIndex; j++)
                        fprintf(fp1, "%-15.6f ", VPASurveyFactors[i][j]);
                    fprintf(fp1, "\n");
                }
            }

            if (CMultFlag) {
                fprintf(fp1, "CMULT\n");
                fprintf(fp1, "%d  %d\n", CMStart, CMEnd);
                fprintf(fp1, "%10.4f  %10.4f\n", CMxlb, CMxub);
            }


            fprintf(fp1, "SURVEY INDEX\n");
            for (i = 0; i < NIndex; i++)
                fprintf(fp1, "%-15s ", IndexData[i].Name);
            fprintf(fp1, "\n");

            for (i = 0; i < NIndex; i++)
                if (IndexData[i].StartAge == IndexData[i].EndAge)
                    fprintf(fp1, "%-15d ", IndexData[i].StartAge);
                else if (IndexData[i].StartAge < NAgePlus && IndexData[i].EndAge <= NAgePlus)
                    fprintf(fp1, "%7d:%-7d ", IndexData[i].StartAge, IndexData[i].EndAge);
                else
                    fprintf(fp1, "%-15d ", NAgePlus);
            fprintf(fp1, "\n");

            for (i = 0; i < NIndex; i++)
                fprintf(fp1, "1-Jan           ");
            fprintf(fp1, "\n");

            for (i = 0; i < NIndex; i++)
                if (IndexData[i].Type == 0)
                    fprintf(fp1, "Number          ");
                else
                    fprintf(fp1, "Weight          ");
            fprintf(fp1, "\n");

            for (i = 0; i < NYears + 1; i++) {
                for (j = 0; j < NIndex; j++)
                    fprintf(fp1, "%-18.8E ", IndexValues[i][j]);

                fprintf(fp1, "\n");
            }

            fprintf(fp1, "CHECKED INDEX\n");
            for (i = 0; i < NIndex; i++)
                fprintf(fp1, "%-d ", i + 1);
            fprintf(fp1, "\n");

            fprintf(fp1, "CHECKED RETRO\n");
            for (i = 0; i < NIndex; i++)
                fprintf(fp1, "%-d ", 1);
            fprintf(fp1, "\n");

            fprintf(fp1, "OPTIONS\n");
            for (i = 0; i < NOPT; i++)
                fprintf(fp1, "%d  ", OptionFlags[i]);
            fprintf(fp1, "\n");

            if (AgeEst[0] != NFAge) {

                fprintf(fp1, "YEAR_RANGE\n");
                fprintf(fp1, "%d-%d\n", NRange[0], NRange[1]);
            }

            fprintf(fp1, "BOOTFAC\n");
            fprintf(fp1, "1000\n");

            fprintf(fp1, "XPARM\n");
            fprintf(fp1, "80   1\n");

            fclose(fp1);

        }

        void RemoveOutputFiles(char *fn) {
            char xname[MAXBUF];
            char *c;

            strcpy(xname, fn);

            if ((c = strrchr(xname, '.')) != NULL)
                *c = '\0';

            strcat(xname, "_vpa.out");

            remove(xname);

            if ((c = strrchr(xname, '.')) != NULL)
                *c = '\0';

            strcat(xname, ".aux");

            remove(xname);

            *c = '\0';

            strcat(xname, ".pp2");

            remove(xname);

            *c = '\0';

            strcat(xname, ".bsn");

            remove(xname);

            *c = '\0';

            strcat(xname, ".bsf");

            remove(xname);

            strcpy(xname, fn);

            if ((c = strrchr(xname, '.')) != NULL)
                *c = '\0';

            strcat(xname, "_apro.in");

            remove(xname);

            if ((c = strrchr(xname, '.')) != NULL)
                *c = '\0';

            strcat(xname, ".out");

            remove(xname);

            /* Aspic */

            strcpy(xname, fn);

            if ((c = strrchr(xname, '.')) != NULL)
                *c = '\0';

            strcat(xname, "_aspic.bot");

            remove(xname);

            if ((c = strrchr(xname, '.')) != NULL)
                *c = '\0';

            strcat(xname, ".bio");

            remove(xname);
            if ((c = strrchr(xname, '.')) != NULL)
                *c = '\0';

            strcat(xname, ".prb");

            remove(xname);

            if ((c = strrchr(xname, '.')) != NULL)
                *c = '\0';

            strcat(xname, ".prj");

            remove(xname);


        }

        void ScanVPAResults(char *fn) {
            char buffer[MAXBUF];
            char xname[FILBUF];
            char *tok, *c;
            short i, j;

            /* Open VPA Auxilliary Output File */

            strcpy(xname, fn);

            if ((c = strrchr(xname, '.')) != NULL)
                *c = '\0';

            strcat(xname, "_vpa.pp2");

            if ((fp1 = fopen(xname, "r")) == NULL) {
                fprintf(stderr, "Unable to Open VPA Output File: %s\n", xname);
                exit(1);
            }

            for (i = 0; i < KAges; i++) {
                fgets(buffer, MAXBUF - 1, fp1);
                tok = strtok(buffer, " \t\r\n");
                for (j = 0; j < NYears + 1; j++) {
                    StockVPA[i][j] = atof(tok);
                    tok = strtok(NULL, " \t\r\n");
                }
            }

            fclose(fp1);

            if ((c = strrchr(xname, '.')) != NULL)
                *c = '\0';

            strcat(xname, ".aux");

            if ((fp1 = fopen(xname, "r")) == NULL) {
                fprintf(stderr, "Unable to Open VPA Output File: %s\n", xname);
                exit(1);
            }

            fgets(buffer, MAXBUF - 1, fp1);
            for (i = 0; i < KAges; i++) {
                fgets(buffer, MAXBUF - 1, fp1);
                tok = strtok(buffer, " \t\r\n");
                for (j = 0; j < NYears; j++) {
                    FVPA[i][j] = atof(tok);
                    tok = strtok(NULL, " \t\r\n");
                }
            }

            /* Skip JAN-1 Biomass */

            for (i = 0; i < KAges + 2; i++)
                fgets(buffer, MAXBUF - 1, fp1);

            for (i = 0; i < KAges; i++) {
                fgets(buffer, MAXBUF - 1, fp1);
                tok = strtok(buffer, " \t\r\n");
                for (j = 0; j < NYears; j++) {
                    SSBVPA[i][j] = atof(tok);
                    tok = strtok(NULL, " \t\r\n");
                }
            }

            fclose(fp1);

            if (NYears == MaxYears) {
                for (i = 0; i < KAges; i++) {

                    for (j = 0; j < NYears; j++)
                        fprintf(fp5, "%10.4f ", FVPA[i][j]);
                    fprintf(fp5, "\n");
                }
                fprintf(fp5, "\n");
            }


        }

        long LaunchVPA(char *fn) {
            std::cout << __func__ << "\n";
            char xname[FILBUF];
            char cmdline[MAXBUF];
            char *c;
            long k;
            char qt[2];

            qt[0] = '\x22';
            qt[1] = '\0';

            strcpy(xname, fn);

            if ((c = strrchr(xname, '.')) != NULL)
                *c = '\0';

            strcat(xname, "_vpa.dat");

            //            VPAPath  ="c:\nft\vpav301\vpa2.exe";
            strcpy(cmdline, VPAPath);
            strcat(cmdline, " ");
            strcat(cmdline, qt);
            strcat(cmdline, xname);
            strcat(cmdline, qt);
            LowerCase(cmdline);

            k = CreateConsoleProcess(cmdline);

            return k;
        }

        long LaunchAgePro() {
            std::cout << __func__ << "\n";

            char cmdline[FILBUF];
            long k;


            //            strcpy(cmdline, "cmd.exe  /C ");
            strcat(cmdline, AgeProPath);
            strcat(cmdline, " < AgePro34.ctrl");
            LowerCase(cmdline);

            k = CreateConsoleProcess(cmdline);

            return k;
        }

        long LaunchAspic(char *fn) {
            std::cout << __func__ << "\n";

            char xname[FILBUF];
            char cmdline[MAXBUF];
            char *c;
            long k;
            char qt[2];

            qt[0] = '\x22';
            qt[1] = '\0';

            strcpy(xname, fn);

            if ((c = strrchr(xname, '.')) != NULL)
                *c = '\0';

            strcat(xname, "_aspic.inp");


            strcpy(cmdline, AspicPath);
            strcat(cmdline, " ");
            strcat(cmdline, qt);
            strcat(cmdline, xname);
            strcat(cmdline, qt);
            LowerCase(cmdline);

            k = CreateConsoleProcess(cmdline);

            return k;
        }

        void CalcMarketSamples(short k) {
            short j;
            long xs;

            for (j = 0; j < NMarket; j++) {

                xs = (long) (MarketSampleRate[j][k] * MarketLand[j][k] / 100. + 0.5);
                MarketSampleSize[j][k] = (double) xs;
            }


        }

        double GetTargetF(double *b) {
            short i, k;
            double x, f;

            k = NYears - 1;

            x = 0.;
            for (i = 0; i < KAges; i++)
                x += SSBVPA[i][k];

            x = x / 1000.0;

            if (x <= Blim)
                f = 0.0;
            else if (x < Bmsy)
                f = x * Fmsy / Bmsy;
            else
                f = Fmsy;


            *b = x;

            return f;
        }

        double GetAspicTargetF(double *b) {
            double x, f;

            x = AspicB[NYears - 1];

            if (x <= Blim)
                f = 0.0;
            else if (x < Bmsy)
                f = x * Fmsy / Bmsy;
            else
                f = Fmsy;


            *b = x;

            return f;
        }

        void SetupAgepro(char *fn) {
            std::cout << __func__ << "\n";
            char xname[FILBUF];
            char *c;

            if ((fp6 = fopen("AgePro34.ctrl", "w")) == NULL) {
                fprintf(stderr, "Unable to Open AgePro Text File: %s\n", xname);
                exit(1);
            }


            strcpy(xname, fn);

            if ((c = strrchr(xname, '.')) != NULL)
                *c = '\0';

            strcat(xname, "_apro.in");

            fprintf(fp6, "%-s\n", xname);

            if ((c = strrchr(xname, '.')) != NULL)
                *c = '\0';

            strcat(xname, ".out");

            fprintf(fp6, "%-s\n", xname);

            fclose(fp6);

        }

        void ReadAgeproTemplate(char *fn) {
            char xname[FILBUF];
            char buffer[MAXBUF];
            char *c;
            char *tok;

            strcpy(xname, fn);

            if ((c = strrchr(xname, '.')) != NULL)
                *c = '\0';

            strcat(xname, ".tpx");

            if ((fp6 = fopen(xname, "r")) == NULL) {
                fprintf(stderr, "Unable to Open AgePro Template File: %s\n", xname);
                exit(1);
            }

            fgets(AgeProPath, MAXBUF, fp6);

            if ((c = strchr(AgeProPath, '\n')) != NULL)
                *c = '\0';

            fgets(buffer, MAXBUF, fp6);
            tok = strtok(buffer, " \t\r\n");
            NAgeProSim = atoi(tok);

            fgets(buffer, MAXBUF, fp6);
            tok = strtok(buffer, " \t\r\n");
            AgeProCV = atof(tok);

            tok = strtok(NULL, " \t\r\n");
            AgeProBias = atof(tok);


            fclose(fp6);
        }

        void WriteAgeProInputFile(char *fn) {
            char xname[FILBUF];
            char *c;
            short i;
            double x, xx;
            double sd, zx, err, tl;

            strcpy(xname, fn);

            if ((c = strrchr(xname, '.')) != NULL)
                *c = '\0';

            strcat(xname, "_apro.in");

            if ((fp1 = fopen(xname, "w")) == NULL) {
                fprintf(stderr, "Unable to Open AgePro Input File: %s\n", xname);
                exit(1);
            }

            fprintf(fp1, "AGEPRO VERSION 3.3\n");
            fprintf(fp1, "MSE\n");
            fprintf(fp1, "%d\n", NFYear + NYears);
            fprintf(fp1, "2\n");
            fprintf(fp1, "%d\n", NAgeProSim);
            fprintf(fp1, "1\n");
            fprintf(fp1, "%-ld\n", BaseSeed);
            fprintf(fp1, "1\n");
            fprintf(fp1, "1\n");
            fprintf(fp1, "%d\n", DiscFlag);
            fprintf(fp1, "0\n");
            fprintf(fp1, "0\n");
            fprintf(fp1, "0\n");
            fprintf(fp1, "0\n");
            fprintf(fp1, "0\n");
            fprintf(fp1, "0\n");
            fprintf(fp1, "0\n");
            fprintf(fp1, "0\n");
            fprintf(fp1, "1\n");
            fprintf(fp1, "0\n");
            fprintf(fp1, "0\n");
            fprintf(fp1, "0\n");
            fprintf(fp1, "1\n");
            fprintf(fp1, "0\n");
            fprintf(fp1, "%d  1  %d\n", KAges, KAges);
            fprintf(fp1, "%-10.4f\n", NMortEst);
            for (i = 0; i < KAges; i++) {
                x = (NMSelect[i][0] + NMSelect[i][1]) / 2.0;
                fprintf(fp1, "%-10.4f ", x);
            }
            fprintf(fp1, "\n");

            for (i = 0; i < KAges; i++)
                fprintf(fp1, "%-10.3f ", StockWeights[i][NYears - 1]);
            fprintf(fp1, "\n");

            for (i = 0; i < KAges; i++)
                fprintf(fp1, "%-10.3f ", CatchWeights[i][NYears - 1]);
            fprintf(fp1, "\n");

            for (i = 0; i < KAges; i++)
                fprintf(fp1, "%-10.3f ", SpStockWeights[i][NYears - 1]);
            fprintf(fp1, "\n");

            for (i = 0; i < KAges; i++)
                fprintf(fp1, "%-10.3f ", CatchWeights[i][NYears - 1]);
            fprintf(fp1, "\n");

            if (DiscFlag) {
                for (i = 0; i < KAges; i++)
                    fprintf(fp1, "%-10.3f ", CatchWeights[i][NYears - 1]);
                fprintf(fp1, "\n");
            }

            for (i = 0; i < KAges; i++)
                fprintf(fp1, "%-10.3f ", Mature[i][NYears - 1]);
            fprintf(fp1, "\n");


            /* Fraction of Z before Spawning */


            xx = (NMortEst * TM + FFull[NYears - 1] * TF) / (NMortEst + FFull[NYears - 1]);

            fprintf(fp1, "%10.4f\n", xx);

            /* Recruitment */

            fprintf(fp1, "14\n");
            fprintf(fp1, "%d\n", NYears);
            for (i = 0; i < NYears; i++)
                fprintf(fp1, "%18.3f ", StockVPA[0][i] / 1000.);
            fprintf(fp1, "\n");

            fprintf(fp1, "1.0\n");
            fprintf(fp1, "1.0\n");

            fprintf(fp1, "%-ld\n", NBoot);

            strcpy(xname, fn);

            if ((c = strrchr(xname, '.')) != NULL)
                *c = '\0';

            strcat(xname, "_vpa.bsn");

            fprintf(fp1, "%-s\n", xname);
            fprintf(fp1, "1000\n");

            PRVec[KAges - 1] = 1.000;
            for (i = 0; i < KAges; i++)
                fprintf(fp1, "%-10.4f ", PRVec[i]);
            fprintf(fp1, "\n");

            if (DiscFlag) {

                for (i = 0; i < KAges; i++)
                    fprintf(fp1, "%-10.4f ", DFrac[i]);
                fprintf(fp1, "\n");
            }

            sd = sqrt(log(AgeProCV * AgeProCV + 1.0));
            zx = AgeProDev[NYears - 1];
            err = exp(zx * sd);

            tl = TargetLand * (1.0 + AgeProBias) * err;

            fprintf(fp1, "1  0\n");
            fprintf(fp1, "%10.0f  -1\n", tl * 1000.);
            fprintf(fp1, "-1   %10.5f\n", TargetF);

            fclose(fp1);

        }

        void ScanAgeProReport(char *fn) {
            char xname[FILBUF];
            char buffer[MAXBUF];
            char *c;
            char *tok;
            short i;
            short flag;

            strcpy(xname, fn);

            c = strrchr(xname, '.');
            *c = '\0';

            strcat(xname, "_apro.out");

            if ((fp1 = fopen(xname, "r")) == NULL) {
                fprintf(stderr, "Unable to Open Agepro Output File: %s\n", xname);
                exit(1);
            }

            flag = 0;
            while (!feof(fp1)) {
                fgets(buffer, MAXBUF, fp1);
                if (strstr(buffer, "PERCENTILES OF LANDINGS")) {
                    flag = 1;
                    break;
                }
            }



            if (!flag) {
                fprintf(stderr, "Missing Data From Agepro Report\n");
                exit(1);
            }

            for (i = 0; i < 3; i++)
                fgets(buffer, MAXBUF, fp1);


            tok = strtok(buffer, " \t\r\n");
            for (i = 0; i < 5; i++)
                tok = strtok(NULL, " \t\r\n");

            ProjLand = atof(tok) * 1000.0;

            if (DiscFlag) {
                flag = 0;
                while (!feof(fp1)) {
                    fgets(buffer, MAXBUF, fp1);
                    if (strstr(buffer, "PERCENTILES OF DISCARDS")) {
                        flag = 1;
                        break;
                    }
                }



                if (!flag) {
                    fprintf(stderr, "Missing Data From Agepro Report\n");
                    exit(1);
                }

                for (i = 0; i < 3; i++)
                    fgets(buffer, MAXBUF, fp1);


                tok = strtok(buffer, " \t\r\n");

                for (i = 0; i < 5; i++)
                    tok = strtok(NULL, " \t\r\n");

                ProjLand += atof(tok) * 1000.0;
            }

            fclose(fp1);
        }

        double CalculateProjectedF(double tl) {
            short iter;
            short i, j, k;
            double m, f, z, ff;
            double xl;
            double ad, bd;
            double xc, xd, fd;
            double alpha, beta;
            double xw, xxl, xxx;
            double xx, err1, err2, slope;

            tl = tl * 1000.;

            k = NYears - 1;

            ad = DiscFrac[k][0];
            bd = DiscFrac[k][1];


            /* Initial Guess */

            ff = 0.5;

            for (iter = 0; iter < MXITER; iter++) {

                /* Males */

                alpha = LenWtCoeff[k][0];
                beta = LenWtCoeff[k][1];

                /* Apply Mortality to Males */

                xxl = 0.0;
                for (i = 0; i < NBins; i++) {
                    xl = (double) (MinLen + i);
                    xw = exp(alpha + beta * log(xl));

                    fd = DiscardFraction(xl, ad, bd);
                    xd = 1.0 - fd;

                    f = ff * FSelect[i][k];


                    for (j = 0; j < NAges; j++) {

                        m = NMort[k][0] * NMSelect[j][0];

                        z = f + m;

                        xc = f * (1.0 - exp(-z)) / z;

                        xxx = MaleBins[i][j] * xc;

                        xxl += xxx * xw;

                    }
                }

                /* Females */

                alpha = LenWtCoeff[k][2];
                beta = LenWtCoeff[k][3];

                /* Apply Mortality to Females */

                for (i = 0; i < NBins; i++) {
                    xl = (double) (MinLen + i);
                    xw = exp(alpha + beta * log(xl));

                    fd = DiscardFraction(xl, ad, bd);
                    xd = 1.0 - fd;

                    f = ff * FSelect[i][k];


                    for (j = 0; j < NAges; j++) {
                        m = NMort[k][1] * NMSelect[j][1];

                        z = f + m;

                        xc = f * (1.0 - exp(-z)) / z;

                        xxx = FemaleBins[i][j] * xc;

                        xxl += xxx * xw;

                    }
                }



                err1 = (xxl - tl) / tl;
                if (fabs(err1) < XTOL)
                    break;

                if (iter == 0) {
                    xx = ff;
                    err2 = err1;
                    if (xxl > tl)
                        ff = ff * 0.9;
                    else
                        ff = ff * 1.1;
                } else {
                    slope = (err2 - err1) / (xx - ff);
                    xx = ff;
                    err2 = err1;
                    ff = (slope * xx - err2) / slope;
                }

                if (ff > XUBOUND)
                    ff = XUBOUND;

                if (ff < XLBOUND)
                    ff = XLBOUND;
            }

            return ff;

        }

        void ReadAspicTemplateFile(char *fn) {
            char buffer[MAXBUF];
            char *c;
            short i;

            if ((c = strrchr(fn, '.')) != NULL)
                *c = '\0';

            strcat(fn, ".tp3");

            if ((fp1 = fopen(fn, "r")) == NULL) {
                fprintf(stderr, "Unable to Open Template File: %s\n", fn);
                exit(1);
            }


            AspicF = (double *) calloc(MaxYears, sizeof (double));
            AspicB = (double *) calloc(MaxYears, sizeof (double));
            AspicY = (double *) calloc(MaxYears, sizeof (double));


            fgets(buffer, MAXBUF, fp1);

            if ((c = strchr(buffer, '\n')) != NULL)
                *c = '\0';

            strcpy(AspicPath, buffer);
            strcpy(AspicpPath, AspicPath);

            if ((c = strrchr(AspicpPath, '.')) != NULL)
                *c = '\0';

            strcat(AspicpPath, "p.exe");

            for (i = 0; i < NSTRING; i++) {

                fgets(buffer, MAXBUF, fp1);
                PString[i] = _strdup(buffer);
            }



            fclose(fp1);

        }

        void WriteAspicInputFile(char *fn) {
            char xname[FILBUF];
            char *c;
            char quote = '\x22';
            short i, j, k;
            double xt, fd;

            strcpy(xname, fn);

            if ((c = strrchr(xname, '.')) != NULL)
                *c = '\0';

            strcat(xname, "_aspic.inp");

            if ((fp1 = fopen(xname, "w")) == NULL) {
                fprintf(stderr, "Unable to Open ASPIC Input File: %s\n", xname);
                exit(1);
            }

            for (i = 0; i < NSTRING; i++)
                fputs(PString[i], fp1);

            fprintf(fp1, "%d                     ## Number of Years In Catch Data\n", NYears);

            fprintf(fp1, "%cCPUE  & Yield%c\n", quote, quote);

            fprintf(fp1, "CC\n");

            for (i = 0; i < NYears; i++) {
                xt = 0.0;
                for (j = 0; j < KAges; j++) {
                    if (DiscFlag) {
                        fd = DiscAge[j][i] / (DiscAge[j][i] + CatchAge[j][i]);
                        xt += ExpandCatch[j][i] * ExpandWeights[j][i] / 1000. / (1.0 - fd);
                    } else
                        xt += ExpandCatch[j][i] * ExpandWeights[j][i] / 1000.;
                }

                fprintf(fp1, "%d    %18.8E  %18.8E\n", NFYear + i, IndexValues[i][0] / 1000., xt);
            }
            for (k = 1; k < NSurvey; k++) {
                fprintf(fp1, "%cSurvey #%d Biomass Index%c\n", quote, k + 1, quote);
                fprintf(fp1, "I0\n");

                for (i = 0; i < NYears; i++)
                    fprintf(fp1, "%d   %18.8E\n", NFYear + i, IndexValues[i][k] / 1000.);
            }

            fclose(fp1);
        }

        void CaptureAspicResults(char *fn) {
            char buffer[MAXBUF];
            char xname[FILBUF];
            char *tok, *c;
            short i;

            /* Open ASPIC Output File */

            strcpy(xname, fn);

            if ((c = strrchr(xname, '.')) != NULL)
                *c = '\0';

            strcat(xname, "_aspic.bot");

            if ((fp1 = fopen(xname, "r")) == NULL) {
                fprintf(stderr, "Unable to Open ASPIC Output File: %s\n", xname);
                exit(1);
            }

            while (!feof(fp1)) {

                fgets(buffer, MAXBUF - 1, fp1);

                if (strstr(buffer, "ESTIMATED POPULATION TRAJECTORY")) {
                    for (i = 0; i < 6; i++)
                        fgets(buffer, MAXBUF - 1, fp1);

                    for (i = 0; i < NYears; i++) {

                        fgets(buffer, MAXBUF - 1, fp1);
                        tok = strtok(buffer, " \t\r\n");
                        tok = strtok(NULL, " \t\r\n");
                        tok = strtok(NULL, " \t\r\n");
                        AspicF[i] = atof(tok);
                        tok = strtok(NULL, " \t\r\n");
                        AspicB[i] = atof(tok);
                        tok = strtok(NULL, " \t\r\n");
                        tok = strtok(NULL, " \t\r\n");
                        tok = strtok(NULL, " \t\r\n");
                        AspicY[i] = atof(tok);
                    }
                }

            }


            fclose(fp1);
        }

        void WriteAspicPFile(char *fn) {
            char xname[FILBUF];
            char *c;
            char quote = '\x22';
            double sd, zx, err, tl;

            strcpy(xname, fn);

            if ((c = strrchr(xname, '.')) != NULL)
                *c = '\0';

            strcat(xname, "_aspic.ctl");

            if ((fp1 = fopen(xname, "w")) == NULL) {
                fprintf(stderr, "Unable to Open ASPIC Projection Input File: %s\n", xname);
                exit(1);
            }

            if ((c = strrchr(xname, '.')) != NULL)
                *c = '\0';

            strcat(xname, ".bio");

            fprintf(fp1, "%cMSE Version 3.2%c\n", quote, quote);

            fprintf(fp1, "%-s\n", xname);

            fprintf(fp1, "XX\n");

            fprintf(fp1, "PC  1\n");

            fprintf(fp1, "0\n");

            fprintf(fp1, "2\n");

            sd = sqrt(log(AgeProCV * AgeProCV + 1.0));
            zx = AgeProDev[NYears - 1];
            err = exp(zx * sd);

            tl = TargetLand * (1.0 + AgeProBias) * err;

            fprintf(fp1, "%10.3f   Y\n", tl);

            fprintf(fp1, "%10.4f   F\n", TargetF / FFull[NYears - 1]);

            fclose(fp1);
        }

        long LaunchAspicProj(char *fn) {
            std::cout << __func__ << "\n";

            char xname[FILBUF];
            char cmdline[MAXBUF];
            char *c;
            long k;
            char qt[2];

            qt[0] = '\x22';
            qt[1] = '\0';

            strcpy(xname, fn);

            if ((c = strrchr(xname, '.')) != NULL)
                *c = '\0';

            strcat(xname, "_aspic.ctl");


            strcpy(cmdline, AspicpPath);
            strcat(cmdline, " ");
            strcat(cmdline, qt);
            strcat(cmdline, xname);
            strcat(cmdline, qt);
            LowerCase(cmdline);

            k = CreateConsoleProcess(cmdline);

            return k;
        }

        void ScanAspicProj(char *fn) {
            char xname[FILBUF];
            char buffer[MAXBUF];
            char *c;
            char *tok;
            short i, flag;

            strcpy(xname, fn);

            if ((c = strrchr(xname, '.')) != NULL)
                *c = '\0';

            strcat(xname, "_aspic.prj");

            if ((fp1 = fopen(xname, "r")) == NULL) {
                fprintf(stderr, "Unable to Open ASPIC Projection Input File: %s\n", xname);
                exit(1);
            }

            flag = 0;
            while (!feof(fp1)) {
                fgets(buffer, MAXBUF, fp1);
                if (strstr(buffer, "TABLE OF PROJECTED YIELDS")) {
                    flag = 1;
                    break;
                }
            }

            if (!flag) {
                fprintf(stderr, "Missing Data From Aspic Projection Report\n");
                exit(1);
            }

            for (i = 0; i < 3; i++)
                fgets(buffer, MAXBUF, fp1);

            tok = strtok(buffer, " \t\r\n");
            tok = strtok(NULL, " \t\r\n");

            ProjLand = atof(tok);

            fclose(fp1);


        }

        void CopyMisAgeMatrix(short kx) {
            short i, j, k, n, nx, ix;

            ix = NFYear + kx;

            nx = 0;

            for (k = 0; k < NAgeErr; k++) {
                if (ix >= AgeErrStart[k] && ix <= AgeErrEnd[k]) {
                    nx = k + 1;
                    break;
                }
            }

            if (nx) {
                for (i = 0; i < NAges; i++) {
                    n = (nx - 1) * NAges + i;
                    for (j = 0; j < NAges; j++) {
                        MisAge[i][j] = AgeErr[n][j];
                    }
                }
            }
        }

    };

}

#endif /* MSE_HPP */

