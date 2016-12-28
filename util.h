/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/*
 * File:   util.h
 * Author: matthewsupernaw
 *
 * Created on January 19, 2016, 11:21 AM
 */

#ifndef UTIL_H
#define UTIL_H

#ifdef __cplusplus
extern "C" {
#endif

    /* Function Prototypes */

    char *UpperCase(char *);
    char *LowerCase(char *);
    char *TrimString(char *);

    void DateString(char *);
    void TimeString(char *);

    short ParseMultiToken(char *,double *);

    void ZeroVector(double *,long);

    void ZeroMatrix(double **, long, long);

    void ZeroMat3(double ***,long,long,long);

    void MultMatrix(long, long, long, double **, double **, double **);

    void BubbleUpSort(long, double *);

    double MeanValue(long, double *);

    double StdDev(long, double, double *);

    double StdDevP(long, double, double *);


    double **AllocMatrix(long, long);

    double ***AllocMat3(long,long,long);

    int compdouble(const void *,const void *);


    /* Date Time Strings */

    char dstring[12];
    char tstring[6];


#ifdef __cplusplus
}
#endif

#endif /* UTIL_H */

