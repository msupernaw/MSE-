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

    short ParseMultiToken(char *, double *);

    void ZeroMatrix(double **, short, short);

    void MultMatrix(short, short, short, double **, double **, double **);

    void BubbleUpSort(long, double *);

    double MeanValue(long, double *);

    double StdDev(long, double, double *);

    double StdDevP(long, double, double *);


    double **AllocMatrix(short, short);


    /* Date Time Strings */

    char dstring[12];
    char tstring[6];



#ifdef __cplusplus
}
#endif

#endif /* UTIL_H */

