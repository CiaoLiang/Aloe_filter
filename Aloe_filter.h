/**************************************************************************

* FILE


*


*    ALOE_Filter.h


*


* DESCRIPTION


*


*    This file define the functions that used to complemment CCI cancellation.


*


* NOTE


*


*    None


*


******************************************************************************/


 

 

 

 

/******************************************************************************


*  HISTORY OF CHANGES


*******************************************************************************


*


*    <date>        <author>       <CR_ID>        <DESCRIPTION>


*    


*


******************************************************************************/


#ifdef __cplusplus


extern "C"{


#endif


 

 

#ifndef ALOE_FILTER_H


#define ALOE_FILTER_H


 

 

/******************************************************************************


*  INCLUDING FILES


******************************************************************************/


#include <stdlib.h>


#include <stdio.h>


#include "lib_typedef_f.h"


#include "lib_arith_f.h"

 

 

/******************************************************************************


*  FUNCTION PROTOTYPE DECLARATION


******************************************************************************/

void ALOE_Filter(int     flag,

                 cmplx_t y1[156],

                            cmplx_t y2[156],

                            uint8   TSC_LEN,

                            real_t  t[156],

                            uint8   FILTER_LEN,

                            uint8   Delay,

                            cmplx_t bk[156],

                            real_t  data_out[156],

                 real_t  debug_port[] );   //test by gu 2015-08-19


 

void CalculateZ(cmplx_t data[156],

                           uint8   TSC_LEN,

                           uint8   FILTER_LEN,

                           uint8   DELAY,

                           real_t  z[156][20]);

 

void CalculateRzz(real_t z[156][20],

                  uint8  TSC_LEN,

                             uint8  FILTER_LEN,

                             uint8  DELAY,

                  real_t rzz[20][20],

                  real_t scale[1]); //for fix


 

void CalculateP(real_t z[156][20],

                           real_t t[156],

                           uint8  TSC_LEN,

                           uint8  FILTER_LEN,

                           uint8  DELAY,

                           real_t p[20]);

 

void CalculateW(real_t invrzz[20][20],

                           real_t p[20],

                           uint8  FILTER_LEN,

                           real_t w[20],

                real_t scale);


 

void EstimateBk(cmplx_t data[156],

                real_t  w[20],

                           uint8   FILTER_LEN,

                           uint8   DELAY,

                cmplx_t bk[156]);

 

void Inv_complex_matrix(real_t Fac[][20],

                                        real_t FacInvMat[][20],

                        uint8  MatDimension);


 

 

#endif


#ifdef __cplusplus


}


#endif


/*$Log$*/