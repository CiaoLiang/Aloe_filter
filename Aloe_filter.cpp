/***********************************************************************************************************

* FileName::    <ALOE_Filter.c>

* version::     <1.0.0>

* Purpose:      <This file includes the processing functions for linear equalization filter.>

* Authors:      <Liang Wenqiao, 2015-08-12>

* Notes:                   First  is the situation that the two way of data processed in different way than added.

*               second is the situation that the extended TS method.

***********************************************************************************************************/

/***********************************************************************************************************

*  HISTORY OF CHANGES

*   <Date>          <Author>        <Version>       <DESCRIPTION>

*  2015-08-12  Liang Wenqiao           V1.0.0          original

***********************************************************************************************************/

/*-----------including external files -----------------------------*/

#include "math.h"

#include "ALOE_Filter.h"

/*-----------external variable declaration--------------------------*/

 

/*-----------file-local constant and type definition----------------*/

 

/*-----------file-local macro definition----------------------------*/

#define new_inv

/*-----------file-local variables definition------------------------*/

 

/*----------- function definition-----------------------------------*/

/******************************************************************************

*

* FUNCTION

*

*   CalculateZ()

*

* DESCRIPTION

*

*   This function is calculating the matrix Z of received data.

*

* PARAMETERS

*

*       data[]                        in,              input data

*       TSC_LEN                  in,              length of training sequence    

*       z[][]                            out,           received signal matrix

*

* VERSION

*

*    <date>        <author>       <Version>        <DESCRIPTION>

*  2015-08-12       Liang Wenqiao         V1.0.0            original

*

**********************************************************************************/

void CalculateZ(cmplx_t data[156],

                           uint8   TSC_LEN,

                           uint8   FILTER_LEN,

                           uint8   DELAY,

                           real_t  z[156][20])

{

         int i,j;

         if(TSC_LEN==26)

         {

                   for(i=0;i<TSC_LEN;i++)

                   {

                            for(j=0;j<FILTER_LEN;j++)

                            {

                                     z[i][j] = data[61+DELAY+i-j].re;

                                     z[i][j+FILTER_LEN] = 0-data[61+DELAY+i-j].im;

                            }

                   }

         }

         else if(TSC_LEN==156)

         {

                   for(i=DELAY;i<TSC_LEN-DELAY;i++)

                   {

                            for(j=0;j<FILTER_LEN;j++)

                            {

                                     z[i][j] = data[i-j+DELAY].re;

                                     z[i][j+FILTER_LEN] = 0-data[i-j+DELAY].im;

                            }

                   }

         }

}

 

/******************************************************************************

*

* FUNCTION

*

*   CalculateP()

*

* DESCRIPTION

*

*   This function is calculating the crosscorrelation between the burst sequence

*   and the received data.

*bk[]                            out,           soft bits

* PARAMETERS

*

*       z[][]                            in,              received signal matrix

*       t[]                     in,              training sequence

*       TSC_LEN                  in,              length of training sequence

*       p[]                     out,           crosscorrelation vector

*

* VERSION

*

*    <date>        <author>       <Version>        <DESCRIPTION>

*  2015-08-12       Liang Wenqiao         V1.0.0            original

*

**********************************************************************************/

void CalculateP(real_t z[156][20],

                           real_t t[156],

                           uint8  TSC_LEN,

                           uint8  FILTER_LEN,

                           uint8  DELAY,

                           real_t p[20])

{

         int i,j;

    int StartIdx;

         real_t temp = 0.0;

    if(TSC_LEN == 26)

    {

        StartIdx = 0;

    }

    else

    {

        StartIdx = DELAY;

    }

 

         for(i=0;i<2*FILTER_LEN;i++)

         {

                   for(j=StartIdx;j<TSC_LEN-StartIdx;j++)

                   {

                            temp += z[j][i]*t[j];

                   }

                   p[i] = temp;

                   temp = 0;

         }

}

/******************************************************************************

*

* FUNCTION

*

*   CalculateRzz()

*

* DESCRIPTION

*

*   This function is calculating the autocorrelation of received data.

*

* PARAMETERS

*

*       z[][]                            in,              received signal matrix

*       TSC_LEN                  in,              length of training sequence

*       rzz[][]               out,           autocorrelation matrix

*

* VERSION

*

*    <date>        <author>       <Version>        <DESCRIPTION>

*  2015-08-12       Liang Wenqiao         V1.0.0            original

*

**********************************************************************************/

void CalculateRzz(real_t z[156][20],

                  uint8  TSC_LEN,

                             uint8  FILTER_LEN,

                             uint8  DELAY,

                  real_t rzz[20][20],

                  real_t scale[1])

{

         int i,j,k;

    int StartIdx;

         real_t temp = 0.0;

    if(TSC_LEN == 26)

    {

        StartIdx = 0;

    }

    else

    {

        StartIdx = DELAY;

    }

 

         for(i=0; i<2*FILTER_LEN;i++)

         {

                   for(j=0;j<2*FILTER_LEN;j++)

                   {

                            for(k=StartIdx;k<TSC_LEN-StartIdx;k++)

                            {

                                     temp += z[k][i]*z[k][j];

                            }

                            rzz[i][j] = temp;                          

                            temp = 0;

                   }

         }

    temp = 0;

         for(i=0; i<2*FILTER_LEN;i++)

         {

        if(rzz[i][i] > temp)

        {

            temp = rzz[i][i];

        }

    }

    scale[0] = temp;

         for(i=0; i<2*FILTER_LEN;i++)

         {

                   for(j=0;j<2*FILTER_LEN;j++)

                   {

#if defined (new_inv)

                            rzz[i][j] = rzz[i][j]/temp;

#else

            rzz[i][j] = rzz[i][j];

#endif

                   }

         }

}

 

/******************************************************************************

*

* FUNCTION

*

*   Inv_complex_matrix()

*

* DESCRIPTION

*

*   This function is calculating the filter coefficient.

*

* PARAMETERS

*

*       Fac[][]              in,              input matrix

*       FacInvMat[][]          out,  inv matrix

*       MatDimension       in,              matrix dimension

*

* VERSION

*

*    <date>        <author>       <Version>        <DESCRIPTION>

*  2015-08-12       Liang Wenqiao     bk[]                    from LTE

*

**********************************************************************************/

void Inv_complex_matrix(

            real_t                 Fac[][20],

                            real_t                 FacInvMat[][20],

                            uint8                  MatDimension

                            )

{

         uint8         RowIndex,ColIndex;                                                                         //row index,column index

         uint8         Loop_Index;                                                                                                  //loop index

         real_t       FacInvTmp[30][30];                                                                            //temp matrix for calculating FacInv

         real_t       TempVarM_1,TempVarM_2,TempVarM_3,TempVarM_4;                //temp variant of type real_t

 

         for (Loop_Index=0; Loop_Index<MatDimension; Loop_Index++)

         {       

                   //calculate 1/Fac(0,0)

                   //first iteration

                   if (0 == Loop_Index)

                   {                          

                            //Q(N,3)->Q(M,3)

                            //TempVarM_1 = shl_fr1xM((real_t)Fac[0][0], bitWidth());

            TempVarM_1 = Fac[0][0];

                            //FacInvTmp[MatDimension - 1][MatDimension - 1] = divs_i1xM((real_t)(1 << (2*bitWidth() - 2)), TempVarM_1, 10); //Q(M.9) 

                            FacInvTmp[MatDimension - 1][MatDimension - 1] = 1.0 / TempVarM_1;

                   }

                   //otherwise

                   else

                   {

                           

                            //get FacInvMat(0,0)

                            TempVarM_1 = FacInvMat[0][0];

                            //Modified by liuqingwei 2009-06-29

                            ///FacInvTmp[MatDimension - 1][MatDimension - 1] = divs_i1xM((real_t)(1 << (2*bitWidth() - 2)), TempVarM_1, 16); //Q(M.9)

                            FacInvTmp[MatDimension - 1][MatDimension - 1] = 1.0 / TempVarM_1;

           

 

                   }

 

                   //calculate elements of upper triangle of the inverse matrix

                   for (RowIndex=0; RowIndex<MatDimension; RowIndex++)

                   {

                            for (ColIndex=RowIndex; ColIndex<MatDimension; ColIndex++)

                            {

                                     //first iteration

                                     if (0 == Loop_Index)

                                     {

                                               if ((RowIndex != MatDimension - 1) && (ColIndex == MatDimension - 1))

                                               {

                                                        //get Fac(i+1,MatDimension-1)

                                                        //TempVarN_1 = Fac[RowIndex + 1][0];

                        TempVarM_1 = Fac[RowIndex + 1][0];

                        //FacInvTmp[RowIndex][MatDimension - 1] = mlt_frMXfrM_frM(TempVarM_1,FacInvTmp[MatDimension - 1][MatDimension - 1]);

                                                        FacInvTmp[RowIndex][MatDimension - 1] = TempVarM_1 * FacInvTmp[MatDimension - 1][MatDimension - 1];

                                               }

                                               else if ((RowIndex != MatDimension - 1) && (RowIndex == ColIndex))

                                               {

                                                       

                                                        //get Fac(i+1,0)

                                                        //TempVarN_1 = Fac[RowIndex + 1][0];

                        TempVarM_3 = Fac[RowIndex + 1][0];

                                                        //get Fac(0,j+1)

                                                        //TempVarN_2 = Fac[0][ColIndex + 1];

                        TempVarM_4 = Fac[0][ColIndex + 1];

                                                        //get Fac(i+1,j+1)

                        //TempVarM_1 = shr_fr1xM(Fac[RowIndex + 1][ColIndex + 1], 6); //Q(M,3)->Q(M,9)

                                                        TempVarM_1 = Fac[RowIndex + 1][ColIndex + 1];

                        //TempVarM_2 = mlt_frMXfrM_frM(TempVarM_3,TempVarM_4);

                                                        TempVarM_2 = TempVarM_3 * TempVarM_4 * FacInvTmp[MatDimension - 1][MatDimension - 1];

                                                        //FacInvTmp[RowIndex][ColIndex] = sub_fr1xM(TempVarM_1,TempVarM_2);

                                                        FacInvTmp[RowIndex][ColIndex] = TempVarM_1 - TempVarM_2;                                               

 

                                               }

                                               else if ((RowIndex != MatDimension - 1) || (ColIndex != MatDimension - 1))

                                               {

                                                       

                                                        //get Fac(i+1,0)

                                                        //TempVarN_1 = Fac[RowIndex + 1][0];

                        TempVarM_3 = Fac[RowIndex + 1][0];

                                                        //get Fac(0,j+1)

                        TempVarM_4 = Fac[0][ColIndex + 1];

                                                        //get Fac(i+1,j+1)

                        //TempVarM_2 = mlt_frMXfrM_frM(TempVarM_2,FacInvTmp[MatDimension - 1][MatDimension - 1]);

                                                        //FacInvTmp[RowIndex][ColIndex] = sub_fr1xM(TempVarM_1,TempVarM_2);    

                                                        TempVarM_1 = Fac[RowIndex + 1][ColIndex + 1];

                                                        TempVarM_2 = TempVarM_3 * TempVarM_4 * FacInvTmp[MatDimension - 1][MatDimension - 1];

                                                        FacInvTmp[RowIndex][ColIndex] = TempVarM_1 - TempVarM_2;

                                                       

                                               }

                                     }

                                     //otherwise

                                     else

                                     {

                                               if ((RowIndex != MatDimension - 1) && (ColIndex == MatDimension - 1))

                                               {

                                                        //get Fac(i+1,MatDimension-1)                       

                                                        TempVarM_1 = FacInvMat[RowIndex + 1][0];

                        TempVarM_2 = FacInvTmp[MatDimension - 1][MatDimension - 1];// add                       

                                                       

                        FacInvTmp[RowIndex][ColIndex] = TempVarM_1 * TempVarM_2;

                                               }

                                               else if ((RowIndex != MatDimension - 1) || (ColIndex != MatDimension - 1))

                                               {

                                                       

                                                        //get Fac(i+1,0)

                                                        TempVarM_1 = FacInvMat[RowIndex + 1][0];

                                                        //get Fac(0,j+1)

                                                        TempVarM_2 = FacInvMat[0][ColIndex + 1];

                                                        //get Fac(i+1,j+1)

                                                        TempVarM_3 = FacInvMat[RowIndex + 1][ColIndex + 1];                                                     

                                                        TempVarM_4 = TempVarM_1 * TempVarM_2 * FacInvTmp[MatDimension - 1][MatDimension - 1];

                                                        FacInvTmp[RowIndex][ColIndex] = TempVarM_3 - TempVarM_4;                                               

                                                       

                                               }

                                     }

                            }

                   }

                   //calculate elements of lower triangle of the inverse matrix

                   for (RowIndex=0; RowIndex<MatDimension; RowIndex++)

                   {

                            for (ColIndex=0; ColIndex<RowIndex; ColIndex++)

                            {

                                     if ((RowIndex <= MatDimension - Loop_Index - 2) || (ColIndex >= MatDimension - Loop_Index - 1))

                                     {

                                               FacInvTmp[RowIndex][ColIndex] = FacInvTmp[ColIndex][RowIndex];

                                     }

                                     else

                                     {

                                               //FacInvTmp[RowIndex][ColIndex] = negate_fr1xM(FacInvTmp[ColIndex][RowIndex]);

                                               FacInvTmp[RowIndex][ColIndex] = - FacInvTmp[ColIndex][RowIndex];

                                     }

                            }

                   }

 

                   //update FacInvMat

                   for (RowIndex=0; RowIndex<MatDimension; RowIndex++)

                   {

                            for (ColIndex=0; ColIndex<MatDimension; ColIndex++)

                            {

                                     FacInvMat[RowIndex][ColIndex] = FacInvTmp[RowIndex][ColIndex];

                            }

                   }

 

         } //end for (Loop_Index=0; Loop_Index<MatDimension; Loop_Index++)

}

 

/******************************************************************************

*

* FUNCTION

*

*   EstimateBk()

*

* DESCRIPTION

*

*   This function is calculating the filtered received data.

*

* PARAMETERS

*

*       data[]                        in,              input data

*       w[]                    in,              filter coefficent        

*       bk[]                            out,           soft bits

*

* VERSION

*

*    <date>        <author>       <Version>        <DESCRIPTION>

*  2015-08-12       Liang Wenqiao         V1.0.0            original

*      

*

***************************************************************************/

void EstimateBk(cmplx_t data[156],

                real_t w[20],

                           uint8  FILTER_LEN,

                           uint8  DELAY,

                cmplx_t bk[156])

{

    cmplx_t w_temp[20];

         cmplx_t data_temp[200]={0};

         int i,j;

         for(i=0;i<FILTER_LEN;i++)

         {

                   w_temp[i].re = w[i];

        w_temp[i].im = w[i + FILTER_LEN];

         }

    algo_gsm_eq_gmsk_eic_flt_complex_fltp(data, 156, w_temp,  5,  1,  0, data_temp);

       for(i=0;i<156;i++)

         {

                   bk[i].re = data_temp[i+DELAY-2].re;

        bk[i].im = data_temp[i+DELAY-2].im;

         }  

}

 

/******************************************************************************

*

* FUNCTION

*

*   CalculateW()

*

* DESCRIPTION

*

*   This function is calculating the filter coefficient.

*

* PARAMETERS

*

*       invrzz[][]                   in,              inv of received signal matrix

*       p[]                     in,              crosscorrelation vector

*       w[]                    out,           filter coefficient

*

* VERSION

*

*    <date>        <author>       <Version>        <DESCRIPTION>

*  2015-08-12       Liang Wenqiao         V1.0.0            original

*

***************************************************************************/

void CalculateW(real_t invrzz[20][20],

                           real_t p[20],

                           uint8  FILTER_LEN,

                           real_t w[20],

                real_t scale)

{

         int i,j;

         real_t temp = 0.0;

         for(i=0;i<2*FILTER_LEN;i++)

         {

                   for(j=0;j<2*FILTER_LEN;j++)

                   {

                            temp += invrzz[i][j]*p[j];

                   }

                   w[i] = temp/scale;

                   temp = 0;

         }

}

 

/******************************************************************************

*

* FUNCTION

*

*   ALOE_Filter()

*

* DESCRIPTION

*

*   This function is calculating soft bits and hard decision of filtered received data.  

*

* PARAMETERS

*

*       y1[]                            in,              received signal odd

*       y2[]                            in,              received signal even

*       t[]                     in,              training sequence

*       bk[]                            out,           soft bits

*       data_out[][]            out,           hard decision

*

* VERSION

*

*    <date>        <author>       <Version>        <DESCRIPTION>

*  2015-08-12       Liang Wenqiao         V1.0.0            original

*

***************************************************************************/

void ALOE_Filter(int     flag,

                 cmplx_t y1[156],

                            cmplx_t y2[156],

                            uint8   TSC_LEN,

                            real_t  t[156],

                            uint8   Filter_LEN,

                            uint8   Delay,

                            cmplx_t bk[156],

                            real_t  data_out[156],

                 real_t  debug_port[]) //test

{                

         int i,j;       

         real_t imag[20][20]={0};

        real_t tsc[156]={0};

         real_t z1[156][20]={0};

         real_t rzz1[20][20]={0};            

         real_t invrzz1[20][20]={0};      

         real_t p1[20]={0};      

         real_t w1[20]={0};

         cmplx_t b1[156]={0};

         real_t z2[156][20]={0};

         real_t rzz2[20][20]={0};          

         real_t invrzz2[20][20]={0};      

         real_t p2[20]={0};         

         real_t w2[20]={0};

         cmplx_t b2[156]={0};

         real_t scale[2]={0};

         real_t sum1,sum2,index1,index2,temp1,temp2;

         //printf("%d\n",TSC_LEN);

         for(i=0; i<TSC_LEN; i++)

       {

                   tsc[i] =1-2*t[i];

                   //printf("tsc[%d]=%f\n",i,tsc[i]);

         }

         //odd                 

         CalculateZ(y1,TSC_LEN,Filter_LEN,Delay,z1);

         for(i=0; i<TSC_LEN; i++)

       {

                   for(j=0; j<2*Filter_LEN; j++)

                {

                            //printf("z1[%d][%d]=%f\n",i,j,z1[i][j]);

                   }

         }

         CalculateRzz(z1,TSC_LEN,Filter_LEN,Delay,rzz1,&(scale[0]));

         for(i=0; i<2*Filter_LEN; i++)

       {

                   for(j=0; j<2*Filter_LEN; j++)

                {

                            printf("rzz1[%d][%d]=%f\n",i,j,rzz1[i][j]);

                   }

         }

      

         CalculateP(z1,tsc,TSC_LEN,Filter_LEN,Delay,p1);

         Inv_complex_matrix(rzz1,invrzz1,2*Filter_LEN);

         CalculateW(invrzz1,p1,Filter_LEN,w1,scale[0]);

         for(i=0; i<2*Filter_LEN; i++)

       {

                   for(j=0; j<2*Filter_LEN; j++)

                {

                          printf("invrzz1[%d][%d]=%f\n",i,j,invrzz1[i][j]);

                   }

         }

         EstimateBk(y1,w1,Filter_LEN,Delay,b1);

         //even

         CalculateZ(y2,TSC_LEN,Filter_LEN,Delay,z2);

         CalculateRzz(z2,TSC_LEN,Filter_LEN,Delay,rzz2,&(scale[1]));

         CalculateP(z2,tsc,TSC_LEN,Filter_LEN,Delay,p2);

         Inv_complex_matrix(rzz2,invrzz2,2*Filter_LEN);

         CalculateW(invrzz2,p2,Filter_LEN,w2,scale[1]);

         EstimateBk(y2,w2,Filter_LEN,Delay,b2);

         for(i=0; i<2*Filter_LEN; i++)

       {

                   printf("w1[%d]=%f\n",i,w1[i]);

                   printf("w2[%d]=%f\n",i,w2[i]);

         }

     //testport Rzz1 by Liang Wenqiao

         for(i=0; i<2*Filter_LEN; i++)

       {

                   for(j=0; j<2*Filter_LEN; j++)

                {

                     debug_port[j + i*2*Filter_LEN] = rzz1[i][j];

                   }

         }

         for(i=0; i<2*Filter_LEN; i++)

       {

                   for(j=0; j<2*Filter_LEN; j++)

                {

                     debug_port[2*Filter_LEN * 2*Filter_LEN + j + i*2*Filter_LEN] = invrzz1[i][j];

                   }

         }

         for(i=0; i<2*Filter_LEN; i++)

       {

         debug_port[2*Filter_LEN * 2*Filter_LEN * 2 + i] = p1[i];

         }

         for(i=0; i<2*Filter_LEN; i++)

       {

         debug_port[2*Filter_LEN * 2*Filter_LEN * 2 + 2*Filter_LEN + i] = w1[i];

         }

 

         sum1 = 0;

         sum2 = 0;

         for (i=62; i<88; i++)

         {

                   temp1 = y1[i].re-tsc[i];

                   temp2 = y2[i].re-tsc[i];

                   sum1 = sum1 + temp1*temp1;

                   sum2 = sum2 + temp2*temp2; ;

    }

         index1 = sum2/(sum1+sum2);

         index2 = sum1/(sum1+sum2);

 

         //hard decision

    if(flag == 0) //even

    {

             for(i=0;i<156;i++)

             {

                      bk[i].re = b1[i].re;

                      bk[i].im = b1[i].im;

                      //bk[i].im = 0;                                      

             }       

             for(i=0; i<61; i++)

           {                                         

                      if(bk[i].re >= 0.0)

                                data_out[i] = 0;                 

                      else                  

                                data_out[i] = 1;       

             }                                   

             for(i=61; i<87; i++)

           {

                      data_out[i] = t[i-61];                 

             }

             for(i=87; i<156; i++)

           {                                         

                      if(bk[i].re >= 0.0)

                                data_out[i] = 0;                 

                      else                  

                                data_out[i] = 1;       

             }

    }

    else if(flag == 1) //odd

    {

             for(i=0;i<156;i++)

             {

                      bk[i].re = b2[i].re;

                      bk[i].im = b2[i].im;

                      //bk[i].im = 0;                                      

             }       

             for(i=0; i<61; i++)

           {                                         

                      if(bk[i].re >= 0.0)

                                data_out[i] = 0;                 

                      else                  

                                data_out[i] = 1;       

             }                                   

             for(i=61; i<87; i++)

           {

                      data_out[i] = t[i-61];                 

             }

             for(i=87; i<156; i++)

           {                                         

                      if(bk[i].re >= 0.0)

                                data_out[i] = 0;                 

                      else                  

                                data_out[i] = 1;       

             }

    }

    else //2x

    {

             for(i=0;i<156;i++)

             {

                      bk[i].re = b1[i].re*index1+b2[i].re*index2;

                      bk[i].im = b1[i].im*index1+b2[i].im*index2;

                      //bk[i].im = 0;                                      

             }       

             for(i=0; i<61; i++)

           {                                         

                      if(bk[i].re >= 0.0)

                                data_out[i] = 0;                 

                      else                  

                                data_out[i] = 1;       

             }                                   

             for(i=61; i<87; i++)

           {

                      data_out[i] = t[i-61];                 

             }

             for(i=87; i<156; i++)

           {                                         

                      if(bk[i].re >= 0.0)

                                data_out[i] = 0;                 

                      else                  

                                data_out[i] = 1;       

             }

    }

}

/*$Log$*/