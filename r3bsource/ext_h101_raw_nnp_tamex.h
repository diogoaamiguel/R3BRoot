/********************************************************
 *
 * Structure for ext_data_fetch_event() filling.
 *
 * Do not edit - automatically generated.
 */

#ifndef __GUARD_H101_RAW_NNP_TAMEX_EXT_H101_RAW_NNP_TAMEX_H__
#define __GUARD_H101_RAW_NNP_TAMEX_EXT_H101_RAW_NNP_TAMEX_H__

#ifndef __CINT__
#include <stdint.h>
#else
/* For CINT (old version trouble with stdint.h): */
#ifndef uint32_t
typedef unsigned int uint32_t;
typedef int int32_t;
#endif
#endif
#ifndef EXT_STRUCT_CTRL
#define EXT_STRUCT_CTRL(x)
#endif

/********************************************************
 *
 * Plain structure (layout as ntuple/root file):
 */

typedef struct EXT_STR_h101_raw_nnp_tamex_t
{
    /* RAW */
    uint32_t NN_P1tcl_T1BM /* [1,50] */;
    uint32_t NN_P1tcl_T1BMI[50 EXT_STRUCT_CTRL(NN_P1tcl_T1BM)] /* [1,50] */;
    uint32_t NN_P1tcl_T1BME[50 EXT_STRUCT_CTRL(NN_P1tcl_T1BM)] /* [1,10000] */;
    uint32_t NN_P1tcl_T1B /* [0,10000] */;
    uint32_t NN_P1tcl_T1Bv[10000 EXT_STRUCT_CTRL(NN_P1tcl_T1B)] /* [0,65535] */;
    uint32_t NN_P1tcl_T2BM /* [1,50] */;
    uint32_t NN_P1tcl_T2BMI[50 EXT_STRUCT_CTRL(NN_P1tcl_T2BM)] /* [1,50] */;
    uint32_t NN_P1tcl_T2BME[50 EXT_STRUCT_CTRL(NN_P1tcl_T2BM)] /* [1,10000] */;
    uint32_t NN_P1tcl_T2B /* [0,10000] */;
    uint32_t NN_P1tcl_T2Bv[10000 EXT_STRUCT_CTRL(NN_P1tcl_T2B)] /* [0,65535] */;
    uint32_t NN_P1tfl_T1BM /* [1,50] */;
    uint32_t NN_P1tfl_T1BMI[50 EXT_STRUCT_CTRL(NN_P1tfl_T1BM)] /* [1,50] */;
    uint32_t NN_P1tfl_T1BME[50 EXT_STRUCT_CTRL(NN_P1tfl_T1BM)] /* [1,10000] */;
    uint32_t NN_P1tfl_T1B /* [0,10000] */;
    uint32_t NN_P1tfl_T1Bv[10000 EXT_STRUCT_CTRL(NN_P1tfl_T1B)] /* [0,65535] */;
    uint32_t NN_P1tfl_T2BM /* [1,50] */;
    uint32_t NN_P1tfl_T2BMI[50 EXT_STRUCT_CTRL(NN_P1tfl_T2BM)] /* [1,50] */;
    uint32_t NN_P1tfl_T2BME[50 EXT_STRUCT_CTRL(NN_P1tfl_T2BM)] /* [1,10000] */;
    uint32_t NN_P1tfl_T2B /* [0,10000] */;
    uint32_t NN_P1tfl_T2Bv[10000 EXT_STRUCT_CTRL(NN_P1tfl_T2B)] /* [0,65535] */;
    uint32_t NN_P1tct_T1BM /* [1,50] */;
    uint32_t NN_P1tct_T1BMI[50 EXT_STRUCT_CTRL(NN_P1tct_T1BM)] /* [1,50] */;
    uint32_t NN_P1tct_T1BME[50 EXT_STRUCT_CTRL(NN_P1tct_T1BM)] /* [1,10000] */;
    uint32_t NN_P1tct_T1B /* [0,10000] */;
    uint32_t NN_P1tct_T1Bv[10000 EXT_STRUCT_CTRL(NN_P1tct_T1B)] /* [0,65535] */;
    uint32_t NN_P1tct_T2BM /* [1,50] */;
    uint32_t NN_P1tct_T2BMI[50 EXT_STRUCT_CTRL(NN_P1tct_T2BM)] /* [1,50] */;
    uint32_t NN_P1tct_T2BME[50 EXT_STRUCT_CTRL(NN_P1tct_T2BM)] /* [1,10000] */;
    uint32_t NN_P1tct_T2B /* [0,10000] */;
    uint32_t NN_P1tct_T2Bv[10000 EXT_STRUCT_CTRL(NN_P1tct_T2B)] /* [0,65535] */;
    uint32_t NN_P1tft_T1BM /* [1,50] */;
    uint32_t NN_P1tft_T1BMI[50 EXT_STRUCT_CTRL(NN_P1tft_T1BM)] /* [1,50] */;
    uint32_t NN_P1tft_T1BME[50 EXT_STRUCT_CTRL(NN_P1tft_T1BM)] /* [1,10000] */;
    uint32_t NN_P1tft_T1B /* [0,10000] */;
    uint32_t NN_P1tft_T1Bv[10000 EXT_STRUCT_CTRL(NN_P1tft_T1B)] /* [0,65535] */;
    uint32_t NN_P1tft_T2BM /* [1,50] */;
    uint32_t NN_P1tft_T2BMI[50 EXT_STRUCT_CTRL(NN_P1tft_T2BM)] /* [1,50] */;
    uint32_t NN_P1tft_T2BME[50 EXT_STRUCT_CTRL(NN_P1tft_T2BM)] /* [1,10000] */;
    uint32_t NN_P1tft_T2B /* [0,10000] */;
    uint32_t NN_P1tft_T2Bv[10000 EXT_STRUCT_CTRL(NN_P1tft_T2B)] /* [0,65535] */;
    uint32_t NN_P2tcl_T1BM /* [1,50] */;
    uint32_t NN_P2tcl_T1BMI[50 EXT_STRUCT_CTRL(NN_P2tcl_T1BM)] /* [1,50] */;
    uint32_t NN_P2tcl_T1BME[50 EXT_STRUCT_CTRL(NN_P2tcl_T1BM)] /* [1,10000] */;
    uint32_t NN_P2tcl_T1B /* [0,10000] */;
    uint32_t NN_P2tcl_T1Bv[10000 EXT_STRUCT_CTRL(NN_P2tcl_T1B)] /* [0,65535] */;
    uint32_t NN_P2tcl_T2BM /* [1,50] */;
    uint32_t NN_P2tcl_T2BMI[50 EXT_STRUCT_CTRL(NN_P2tcl_T2BM)] /* [1,50] */;
    uint32_t NN_P2tcl_T2BME[50 EXT_STRUCT_CTRL(NN_P2tcl_T2BM)] /* [1,10000] */;
    uint32_t NN_P2tcl_T2B /* [0,10000] */;
    uint32_t NN_P2tcl_T2Bv[10000 EXT_STRUCT_CTRL(NN_P2tcl_T2B)] /* [0,65535] */;
    uint32_t NN_P2tfl_T1BM /* [1,50] */;
    uint32_t NN_P2tfl_T1BMI[50 EXT_STRUCT_CTRL(NN_P2tfl_T1BM)] /* [1,50] */;
    uint32_t NN_P2tfl_T1BME[50 EXT_STRUCT_CTRL(NN_P2tfl_T1BM)] /* [1,10000] */;
    uint32_t NN_P2tfl_T1B /* [0,10000] */;
    uint32_t NN_P2tfl_T1Bv[10000 EXT_STRUCT_CTRL(NN_P2tfl_T1B)] /* [0,65535] */;
    uint32_t NN_P2tfl_T2BM /* [1,50] */;
    uint32_t NN_P2tfl_T2BMI[50 EXT_STRUCT_CTRL(NN_P2tfl_T2BM)] /* [1,50] */;
    uint32_t NN_P2tfl_T2BME[50 EXT_STRUCT_CTRL(NN_P2tfl_T2BM)] /* [1,10000] */;
    uint32_t NN_P2tfl_T2B /* [0,10000] */;
    uint32_t NN_P2tfl_T2Bv[10000 EXT_STRUCT_CTRL(NN_P2tfl_T2B)] /* [0,65535] */;
    uint32_t NN_P2tct_T1BM /* [1,50] */;
    uint32_t NN_P2tct_T1BMI[50 EXT_STRUCT_CTRL(NN_P2tct_T1BM)] /* [1,50] */;
    uint32_t NN_P2tct_T1BME[50 EXT_STRUCT_CTRL(NN_P2tct_T1BM)] /* [1,10000] */;
    uint32_t NN_P2tct_T1B /* [0,10000] */;
    uint32_t NN_P2tct_T1Bv[10000 EXT_STRUCT_CTRL(NN_P2tct_T1B)] /* [0,65535] */;
    uint32_t NN_P2tct_T2BM /* [1,50] */;
    uint32_t NN_P2tct_T2BMI[50 EXT_STRUCT_CTRL(NN_P2tct_T2BM)] /* [1,50] */;
    uint32_t NN_P2tct_T2BME[50 EXT_STRUCT_CTRL(NN_P2tct_T2BM)] /* [1,10000] */;
    uint32_t NN_P2tct_T2B /* [0,10000] */;
    uint32_t NN_P2tct_T2Bv[10000 EXT_STRUCT_CTRL(NN_P2tct_T2B)] /* [0,65535] */;
    uint32_t NN_P2tft_T1BM /* [1,50] */;
    uint32_t NN_P2tft_T1BMI[50 EXT_STRUCT_CTRL(NN_P2tft_T1BM)] /* [1,50] */;
    uint32_t NN_P2tft_T1BME[50 EXT_STRUCT_CTRL(NN_P2tft_T1BM)] /* [1,10000] */;
    uint32_t NN_P2tft_T1B /* [0,10000] */;
    uint32_t NN_P2tft_T1Bv[10000 EXT_STRUCT_CTRL(NN_P2tft_T1B)] /* [0,65535] */;
    uint32_t NN_P2tft_T2BM /* [1,50] */;
    uint32_t NN_P2tft_T2BMI[50 EXT_STRUCT_CTRL(NN_P2tft_T2BM)] /* [1,50] */;
    uint32_t NN_P2tft_T2BME[50 EXT_STRUCT_CTRL(NN_P2tft_T2BM)] /* [1,10000] */;
    uint32_t NN_P2tft_T2B /* [0,10000] */;
    uint32_t NN_P2tft_T2Bv[10000 EXT_STRUCT_CTRL(NN_P2tft_T2B)] /* [0,65535] */;
    uint32_t NN_P3tcl_T1BM /* [1,50] */;
    uint32_t NN_P3tcl_T1BMI[50 EXT_STRUCT_CTRL(NN_P3tcl_T1BM)] /* [1,50] */;
    uint32_t NN_P3tcl_T1BME[50 EXT_STRUCT_CTRL(NN_P3tcl_T1BM)] /* [1,10000] */;
    uint32_t NN_P3tcl_T1B /* [0,10000] */;
    uint32_t NN_P3tcl_T1Bv[10000 EXT_STRUCT_CTRL(NN_P3tcl_T1B)] /* [0,65535] */;
    uint32_t NN_P3tcl_T2BM /* [1,50] */;
    uint32_t NN_P3tcl_T2BMI[50 EXT_STRUCT_CTRL(NN_P3tcl_T2BM)] /* [1,50] */;
    uint32_t NN_P3tcl_T2BME[50 EXT_STRUCT_CTRL(NN_P3tcl_T2BM)] /* [1,10000] */;
    uint32_t NN_P3tcl_T2B /* [0,10000] */;
    uint32_t NN_P3tcl_T2Bv[10000 EXT_STRUCT_CTRL(NN_P3tcl_T2B)] /* [0,65535] */;
    uint32_t NN_P3tfl_T1BM /* [1,50] */;
    uint32_t NN_P3tfl_T1BMI[50 EXT_STRUCT_CTRL(NN_P3tfl_T1BM)] /* [1,50] */;
    uint32_t NN_P3tfl_T1BME[50 EXT_STRUCT_CTRL(NN_P3tfl_T1BM)] /* [1,10000] */;
    uint32_t NN_P3tfl_T1B /* [0,10000] */;
    uint32_t NN_P3tfl_T1Bv[10000 EXT_STRUCT_CTRL(NN_P3tfl_T1B)] /* [0,65535] */;
    uint32_t NN_P3tfl_T2BM /* [1,50] */;
    uint32_t NN_P3tfl_T2BMI[50 EXT_STRUCT_CTRL(NN_P3tfl_T2BM)] /* [1,50] */;
    uint32_t NN_P3tfl_T2BME[50 EXT_STRUCT_CTRL(NN_P3tfl_T2BM)] /* [1,10000] */;
    uint32_t NN_P3tfl_T2B /* [0,10000] */;
    uint32_t NN_P3tfl_T2Bv[10000 EXT_STRUCT_CTRL(NN_P3tfl_T2B)] /* [0,65535] */;
    uint32_t NN_P3tct_T1BM /* [1,50] */;
    uint32_t NN_P3tct_T1BMI[50 EXT_STRUCT_CTRL(NN_P3tct_T1BM)] /* [1,50] */;
    uint32_t NN_P3tct_T1BME[50 EXT_STRUCT_CTRL(NN_P3tct_T1BM)] /* [1,10000] */;
    uint32_t NN_P3tct_T1B /* [0,10000] */;
    uint32_t NN_P3tct_T1Bv[10000 EXT_STRUCT_CTRL(NN_P3tct_T1B)] /* [0,65535] */;
    uint32_t NN_P3tct_T2BM /* [1,50] */;
    uint32_t NN_P3tct_T2BMI[50 EXT_STRUCT_CTRL(NN_P3tct_T2BM)] /* [1,50] */;
    uint32_t NN_P3tct_T2BME[50 EXT_STRUCT_CTRL(NN_P3tct_T2BM)] /* [1,10000] */;
    uint32_t NN_P3tct_T2B /* [0,10000] */;
    uint32_t NN_P3tct_T2Bv[10000 EXT_STRUCT_CTRL(NN_P3tct_T2B)] /* [0,65535] */;
    uint32_t NN_P3tft_T1BM /* [1,50] */;
    uint32_t NN_P3tft_T1BMI[50 EXT_STRUCT_CTRL(NN_P3tft_T1BM)] /* [1,50] */;
    uint32_t NN_P3tft_T1BME[50 EXT_STRUCT_CTRL(NN_P3tft_T1BM)] /* [1,10000] */;
    uint32_t NN_P3tft_T1B /* [0,10000] */;
    uint32_t NN_P3tft_T1Bv[10000 EXT_STRUCT_CTRL(NN_P3tft_T1B)] /* [0,65535] */;
    uint32_t NN_P3tft_T2BM /* [1,50] */;
    uint32_t NN_P3tft_T2BMI[50 EXT_STRUCT_CTRL(NN_P3tft_T2BM)] /* [1,50] */;
    uint32_t NN_P3tft_T2BME[50 EXT_STRUCT_CTRL(NN_P3tft_T2BM)] /* [1,10000] */;
    uint32_t NN_P3tft_T2B /* [0,10000] */;
    uint32_t NN_P3tft_T2Bv[10000 EXT_STRUCT_CTRL(NN_P3tft_T2B)] /* [0,65535] */;
    uint32_t NN_P4tcl_T1BM /* [1,50] */;
    uint32_t NN_P4tcl_T1BMI[50 EXT_STRUCT_CTRL(NN_P4tcl_T1BM)] /* [1,50] */;
    uint32_t NN_P4tcl_T1BME[50 EXT_STRUCT_CTRL(NN_P4tcl_T1BM)] /* [1,10000] */;
    uint32_t NN_P4tcl_T1B /* [0,10000] */;
    uint32_t NN_P4tcl_T1Bv[10000 EXT_STRUCT_CTRL(NN_P4tcl_T1B)] /* [0,65535] */;
    uint32_t NN_P4tcl_T2BM /* [1,50] */;
    uint32_t NN_P4tcl_T2BMI[50 EXT_STRUCT_CTRL(NN_P4tcl_T2BM)] /* [1,50] */;
    uint32_t NN_P4tcl_T2BME[50 EXT_STRUCT_CTRL(NN_P4tcl_T2BM)] /* [1,10000] */;
    uint32_t NN_P4tcl_T2B /* [0,10000] */;
    uint32_t NN_P4tcl_T2Bv[10000 EXT_STRUCT_CTRL(NN_P4tcl_T2B)] /* [0,65535] */;
    uint32_t NN_P4tfl_T1BM /* [1,50] */;
    uint32_t NN_P4tfl_T1BMI[50 EXT_STRUCT_CTRL(NN_P4tfl_T1BM)] /* [1,50] */;
    uint32_t NN_P4tfl_T1BME[50 EXT_STRUCT_CTRL(NN_P4tfl_T1BM)] /* [1,10000] */;
    uint32_t NN_P4tfl_T1B /* [0,10000] */;
    uint32_t NN_P4tfl_T1Bv[10000 EXT_STRUCT_CTRL(NN_P4tfl_T1B)] /* [0,65535] */;
    uint32_t NN_P4tfl_T2BM /* [1,50] */;
    uint32_t NN_P4tfl_T2BMI[50 EXT_STRUCT_CTRL(NN_P4tfl_T2BM)] /* [1,50] */;
    uint32_t NN_P4tfl_T2BME[50 EXT_STRUCT_CTRL(NN_P4tfl_T2BM)] /* [1,10000] */;
    uint32_t NN_P4tfl_T2B /* [0,10000] */;
    uint32_t NN_P4tfl_T2Bv[10000 EXT_STRUCT_CTRL(NN_P4tfl_T2B)] /* [0,65535] */;
    uint32_t NN_P4tct_T1BM /* [1,50] */;
    uint32_t NN_P4tct_T1BMI[50 EXT_STRUCT_CTRL(NN_P4tct_T1BM)] /* [1,50] */;
    uint32_t NN_P4tct_T1BME[50 EXT_STRUCT_CTRL(NN_P4tct_T1BM)] /* [1,10000] */;
    uint32_t NN_P4tct_T1B /* [0,10000] */;
    uint32_t NN_P4tct_T1Bv[10000 EXT_STRUCT_CTRL(NN_P4tct_T1B)] /* [0,65535] */;
    uint32_t NN_P4tct_T2BM /* [1,50] */;
    uint32_t NN_P4tct_T2BMI[50 EXT_STRUCT_CTRL(NN_P4tct_T2BM)] /* [1,50] */;
    uint32_t NN_P4tct_T2BME[50 EXT_STRUCT_CTRL(NN_P4tct_T2BM)] /* [1,10000] */;
    uint32_t NN_P4tct_T2B /* [0,10000] */;
    uint32_t NN_P4tct_T2Bv[10000 EXT_STRUCT_CTRL(NN_P4tct_T2B)] /* [0,65535] */;
    uint32_t NN_P4tft_T1BM /* [1,50] */;
    uint32_t NN_P4tft_T1BMI[50 EXT_STRUCT_CTRL(NN_P4tft_T1BM)] /* [1,50] */;
    uint32_t NN_P4tft_T1BME[50 EXT_STRUCT_CTRL(NN_P4tft_T1BM)] /* [1,10000] */;
    uint32_t NN_P4tft_T1B /* [0,10000] */;
    uint32_t NN_P4tft_T1Bv[10000 EXT_STRUCT_CTRL(NN_P4tft_T1B)] /* [0,65535] */;
    uint32_t NN_P4tft_T2BM /* [1,50] */;
    uint32_t NN_P4tft_T2BMI[50 EXT_STRUCT_CTRL(NN_P4tft_T2BM)] /* [1,50] */;
    uint32_t NN_P4tft_T2BME[50 EXT_STRUCT_CTRL(NN_P4tft_T2BM)] /* [1,10000] */;
    uint32_t NN_P4tft_T2B /* [0,10000] */;
    uint32_t NN_P4tft_T2Bv[10000 EXT_STRUCT_CTRL(NN_P4tft_T2B)] /* [0,65535] */;
    uint32_t NN_P5tcl_T1BM /* [1,50] */;
    uint32_t NN_P5tcl_T1BMI[50 EXT_STRUCT_CTRL(NN_P5tcl_T1BM)] /* [1,50] */;
    uint32_t NN_P5tcl_T1BME[50 EXT_STRUCT_CTRL(NN_P5tcl_T1BM)] /* [1,10000] */;
    uint32_t NN_P5tcl_T1B /* [0,10000] */;
    uint32_t NN_P5tcl_T1Bv[10000 EXT_STRUCT_CTRL(NN_P5tcl_T1B)] /* [0,65535] */;
    uint32_t NN_P5tcl_T2BM /* [1,50] */;
    uint32_t NN_P5tcl_T2BMI[50 EXT_STRUCT_CTRL(NN_P5tcl_T2BM)] /* [1,50] */;
    uint32_t NN_P5tcl_T2BME[50 EXT_STRUCT_CTRL(NN_P5tcl_T2BM)] /* [1,10000] */;
    uint32_t NN_P5tcl_T2B /* [0,10000] */;
    uint32_t NN_P5tcl_T2Bv[10000 EXT_STRUCT_CTRL(NN_P5tcl_T2B)] /* [0,65535] */;
    uint32_t NN_P5tfl_T1BM /* [1,50] */;
    uint32_t NN_P5tfl_T1BMI[50 EXT_STRUCT_CTRL(NN_P5tfl_T1BM)] /* [1,50] */;
    uint32_t NN_P5tfl_T1BME[50 EXT_STRUCT_CTRL(NN_P5tfl_T1BM)] /* [1,10000] */;
    uint32_t NN_P5tfl_T1B /* [0,10000] */;
    uint32_t NN_P5tfl_T1Bv[10000 EXT_STRUCT_CTRL(NN_P5tfl_T1B)] /* [0,65535] */;
    uint32_t NN_P5tfl_T2BM /* [1,50] */;
    uint32_t NN_P5tfl_T2BMI[50 EXT_STRUCT_CTRL(NN_P5tfl_T2BM)] /* [1,50] */;
    uint32_t NN_P5tfl_T2BME[50 EXT_STRUCT_CTRL(NN_P5tfl_T2BM)] /* [1,10000] */;
    uint32_t NN_P5tfl_T2B /* [0,10000] */;
    uint32_t NN_P5tfl_T2Bv[10000 EXT_STRUCT_CTRL(NN_P5tfl_T2B)] /* [0,65535] */;
    uint32_t NN_P5tct_T1BM /* [1,50] */;
    uint32_t NN_P5tct_T1BMI[50 EXT_STRUCT_CTRL(NN_P5tct_T1BM)] /* [1,50] */;
    uint32_t NN_P5tct_T1BME[50 EXT_STRUCT_CTRL(NN_P5tct_T1BM)] /* [1,10000] */;
    uint32_t NN_P5tct_T1B /* [0,10000] */;
    uint32_t NN_P5tct_T1Bv[10000 EXT_STRUCT_CTRL(NN_P5tct_T1B)] /* [0,65535] */;
    uint32_t NN_P5tct_T2BM /* [1,50] */;
    uint32_t NN_P5tct_T2BMI[50 EXT_STRUCT_CTRL(NN_P5tct_T2BM)] /* [1,50] */;
    uint32_t NN_P5tct_T2BME[50 EXT_STRUCT_CTRL(NN_P5tct_T2BM)] /* [1,10000] */;
    uint32_t NN_P5tct_T2B /* [0,10000] */;
    uint32_t NN_P5tct_T2Bv[10000 EXT_STRUCT_CTRL(NN_P5tct_T2B)] /* [0,65535] */;
    uint32_t NN_P5tft_T1BM /* [1,50] */;
    uint32_t NN_P5tft_T1BMI[50 EXT_STRUCT_CTRL(NN_P5tft_T1BM)] /* [1,50] */;
    uint32_t NN_P5tft_T1BME[50 EXT_STRUCT_CTRL(NN_P5tft_T1BM)] /* [1,10000] */;
    uint32_t NN_P5tft_T1B /* [0,10000] */;
    uint32_t NN_P5tft_T1Bv[10000 EXT_STRUCT_CTRL(NN_P5tft_T1B)] /* [0,65535] */;
    uint32_t NN_P5tft_T2BM /* [1,50] */;
    uint32_t NN_P5tft_T2BMI[50 EXT_STRUCT_CTRL(NN_P5tft_T2BM)] /* [1,50] */;
    uint32_t NN_P5tft_T2BME[50 EXT_STRUCT_CTRL(NN_P5tft_T2BM)] /* [1,10000] */;
    uint32_t NN_P5tft_T2B /* [0,10000] */;
    uint32_t NN_P5tft_T2Bv[10000 EXT_STRUCT_CTRL(NN_P5tft_T2B)] /* [0,65535] */;
    uint32_t NN_P6tcl_T1BM /* [1,50] */;
    uint32_t NN_P6tcl_T1BMI[50 EXT_STRUCT_CTRL(NN_P6tcl_T1BM)] /* [1,50] */;
    uint32_t NN_P6tcl_T1BME[50 EXT_STRUCT_CTRL(NN_P6tcl_T1BM)] /* [1,10000] */;
    uint32_t NN_P6tcl_T1B /* [0,10000] */;
    uint32_t NN_P6tcl_T1Bv[10000 EXT_STRUCT_CTRL(NN_P6tcl_T1B)] /* [0,65535] */;
    uint32_t NN_P6tcl_T2BM /* [1,50] */;
    uint32_t NN_P6tcl_T2BMI[50 EXT_STRUCT_CTRL(NN_P6tcl_T2BM)] /* [1,50] */;
    uint32_t NN_P6tcl_T2BME[50 EXT_STRUCT_CTRL(NN_P6tcl_T2BM)] /* [1,10000] */;
    uint32_t NN_P6tcl_T2B /* [0,10000] */;
    uint32_t NN_P6tcl_T2Bv[10000 EXT_STRUCT_CTRL(NN_P6tcl_T2B)] /* [0,65535] */;
    uint32_t NN_P6tfl_T1BM /* [1,50] */;
    uint32_t NN_P6tfl_T1BMI[50 EXT_STRUCT_CTRL(NN_P6tfl_T1BM)] /* [1,50] */;
    uint32_t NN_P6tfl_T1BME[50 EXT_STRUCT_CTRL(NN_P6tfl_T1BM)] /* [1,10000] */;
    uint32_t NN_P6tfl_T1B /* [0,10000] */;
    uint32_t NN_P6tfl_T1Bv[10000 EXT_STRUCT_CTRL(NN_P6tfl_T1B)] /* [0,65535] */;
    uint32_t NN_P6tfl_T2BM /* [1,50] */;
    uint32_t NN_P6tfl_T2BMI[50 EXT_STRUCT_CTRL(NN_P6tfl_T2BM)] /* [1,50] */;
    uint32_t NN_P6tfl_T2BME[50 EXT_STRUCT_CTRL(NN_P6tfl_T2BM)] /* [1,10000] */;
    uint32_t NN_P6tfl_T2B /* [0,10000] */;
    uint32_t NN_P6tfl_T2Bv[10000 EXT_STRUCT_CTRL(NN_P6tfl_T2B)] /* [0,65535] */;
    uint32_t NN_P6tct_T1BM /* [1,50] */;
    uint32_t NN_P6tct_T1BMI[50 EXT_STRUCT_CTRL(NN_P6tct_T1BM)] /* [1,50] */;
    uint32_t NN_P6tct_T1BME[50 EXT_STRUCT_CTRL(NN_P6tct_T1BM)] /* [1,10000] */;
    uint32_t NN_P6tct_T1B /* [0,10000] */;
    uint32_t NN_P6tct_T1Bv[10000 EXT_STRUCT_CTRL(NN_P6tct_T1B)] /* [0,65535] */;
    uint32_t NN_P6tct_T2BM /* [1,50] */;
    uint32_t NN_P6tct_T2BMI[50 EXT_STRUCT_CTRL(NN_P6tct_T2BM)] /* [1,50] */;
    uint32_t NN_P6tct_T2BME[50 EXT_STRUCT_CTRL(NN_P6tct_T2BM)] /* [1,10000] */;
    uint32_t NN_P6tct_T2B /* [0,10000] */;
    uint32_t NN_P6tct_T2Bv[10000 EXT_STRUCT_CTRL(NN_P6tct_T2B)] /* [0,65535] */;
    uint32_t NN_P6tft_T1BM /* [1,50] */;
    uint32_t NN_P6tft_T1BMI[50 EXT_STRUCT_CTRL(NN_P6tft_T1BM)] /* [1,50] */;
    uint32_t NN_P6tft_T1BME[50 EXT_STRUCT_CTRL(NN_P6tft_T1BM)] /* [1,10000] */;
    uint32_t NN_P6tft_T1B /* [0,10000] */;
    uint32_t NN_P6tft_T1Bv[10000 EXT_STRUCT_CTRL(NN_P6tft_T1B)] /* [0,65535] */;
    uint32_t NN_P6tft_T2BM /* [1,50] */;
    uint32_t NN_P6tft_T2BMI[50 EXT_STRUCT_CTRL(NN_P6tft_T2BM)] /* [1,50] */;
    uint32_t NN_P6tft_T2BME[50 EXT_STRUCT_CTRL(NN_P6tft_T2BM)] /* [1,10000] */;
    uint32_t NN_P6tft_T2B /* [0,10000] */;
    uint32_t NN_P6tft_T2Bv[10000 EXT_STRUCT_CTRL(NN_P6tft_T2B)] /* [0,65535] */;
    uint32_t NN_P7tcl_T1BM /* [1,50] */;
    uint32_t NN_P7tcl_T1BMI[50 EXT_STRUCT_CTRL(NN_P7tcl_T1BM)] /* [1,50] */;
    uint32_t NN_P7tcl_T1BME[50 EXT_STRUCT_CTRL(NN_P7tcl_T1BM)] /* [1,10000] */;
    uint32_t NN_P7tcl_T1B /* [0,10000] */;
    uint32_t NN_P7tcl_T1Bv[10000 EXT_STRUCT_CTRL(NN_P7tcl_T1B)] /* [0,65535] */;
    uint32_t NN_P7tcl_T2BM /* [1,50] */;
    uint32_t NN_P7tcl_T2BMI[50 EXT_STRUCT_CTRL(NN_P7tcl_T2BM)] /* [1,50] */;
    uint32_t NN_P7tcl_T2BME[50 EXT_STRUCT_CTRL(NN_P7tcl_T2BM)] /* [1,10000] */;
    uint32_t NN_P7tcl_T2B /* [0,10000] */;
    uint32_t NN_P7tcl_T2Bv[10000 EXT_STRUCT_CTRL(NN_P7tcl_T2B)] /* [0,65535] */;
    uint32_t NN_P7tfl_T1BM /* [1,50] */;
    uint32_t NN_P7tfl_T1BMI[50 EXT_STRUCT_CTRL(NN_P7tfl_T1BM)] /* [1,50] */;
    uint32_t NN_P7tfl_T1BME[50 EXT_STRUCT_CTRL(NN_P7tfl_T1BM)] /* [1,10000] */;
    uint32_t NN_P7tfl_T1B /* [0,10000] */;
    uint32_t NN_P7tfl_T1Bv[10000 EXT_STRUCT_CTRL(NN_P7tfl_T1B)] /* [0,65535] */;
    uint32_t NN_P7tfl_T2BM /* [1,50] */;
    uint32_t NN_P7tfl_T2BMI[50 EXT_STRUCT_CTRL(NN_P7tfl_T2BM)] /* [1,50] */;
    uint32_t NN_P7tfl_T2BME[50 EXT_STRUCT_CTRL(NN_P7tfl_T2BM)] /* [1,10000] */;
    uint32_t NN_P7tfl_T2B /* [0,10000] */;
    uint32_t NN_P7tfl_T2Bv[10000 EXT_STRUCT_CTRL(NN_P7tfl_T2B)] /* [0,65535] */;
    uint32_t NN_P7tct_T1BM /* [1,50] */;
    uint32_t NN_P7tct_T1BMI[50 EXT_STRUCT_CTRL(NN_P7tct_T1BM)] /* [1,50] */;
    uint32_t NN_P7tct_T1BME[50 EXT_STRUCT_CTRL(NN_P7tct_T1BM)] /* [1,10000] */;
    uint32_t NN_P7tct_T1B /* [0,10000] */;
    uint32_t NN_P7tct_T1Bv[10000 EXT_STRUCT_CTRL(NN_P7tct_T1B)] /* [0,65535] */;
    uint32_t NN_P7tct_T2BM /* [1,50] */;
    uint32_t NN_P7tct_T2BMI[50 EXT_STRUCT_CTRL(NN_P7tct_T2BM)] /* [1,50] */;
    uint32_t NN_P7tct_T2BME[50 EXT_STRUCT_CTRL(NN_P7tct_T2BM)] /* [1,10000] */;
    uint32_t NN_P7tct_T2B /* [0,10000] */;
    uint32_t NN_P7tct_T2Bv[10000 EXT_STRUCT_CTRL(NN_P7tct_T2B)] /* [0,65535] */;
    uint32_t NN_P7tft_T1BM /* [1,50] */;
    uint32_t NN_P7tft_T1BMI[50 EXT_STRUCT_CTRL(NN_P7tft_T1BM)] /* [1,50] */;
    uint32_t NN_P7tft_T1BME[50 EXT_STRUCT_CTRL(NN_P7tft_T1BM)] /* [1,10000] */;
    uint32_t NN_P7tft_T1B /* [0,10000] */;
    uint32_t NN_P7tft_T1Bv[10000 EXT_STRUCT_CTRL(NN_P7tft_T1B)] /* [0,65535] */;
    uint32_t NN_P7tft_T2BM /* [1,50] */;
    uint32_t NN_P7tft_T2BMI[50 EXT_STRUCT_CTRL(NN_P7tft_T2BM)] /* [1,50] */;
    uint32_t NN_P7tft_T2BME[50 EXT_STRUCT_CTRL(NN_P7tft_T2BM)] /* [1,10000] */;
    uint32_t NN_P7tft_T2B /* [0,10000] */;
    uint32_t NN_P7tft_T2Bv[10000 EXT_STRUCT_CTRL(NN_P7tft_T2B)] /* [0,65535] */;
    uint32_t NN_P8tcl_T1BM /* [1,50] */;
    uint32_t NN_P8tcl_T1BMI[50 EXT_STRUCT_CTRL(NN_P8tcl_T1BM)] /* [1,50] */;
    uint32_t NN_P8tcl_T1BME[50 EXT_STRUCT_CTRL(NN_P8tcl_T1BM)] /* [1,10000] */;
    uint32_t NN_P8tcl_T1B /* [0,10000] */;
    uint32_t NN_P8tcl_T1Bv[10000 EXT_STRUCT_CTRL(NN_P8tcl_T1B)] /* [0,65535] */;
    uint32_t NN_P8tcl_T2BM /* [1,50] */;
    uint32_t NN_P8tcl_T2BMI[50 EXT_STRUCT_CTRL(NN_P8tcl_T2BM)] /* [1,50] */;
    uint32_t NN_P8tcl_T2BME[50 EXT_STRUCT_CTRL(NN_P8tcl_T2BM)] /* [1,10000] */;
    uint32_t NN_P8tcl_T2B /* [0,10000] */;
    uint32_t NN_P8tcl_T2Bv[10000 EXT_STRUCT_CTRL(NN_P8tcl_T2B)] /* [0,65535] */;
    uint32_t NN_P8tfl_T1BM /* [1,50] */;
    uint32_t NN_P8tfl_T1BMI[50 EXT_STRUCT_CTRL(NN_P8tfl_T1BM)] /* [1,50] */;
    uint32_t NN_P8tfl_T1BME[50 EXT_STRUCT_CTRL(NN_P8tfl_T1BM)] /* [1,10000] */;
    uint32_t NN_P8tfl_T1B /* [0,10000] */;
    uint32_t NN_P8tfl_T1Bv[10000 EXT_STRUCT_CTRL(NN_P8tfl_T1B)] /* [0,65535] */;
    uint32_t NN_P8tfl_T2BM /* [1,50] */;
    uint32_t NN_P8tfl_T2BMI[50 EXT_STRUCT_CTRL(NN_P8tfl_T2BM)] /* [1,50] */;
    uint32_t NN_P8tfl_T2BME[50 EXT_STRUCT_CTRL(NN_P8tfl_T2BM)] /* [1,10000] */;
    uint32_t NN_P8tfl_T2B /* [0,10000] */;
    uint32_t NN_P8tfl_T2Bv[10000 EXT_STRUCT_CTRL(NN_P8tfl_T2B)] /* [0,65535] */;
    uint32_t NN_P8tct_T1BM /* [1,50] */;
    uint32_t NN_P8tct_T1BMI[50 EXT_STRUCT_CTRL(NN_P8tct_T1BM)] /* [1,50] */;
    uint32_t NN_P8tct_T1BME[50 EXT_STRUCT_CTRL(NN_P8tct_T1BM)] /* [1,10000] */;
    uint32_t NN_P8tct_T1B /* [0,10000] */;
    uint32_t NN_P8tct_T1Bv[10000 EXT_STRUCT_CTRL(NN_P8tct_T1B)] /* [0,65535] */;
    uint32_t NN_P8tct_T2BM /* [1,50] */;
    uint32_t NN_P8tct_T2BMI[50 EXT_STRUCT_CTRL(NN_P8tct_T2BM)] /* [1,50] */;
    uint32_t NN_P8tct_T2BME[50 EXT_STRUCT_CTRL(NN_P8tct_T2BM)] /* [1,10000] */;
    uint32_t NN_P8tct_T2B /* [0,10000] */;
    uint32_t NN_P8tct_T2Bv[10000 EXT_STRUCT_CTRL(NN_P8tct_T2B)] /* [0,65535] */;
    uint32_t NN_P8tft_T1BM /* [1,50] */;
    uint32_t NN_P8tft_T1BMI[50 EXT_STRUCT_CTRL(NN_P8tft_T1BM)] /* [1,50] */;
    uint32_t NN_P8tft_T1BME[50 EXT_STRUCT_CTRL(NN_P8tft_T1BM)] /* [1,10000] */;
    uint32_t NN_P8tft_T1B /* [0,10000] */;
    uint32_t NN_P8tft_T1Bv[10000 EXT_STRUCT_CTRL(NN_P8tft_T1B)] /* [0,65535] */;
    uint32_t NN_P8tft_T2BM /* [1,50] */;
    uint32_t NN_P8tft_T2BMI[50 EXT_STRUCT_CTRL(NN_P8tft_T2BM)] /* [1,50] */;
    uint32_t NN_P8tft_T2BME[50 EXT_STRUCT_CTRL(NN_P8tft_T2BM)] /* [1,10000] */;
    uint32_t NN_P8tft_T2B /* [0,10000] */;
    uint32_t NN_P8tft_T2Bv[10000 EXT_STRUCT_CTRL(NN_P8tft_T2B)] /* [0,65535] */;
    uint32_t NN_P9tcl_T1BM /* [1,50] */;
    uint32_t NN_P9tcl_T1BMI[50 EXT_STRUCT_CTRL(NN_P9tcl_T1BM)] /* [1,50] */;
    uint32_t NN_P9tcl_T1BME[50 EXT_STRUCT_CTRL(NN_P9tcl_T1BM)] /* [1,10000] */;
    uint32_t NN_P9tcl_T1B /* [0,10000] */;
    uint32_t NN_P9tcl_T1Bv[10000 EXT_STRUCT_CTRL(NN_P9tcl_T1B)] /* [0,65535] */;
    uint32_t NN_P9tcl_T2BM /* [1,50] */;
    uint32_t NN_P9tcl_T2BMI[50 EXT_STRUCT_CTRL(NN_P9tcl_T2BM)] /* [1,50] */;
    uint32_t NN_P9tcl_T2BME[50 EXT_STRUCT_CTRL(NN_P9tcl_T2BM)] /* [1,10000] */;
    uint32_t NN_P9tcl_T2B /* [0,10000] */;
    uint32_t NN_P9tcl_T2Bv[10000 EXT_STRUCT_CTRL(NN_P9tcl_T2B)] /* [0,65535] */;
    uint32_t NN_P9tfl_T1BM /* [1,50] */;
    uint32_t NN_P9tfl_T1BMI[50 EXT_STRUCT_CTRL(NN_P9tfl_T1BM)] /* [1,50] */;
    uint32_t NN_P9tfl_T1BME[50 EXT_STRUCT_CTRL(NN_P9tfl_T1BM)] /* [1,10000] */;
    uint32_t NN_P9tfl_T1B /* [0,10000] */;
    uint32_t NN_P9tfl_T1Bv[10000 EXT_STRUCT_CTRL(NN_P9tfl_T1B)] /* [0,65535] */;
    uint32_t NN_P9tfl_T2BM /* [1,50] */;
    uint32_t NN_P9tfl_T2BMI[50 EXT_STRUCT_CTRL(NN_P9tfl_T2BM)] /* [1,50] */;
    uint32_t NN_P9tfl_T2BME[50 EXT_STRUCT_CTRL(NN_P9tfl_T2BM)] /* [1,10000] */;
    uint32_t NN_P9tfl_T2B /* [0,10000] */;
    uint32_t NN_P9tfl_T2Bv[10000 EXT_STRUCT_CTRL(NN_P9tfl_T2B)] /* [0,65535] */;
    uint32_t NN_P9tct_T1BM /* [1,50] */;
    uint32_t NN_P9tct_T1BMI[50 EXT_STRUCT_CTRL(NN_P9tct_T1BM)] /* [1,50] */;
    uint32_t NN_P9tct_T1BME[50 EXT_STRUCT_CTRL(NN_P9tct_T1BM)] /* [1,10000] */;
    uint32_t NN_P9tct_T1B /* [0,10000] */;
    uint32_t NN_P9tct_T1Bv[10000 EXT_STRUCT_CTRL(NN_P9tct_T1B)] /* [0,65535] */;
    uint32_t NN_P9tct_T2BM /* [1,50] */;
    uint32_t NN_P9tct_T2BMI[50 EXT_STRUCT_CTRL(NN_P9tct_T2BM)] /* [1,50] */;
    uint32_t NN_P9tct_T2BME[50 EXT_STRUCT_CTRL(NN_P9tct_T2BM)] /* [1,10000] */;
    uint32_t NN_P9tct_T2B /* [0,10000] */;
    uint32_t NN_P9tct_T2Bv[10000 EXT_STRUCT_CTRL(NN_P9tct_T2B)] /* [0,65535] */;
    uint32_t NN_P9tft_T1BM /* [1,50] */;
    uint32_t NN_P9tft_T1BMI[50 EXT_STRUCT_CTRL(NN_P9tft_T1BM)] /* [1,50] */;
    uint32_t NN_P9tft_T1BME[50 EXT_STRUCT_CTRL(NN_P9tft_T1BM)] /* [1,10000] */;
    uint32_t NN_P9tft_T1B /* [0,10000] */;
    uint32_t NN_P9tft_T1Bv[10000 EXT_STRUCT_CTRL(NN_P9tft_T1B)] /* [0,65535] */;
    uint32_t NN_P9tft_T2BM /* [1,50] */;
    uint32_t NN_P9tft_T2BMI[50 EXT_STRUCT_CTRL(NN_P9tft_T2BM)] /* [1,50] */;
    uint32_t NN_P9tft_T2BME[50 EXT_STRUCT_CTRL(NN_P9tft_T2BM)] /* [1,10000] */;
    uint32_t NN_P9tft_T2B /* [0,10000] */;
    uint32_t NN_P9tft_T2Bv[10000 EXT_STRUCT_CTRL(NN_P9tft_T2B)] /* [0,65535] */;
    uint32_t NN_P10tcl_T1BM /* [1,50] */;
    uint32_t NN_P10tcl_T1BMI[50 EXT_STRUCT_CTRL(NN_P10tcl_T1BM)] /* [1,50] */;
    uint32_t NN_P10tcl_T1BME[50 EXT_STRUCT_CTRL(NN_P10tcl_T1BM)] /* [1,10000] */;
    uint32_t NN_P10tcl_T1B /* [0,10000] */;
    uint32_t NN_P10tcl_T1Bv[10000 EXT_STRUCT_CTRL(NN_P10tcl_T1B)] /* [0,65535] */;
    uint32_t NN_P10tcl_T2BM /* [1,50] */;
    uint32_t NN_P10tcl_T2BMI[50 EXT_STRUCT_CTRL(NN_P10tcl_T2BM)] /* [1,50] */;
    uint32_t NN_P10tcl_T2BME[50 EXT_STRUCT_CTRL(NN_P10tcl_T2BM)] /* [1,10000] */;
    uint32_t NN_P10tcl_T2B /* [0,10000] */;
    uint32_t NN_P10tcl_T2Bv[10000 EXT_STRUCT_CTRL(NN_P10tcl_T2B)] /* [0,65535] */;
    uint32_t NN_P10tfl_T1BM /* [1,50] */;
    uint32_t NN_P10tfl_T1BMI[50 EXT_STRUCT_CTRL(NN_P10tfl_T1BM)] /* [1,50] */;
    uint32_t NN_P10tfl_T1BME[50 EXT_STRUCT_CTRL(NN_P10tfl_T1BM)] /* [1,10000] */;
    uint32_t NN_P10tfl_T1B /* [0,10000] */;
    uint32_t NN_P10tfl_T1Bv[10000 EXT_STRUCT_CTRL(NN_P10tfl_T1B)] /* [0,65535] */;
    uint32_t NN_P10tfl_T2BM /* [1,50] */;
    uint32_t NN_P10tfl_T2BMI[50 EXT_STRUCT_CTRL(NN_P10tfl_T2BM)] /* [1,50] */;
    uint32_t NN_P10tfl_T2BME[50 EXT_STRUCT_CTRL(NN_P10tfl_T2BM)] /* [1,10000] */;
    uint32_t NN_P10tfl_T2B /* [0,10000] */;
    uint32_t NN_P10tfl_T2Bv[10000 EXT_STRUCT_CTRL(NN_P10tfl_T2B)] /* [0,65535] */;
    uint32_t NN_P10tct_T1BM /* [1,50] */;
    uint32_t NN_P10tct_T1BMI[50 EXT_STRUCT_CTRL(NN_P10tct_T1BM)] /* [1,50] */;
    uint32_t NN_P10tct_T1BME[50 EXT_STRUCT_CTRL(NN_P10tct_T1BM)] /* [1,10000] */;
    uint32_t NN_P10tct_T1B /* [0,10000] */;
    uint32_t NN_P10tct_T1Bv[10000 EXT_STRUCT_CTRL(NN_P10tct_T1B)] /* [0,65535] */;
    uint32_t NN_P10tct_T2BM /* [1,50] */;
    uint32_t NN_P10tct_T2BMI[50 EXT_STRUCT_CTRL(NN_P10tct_T2BM)] /* [1,50] */;
    uint32_t NN_P10tct_T2BME[50 EXT_STRUCT_CTRL(NN_P10tct_T2BM)] /* [1,10000] */;
    uint32_t NN_P10tct_T2B /* [0,10000] */;
    uint32_t NN_P10tct_T2Bv[10000 EXT_STRUCT_CTRL(NN_P10tct_T2B)] /* [0,65535] */;
    uint32_t NN_P10tft_T1BM /* [1,50] */;
    uint32_t NN_P10tft_T1BMI[50 EXT_STRUCT_CTRL(NN_P10tft_T1BM)] /* [1,50] */;
    uint32_t NN_P10tft_T1BME[50 EXT_STRUCT_CTRL(NN_P10tft_T1BM)] /* [1,10000] */;
    uint32_t NN_P10tft_T1B /* [0,10000] */;
    uint32_t NN_P10tft_T1Bv[10000 EXT_STRUCT_CTRL(NN_P10tft_T1B)] /* [0,65535] */;
    uint32_t NN_P10tft_T2BM /* [1,50] */;
    uint32_t NN_P10tft_T2BMI[50 EXT_STRUCT_CTRL(NN_P10tft_T2BM)] /* [1,50] */;
    uint32_t NN_P10tft_T2BME[50 EXT_STRUCT_CTRL(NN_P10tft_T2BM)] /* [1,10000] */;
    uint32_t NN_P10tft_T2B /* [0,10000] */;
    uint32_t NN_P10tft_T2Bv[10000 EXT_STRUCT_CTRL(NN_P10tft_T2B)] /* [0,65535] */;
    uint32_t NN_P11tcl_T1BM /* [1,50] */;
    uint32_t NN_P11tcl_T1BMI[50 EXT_STRUCT_CTRL(NN_P11tcl_T1BM)] /* [1,50] */;
    uint32_t NN_P11tcl_T1BME[50 EXT_STRUCT_CTRL(NN_P11tcl_T1BM)] /* [1,10000] */;
    uint32_t NN_P11tcl_T1B /* [0,10000] */;
    uint32_t NN_P11tcl_T1Bv[10000 EXT_STRUCT_CTRL(NN_P11tcl_T1B)] /* [0,65535] */;
    uint32_t NN_P11tcl_T2BM /* [1,50] */;
    uint32_t NN_P11tcl_T2BMI[50 EXT_STRUCT_CTRL(NN_P11tcl_T2BM)] /* [1,50] */;
    uint32_t NN_P11tcl_T2BME[50 EXT_STRUCT_CTRL(NN_P11tcl_T2BM)] /* [1,10000] */;
    uint32_t NN_P11tcl_T2B /* [0,10000] */;
    uint32_t NN_P11tcl_T2Bv[10000 EXT_STRUCT_CTRL(NN_P11tcl_T2B)] /* [0,65535] */;
    uint32_t NN_P11tfl_T1BM /* [1,50] */;
    uint32_t NN_P11tfl_T1BMI[50 EXT_STRUCT_CTRL(NN_P11tfl_T1BM)] /* [1,50] */;
    uint32_t NN_P11tfl_T1BME[50 EXT_STRUCT_CTRL(NN_P11tfl_T1BM)] /* [1,10000] */;
    uint32_t NN_P11tfl_T1B /* [0,10000] */;
    uint32_t NN_P11tfl_T1Bv[10000 EXT_STRUCT_CTRL(NN_P11tfl_T1B)] /* [0,65535] */;
    uint32_t NN_P11tfl_T2BM /* [1,50] */;
    uint32_t NN_P11tfl_T2BMI[50 EXT_STRUCT_CTRL(NN_P11tfl_T2BM)] /* [1,50] */;
    uint32_t NN_P11tfl_T2BME[50 EXT_STRUCT_CTRL(NN_P11tfl_T2BM)] /* [1,10000] */;
    uint32_t NN_P11tfl_T2B /* [0,10000] */;
    uint32_t NN_P11tfl_T2Bv[10000 EXT_STRUCT_CTRL(NN_P11tfl_T2B)] /* [0,65535] */;
    uint32_t NN_P11tct_T1BM /* [1,50] */;
    uint32_t NN_P11tct_T1BMI[50 EXT_STRUCT_CTRL(NN_P11tct_T1BM)] /* [1,50] */;
    uint32_t NN_P11tct_T1BME[50 EXT_STRUCT_CTRL(NN_P11tct_T1BM)] /* [1,10000] */;
    uint32_t NN_P11tct_T1B /* [0,10000] */;
    uint32_t NN_P11tct_T1Bv[10000 EXT_STRUCT_CTRL(NN_P11tct_T1B)] /* [0,65535] */;
    uint32_t NN_P11tct_T2BM /* [1,50] */;
    uint32_t NN_P11tct_T2BMI[50 EXT_STRUCT_CTRL(NN_P11tct_T2BM)] /* [1,50] */;
    uint32_t NN_P11tct_T2BME[50 EXT_STRUCT_CTRL(NN_P11tct_T2BM)] /* [1,10000] */;
    uint32_t NN_P11tct_T2B /* [0,10000] */;
    uint32_t NN_P11tct_T2Bv[10000 EXT_STRUCT_CTRL(NN_P11tct_T2B)] /* [0,65535] */;
    uint32_t NN_P11tft_T1BM /* [1,50] */;
    uint32_t NN_P11tft_T1BMI[50 EXT_STRUCT_CTRL(NN_P11tft_T1BM)] /* [1,50] */;
    uint32_t NN_P11tft_T1BME[50 EXT_STRUCT_CTRL(NN_P11tft_T1BM)] /* [1,10000] */;
    uint32_t NN_P11tft_T1B /* [0,10000] */;
    uint32_t NN_P11tft_T1Bv[10000 EXT_STRUCT_CTRL(NN_P11tft_T1B)] /* [0,65535] */;
    uint32_t NN_P11tft_T2BM /* [1,50] */;
    uint32_t NN_P11tft_T2BMI[50 EXT_STRUCT_CTRL(NN_P11tft_T2BM)] /* [1,50] */;
    uint32_t NN_P11tft_T2BME[50 EXT_STRUCT_CTRL(NN_P11tft_T2BM)] /* [1,10000] */;
    uint32_t NN_P11tft_T2B /* [0,10000] */;
    uint32_t NN_P11tft_T2Bv[10000 EXT_STRUCT_CTRL(NN_P11tft_T2B)] /* [0,65535] */;
    uint32_t NN_P12tcl_T1BM /* [1,50] */;
    uint32_t NN_P12tcl_T1BMI[50 EXT_STRUCT_CTRL(NN_P12tcl_T1BM)] /* [1,50] */;
    uint32_t NN_P12tcl_T1BME[50 EXT_STRUCT_CTRL(NN_P12tcl_T1BM)] /* [1,10000] */;
    uint32_t NN_P12tcl_T1B /* [0,10000] */;
    uint32_t NN_P12tcl_T1Bv[10000 EXT_STRUCT_CTRL(NN_P12tcl_T1B)] /* [0,65535] */;
    uint32_t NN_P12tcl_T2BM /* [1,50] */;
    uint32_t NN_P12tcl_T2BMI[50 EXT_STRUCT_CTRL(NN_P12tcl_T2BM)] /* [1,50] */;
    uint32_t NN_P12tcl_T2BME[50 EXT_STRUCT_CTRL(NN_P12tcl_T2BM)] /* [1,10000] */;
    uint32_t NN_P12tcl_T2B /* [0,10000] */;
    uint32_t NN_P12tcl_T2Bv[10000 EXT_STRUCT_CTRL(NN_P12tcl_T2B)] /* [0,65535] */;
    uint32_t NN_P12tfl_T1BM /* [1,50] */;
    uint32_t NN_P12tfl_T1BMI[50 EXT_STRUCT_CTRL(NN_P12tfl_T1BM)] /* [1,50] */;
    uint32_t NN_P12tfl_T1BME[50 EXT_STRUCT_CTRL(NN_P12tfl_T1BM)] /* [1,10000] */;
    uint32_t NN_P12tfl_T1B /* [0,10000] */;
    uint32_t NN_P12tfl_T1Bv[10000 EXT_STRUCT_CTRL(NN_P12tfl_T1B)] /* [0,65535] */;
    uint32_t NN_P12tfl_T2BM /* [1,50] */;
    uint32_t NN_P12tfl_T2BMI[50 EXT_STRUCT_CTRL(NN_P12tfl_T2BM)] /* [1,50] */;
    uint32_t NN_P12tfl_T2BME[50 EXT_STRUCT_CTRL(NN_P12tfl_T2BM)] /* [1,10000] */;
    uint32_t NN_P12tfl_T2B /* [0,10000] */;
    uint32_t NN_P12tfl_T2Bv[10000 EXT_STRUCT_CTRL(NN_P12tfl_T2B)] /* [0,65535] */;
    uint32_t NN_P12tct_T1BM /* [1,50] */;
    uint32_t NN_P12tct_T1BMI[50 EXT_STRUCT_CTRL(NN_P12tct_T1BM)] /* [1,50] */;
    uint32_t NN_P12tct_T1BME[50 EXT_STRUCT_CTRL(NN_P12tct_T1BM)] /* [1,10000] */;
    uint32_t NN_P12tct_T1B /* [0,10000] */;
    uint32_t NN_P12tct_T1Bv[10000 EXT_STRUCT_CTRL(NN_P12tct_T1B)] /* [0,65535] */;
    uint32_t NN_P12tct_T2BM /* [1,50] */;
    uint32_t NN_P12tct_T2BMI[50 EXT_STRUCT_CTRL(NN_P12tct_T2BM)] /* [1,50] */;
    uint32_t NN_P12tct_T2BME[50 EXT_STRUCT_CTRL(NN_P12tct_T2BM)] /* [1,10000] */;
    uint32_t NN_P12tct_T2B /* [0,10000] */;
    uint32_t NN_P12tct_T2Bv[10000 EXT_STRUCT_CTRL(NN_P12tct_T2B)] /* [0,65535] */;
    uint32_t NN_P12tft_T1BM /* [1,50] */;
    uint32_t NN_P12tft_T1BMI[50 EXT_STRUCT_CTRL(NN_P12tft_T1BM)] /* [1,50] */;
    uint32_t NN_P12tft_T1BME[50 EXT_STRUCT_CTRL(NN_P12tft_T1BM)] /* [1,10000] */;
    uint32_t NN_P12tft_T1B /* [0,10000] */;
    uint32_t NN_P12tft_T1Bv[10000 EXT_STRUCT_CTRL(NN_P12tft_T1B)] /* [0,65535] */;
    uint32_t NN_P12tft_T2BM /* [1,50] */;
    uint32_t NN_P12tft_T2BMI[50 EXT_STRUCT_CTRL(NN_P12tft_T2BM)] /* [1,50] */;
    uint32_t NN_P12tft_T2BME[50 EXT_STRUCT_CTRL(NN_P12tft_T2BM)] /* [1,10000] */;
    uint32_t NN_P12tft_T2B /* [0,10000] */;
    uint32_t NN_P12tft_T2Bv[10000 EXT_STRUCT_CTRL(NN_P12tft_T2B)] /* [0,65535] */;
    uint32_t NN_P13tcl_T1BM /* [1,50] */;
    uint32_t NN_P13tcl_T1BMI[50 EXT_STRUCT_CTRL(NN_P13tcl_T1BM)] /* [1,50] */;
    uint32_t NN_P13tcl_T1BME[50 EXT_STRUCT_CTRL(NN_P13tcl_T1BM)] /* [1,10000] */;
    uint32_t NN_P13tcl_T1B /* [0,10000] */;
    uint32_t NN_P13tcl_T1Bv[10000 EXT_STRUCT_CTRL(NN_P13tcl_T1B)] /* [0,65535] */;
    uint32_t NN_P13tcl_T2BM /* [1,50] */;
    uint32_t NN_P13tcl_T2BMI[50 EXT_STRUCT_CTRL(NN_P13tcl_T2BM)] /* [1,50] */;
    uint32_t NN_P13tcl_T2BME[50 EXT_STRUCT_CTRL(NN_P13tcl_T2BM)] /* [1,10000] */;
    uint32_t NN_P13tcl_T2B /* [0,10000] */;
    uint32_t NN_P13tcl_T2Bv[10000 EXT_STRUCT_CTRL(NN_P13tcl_T2B)] /* [0,65535] */;
    uint32_t NN_P13tfl_T1BM /* [1,50] */;
    uint32_t NN_P13tfl_T1BMI[50 EXT_STRUCT_CTRL(NN_P13tfl_T1BM)] /* [1,50] */;
    uint32_t NN_P13tfl_T1BME[50 EXT_STRUCT_CTRL(NN_P13tfl_T1BM)] /* [1,10000] */;
    uint32_t NN_P13tfl_T1B /* [0,10000] */;
    uint32_t NN_P13tfl_T1Bv[10000 EXT_STRUCT_CTRL(NN_P13tfl_T1B)] /* [0,65535] */;
    uint32_t NN_P13tfl_T2BM /* [1,50] */;
    uint32_t NN_P13tfl_T2BMI[50 EXT_STRUCT_CTRL(NN_P13tfl_T2BM)] /* [1,50] */;
    uint32_t NN_P13tfl_T2BME[50 EXT_STRUCT_CTRL(NN_P13tfl_T2BM)] /* [1,10000] */;
    uint32_t NN_P13tfl_T2B /* [0,10000] */;
    uint32_t NN_P13tfl_T2Bv[10000 EXT_STRUCT_CTRL(NN_P13tfl_T2B)] /* [0,65535] */;
    uint32_t NN_P13tct_T1BM /* [1,50] */;
    uint32_t NN_P13tct_T1BMI[50 EXT_STRUCT_CTRL(NN_P13tct_T1BM)] /* [1,50] */;
    uint32_t NN_P13tct_T1BME[50 EXT_STRUCT_CTRL(NN_P13tct_T1BM)] /* [1,10000] */;
    uint32_t NN_P13tct_T1B /* [0,10000] */;
    uint32_t NN_P13tct_T1Bv[10000 EXT_STRUCT_CTRL(NN_P13tct_T1B)] /* [0,65535] */;
    uint32_t NN_P13tct_T2BM /* [1,50] */;
    uint32_t NN_P13tct_T2BMI[50 EXT_STRUCT_CTRL(NN_P13tct_T2BM)] /* [1,50] */;
    uint32_t NN_P13tct_T2BME[50 EXT_STRUCT_CTRL(NN_P13tct_T2BM)] /* [1,10000] */;
    uint32_t NN_P13tct_T2B /* [0,10000] */;
    uint32_t NN_P13tct_T2Bv[10000 EXT_STRUCT_CTRL(NN_P13tct_T2B)] /* [0,65535] */;
    uint32_t NN_P13tft_T1BM /* [1,50] */;
    uint32_t NN_P13tft_T1BMI[50 EXT_STRUCT_CTRL(NN_P13tft_T1BM)] /* [1,50] */;
    uint32_t NN_P13tft_T1BME[50 EXT_STRUCT_CTRL(NN_P13tft_T1BM)] /* [1,10000] */;
    uint32_t NN_P13tft_T1B /* [0,10000] */;
    uint32_t NN_P13tft_T1Bv[10000 EXT_STRUCT_CTRL(NN_P13tft_T1B)] /* [0,65535] */;
    uint32_t NN_P13tft_T2BM /* [1,50] */;
    uint32_t NN_P13tft_T2BMI[50 EXT_STRUCT_CTRL(NN_P13tft_T2BM)] /* [1,50] */;
    uint32_t NN_P13tft_T2BME[50 EXT_STRUCT_CTRL(NN_P13tft_T2BM)] /* [1,10000] */;
    uint32_t NN_P13tft_T2B /* [0,10000] */;
    uint32_t NN_P13tft_T2Bv[10000 EXT_STRUCT_CTRL(NN_P13tft_T2B)] /* [0,65535] */;
    uint32_t NN_P14tcl_T1BM /* [1,50] */;
    uint32_t NN_P14tcl_T1BMI[50 EXT_STRUCT_CTRL(NN_P14tcl_T1BM)] /* [1,50] */;
    uint32_t NN_P14tcl_T1BME[50 EXT_STRUCT_CTRL(NN_P14tcl_T1BM)] /* [1,10000] */;
    uint32_t NN_P14tcl_T1B /* [0,10000] */;
    uint32_t NN_P14tcl_T1Bv[10000 EXT_STRUCT_CTRL(NN_P14tcl_T1B)] /* [0,65535] */;
    uint32_t NN_P14tcl_T2BM /* [1,50] */;
    uint32_t NN_P14tcl_T2BMI[50 EXT_STRUCT_CTRL(NN_P14tcl_T2BM)] /* [1,50] */;
    uint32_t NN_P14tcl_T2BME[50 EXT_STRUCT_CTRL(NN_P14tcl_T2BM)] /* [1,10000] */;
    uint32_t NN_P14tcl_T2B /* [0,10000] */;
    uint32_t NN_P14tcl_T2Bv[10000 EXT_STRUCT_CTRL(NN_P14tcl_T2B)] /* [0,65535] */;
    uint32_t NN_P14tfl_T1BM /* [1,50] */;
    uint32_t NN_P14tfl_T1BMI[50 EXT_STRUCT_CTRL(NN_P14tfl_T1BM)] /* [1,50] */;
    uint32_t NN_P14tfl_T1BME[50 EXT_STRUCT_CTRL(NN_P14tfl_T1BM)] /* [1,10000] */;
    uint32_t NN_P14tfl_T1B /* [0,10000] */;
    uint32_t NN_P14tfl_T1Bv[10000 EXT_STRUCT_CTRL(NN_P14tfl_T1B)] /* [0,65535] */;
    uint32_t NN_P14tfl_T2BM /* [1,50] */;
    uint32_t NN_P14tfl_T2BMI[50 EXT_STRUCT_CTRL(NN_P14tfl_T2BM)] /* [1,50] */;
    uint32_t NN_P14tfl_T2BME[50 EXT_STRUCT_CTRL(NN_P14tfl_T2BM)] /* [1,10000] */;
    uint32_t NN_P14tfl_T2B /* [0,10000] */;
    uint32_t NN_P14tfl_T2Bv[10000 EXT_STRUCT_CTRL(NN_P14tfl_T2B)] /* [0,65535] */;
    uint32_t NN_P14tct_T1BM /* [1,50] */;
    uint32_t NN_P14tct_T1BMI[50 EXT_STRUCT_CTRL(NN_P14tct_T1BM)] /* [1,50] */;
    uint32_t NN_P14tct_T1BME[50 EXT_STRUCT_CTRL(NN_P14tct_T1BM)] /* [1,10000] */;
    uint32_t NN_P14tct_T1B /* [0,10000] */;
    uint32_t NN_P14tct_T1Bv[10000 EXT_STRUCT_CTRL(NN_P14tct_T1B)] /* [0,65535] */;
    uint32_t NN_P14tct_T2BM /* [1,50] */;
    uint32_t NN_P14tct_T2BMI[50 EXT_STRUCT_CTRL(NN_P14tct_T2BM)] /* [1,50] */;
    uint32_t NN_P14tct_T2BME[50 EXT_STRUCT_CTRL(NN_P14tct_T2BM)] /* [1,10000] */;
    uint32_t NN_P14tct_T2B /* [0,10000] */;
    uint32_t NN_P14tct_T2Bv[10000 EXT_STRUCT_CTRL(NN_P14tct_T2B)] /* [0,65535] */;
    uint32_t NN_P14tft_T1BM /* [1,50] */;
    uint32_t NN_P14tft_T1BMI[50 EXT_STRUCT_CTRL(NN_P14tft_T1BM)] /* [1,50] */;
    uint32_t NN_P14tft_T1BME[50 EXT_STRUCT_CTRL(NN_P14tft_T1BM)] /* [1,10000] */;
    uint32_t NN_P14tft_T1B /* [0,10000] */;
    uint32_t NN_P14tft_T1Bv[10000 EXT_STRUCT_CTRL(NN_P14tft_T1B)] /* [0,65535] */;
    uint32_t NN_P14tft_T2BM /* [1,50] */;
    uint32_t NN_P14tft_T2BMI[50 EXT_STRUCT_CTRL(NN_P14tft_T2BM)] /* [1,50] */;
    uint32_t NN_P14tft_T2BME[50 EXT_STRUCT_CTRL(NN_P14tft_T2BM)] /* [1,10000] */;
    uint32_t NN_P14tft_T2B /* [0,10000] */;
    uint32_t NN_P14tft_T2Bv[10000 EXT_STRUCT_CTRL(NN_P14tft_T2B)] /* [0,65535] */;
    uint32_t NN_P15tcl_T1BM /* [1,50] */;
    uint32_t NN_P15tcl_T1BMI[50 EXT_STRUCT_CTRL(NN_P15tcl_T1BM)] /* [1,50] */;
    uint32_t NN_P15tcl_T1BME[50 EXT_STRUCT_CTRL(NN_P15tcl_T1BM)] /* [1,10000] */;
    uint32_t NN_P15tcl_T1B /* [0,10000] */;
    uint32_t NN_P15tcl_T1Bv[10000 EXT_STRUCT_CTRL(NN_P15tcl_T1B)] /* [0,65535] */;
    uint32_t NN_P15tcl_T2BM /* [1,50] */;
    uint32_t NN_P15tcl_T2BMI[50 EXT_STRUCT_CTRL(NN_P15tcl_T2BM)] /* [1,50] */;
    uint32_t NN_P15tcl_T2BME[50 EXT_STRUCT_CTRL(NN_P15tcl_T2BM)] /* [1,10000] */;
    uint32_t NN_P15tcl_T2B /* [0,10000] */;
    uint32_t NN_P15tcl_T2Bv[10000 EXT_STRUCT_CTRL(NN_P15tcl_T2B)] /* [0,65535] */;
    uint32_t NN_P15tfl_T1BM /* [1,50] */;
    uint32_t NN_P15tfl_T1BMI[50 EXT_STRUCT_CTRL(NN_P15tfl_T1BM)] /* [1,50] */;
    uint32_t NN_P15tfl_T1BME[50 EXT_STRUCT_CTRL(NN_P15tfl_T1BM)] /* [1,10000] */;
    uint32_t NN_P15tfl_T1B /* [0,10000] */;
    uint32_t NN_P15tfl_T1Bv[10000 EXT_STRUCT_CTRL(NN_P15tfl_T1B)] /* [0,65535] */;
    uint32_t NN_P15tfl_T2BM /* [1,50] */;
    uint32_t NN_P15tfl_T2BMI[50 EXT_STRUCT_CTRL(NN_P15tfl_T2BM)] /* [1,50] */;
    uint32_t NN_P15tfl_T2BME[50 EXT_STRUCT_CTRL(NN_P15tfl_T2BM)] /* [1,10000] */;
    uint32_t NN_P15tfl_T2B /* [0,10000] */;
    uint32_t NN_P15tfl_T2Bv[10000 EXT_STRUCT_CTRL(NN_P15tfl_T2B)] /* [0,65535] */;
    uint32_t NN_P15tct_T1BM /* [1,50] */;
    uint32_t NN_P15tct_T1BMI[50 EXT_STRUCT_CTRL(NN_P15tct_T1BM)] /* [1,50] */;
    uint32_t NN_P15tct_T1BME[50 EXT_STRUCT_CTRL(NN_P15tct_T1BM)] /* [1,10000] */;
    uint32_t NN_P15tct_T1B /* [0,10000] */;
    uint32_t NN_P15tct_T1Bv[10000 EXT_STRUCT_CTRL(NN_P15tct_T1B)] /* [0,65535] */;
    uint32_t NN_P15tct_T2BM /* [1,50] */;
    uint32_t NN_P15tct_T2BMI[50 EXT_STRUCT_CTRL(NN_P15tct_T2BM)] /* [1,50] */;
    uint32_t NN_P15tct_T2BME[50 EXT_STRUCT_CTRL(NN_P15tct_T2BM)] /* [1,10000] */;
    uint32_t NN_P15tct_T2B /* [0,10000] */;
    uint32_t NN_P15tct_T2Bv[10000 EXT_STRUCT_CTRL(NN_P15tct_T2B)] /* [0,65535] */;
    uint32_t NN_P15tft_T1BM /* [1,50] */;
    uint32_t NN_P15tft_T1BMI[50 EXT_STRUCT_CTRL(NN_P15tft_T1BM)] /* [1,50] */;
    uint32_t NN_P15tft_T1BME[50 EXT_STRUCT_CTRL(NN_P15tft_T1BM)] /* [1,10000] */;
    uint32_t NN_P15tft_T1B /* [0,10000] */;
    uint32_t NN_P15tft_T1Bv[10000 EXT_STRUCT_CTRL(NN_P15tft_T1B)] /* [0,65535] */;
    uint32_t NN_P15tft_T2BM /* [1,50] */;
    uint32_t NN_P15tft_T2BMI[50 EXT_STRUCT_CTRL(NN_P15tft_T2BM)] /* [1,50] */;
    uint32_t NN_P15tft_T2BME[50 EXT_STRUCT_CTRL(NN_P15tft_T2BM)] /* [1,10000] */;
    uint32_t NN_P15tft_T2B /* [0,10000] */;
    uint32_t NN_P15tft_T2Bv[10000 EXT_STRUCT_CTRL(NN_P15tft_T2B)] /* [0,65535] */;
    uint32_t NN_P16tcl_T1BM /* [1,50] */;
    uint32_t NN_P16tcl_T1BMI[50 EXT_STRUCT_CTRL(NN_P16tcl_T1BM)] /* [1,50] */;
    uint32_t NN_P16tcl_T1BME[50 EXT_STRUCT_CTRL(NN_P16tcl_T1BM)] /* [1,10000] */;
    uint32_t NN_P16tcl_T1B /* [0,10000] */;
    uint32_t NN_P16tcl_T1Bv[10000 EXT_STRUCT_CTRL(NN_P16tcl_T1B)] /* [0,65535] */;
    uint32_t NN_P16tcl_T2BM /* [1,50] */;
    uint32_t NN_P16tcl_T2BMI[50 EXT_STRUCT_CTRL(NN_P16tcl_T2BM)] /* [1,50] */;
    uint32_t NN_P16tcl_T2BME[50 EXT_STRUCT_CTRL(NN_P16tcl_T2BM)] /* [1,10000] */;
    uint32_t NN_P16tcl_T2B /* [0,10000] */;
    uint32_t NN_P16tcl_T2Bv[10000 EXT_STRUCT_CTRL(NN_P16tcl_T2B)] /* [0,65535] */;
    uint32_t NN_P16tfl_T1BM /* [1,50] */;
    uint32_t NN_P16tfl_T1BMI[50 EXT_STRUCT_CTRL(NN_P16tfl_T1BM)] /* [1,50] */;
    uint32_t NN_P16tfl_T1BME[50 EXT_STRUCT_CTRL(NN_P16tfl_T1BM)] /* [1,10000] */;
    uint32_t NN_P16tfl_T1B /* [0,10000] */;
    uint32_t NN_P16tfl_T1Bv[10000 EXT_STRUCT_CTRL(NN_P16tfl_T1B)] /* [0,65535] */;
    uint32_t NN_P16tfl_T2BM /* [1,50] */;
    uint32_t NN_P16tfl_T2BMI[50 EXT_STRUCT_CTRL(NN_P16tfl_T2BM)] /* [1,50] */;
    uint32_t NN_P16tfl_T2BME[50 EXT_STRUCT_CTRL(NN_P16tfl_T2BM)] /* [1,10000] */;
    uint32_t NN_P16tfl_T2B /* [0,10000] */;
    uint32_t NN_P16tfl_T2Bv[10000 EXT_STRUCT_CTRL(NN_P16tfl_T2B)] /* [0,65535] */;
    uint32_t NN_P16tct_T1BM /* [1,50] */;
    uint32_t NN_P16tct_T1BMI[50 EXT_STRUCT_CTRL(NN_P16tct_T1BM)] /* [1,50] */;
    uint32_t NN_P16tct_T1BME[50 EXT_STRUCT_CTRL(NN_P16tct_T1BM)] /* [1,10000] */;
    uint32_t NN_P16tct_T1B /* [0,10000] */;
    uint32_t NN_P16tct_T1Bv[10000 EXT_STRUCT_CTRL(NN_P16tct_T1B)] /* [0,65535] */;
    uint32_t NN_P16tct_T2BM /* [1,50] */;
    uint32_t NN_P16tct_T2BMI[50 EXT_STRUCT_CTRL(NN_P16tct_T2BM)] /* [1,50] */;
    uint32_t NN_P16tct_T2BME[50 EXT_STRUCT_CTRL(NN_P16tct_T2BM)] /* [1,10000] */;
    uint32_t NN_P16tct_T2B /* [0,10000] */;
    uint32_t NN_P16tct_T2Bv[10000 EXT_STRUCT_CTRL(NN_P16tct_T2B)] /* [0,65535] */;
    uint32_t NN_P16tft_T1BM /* [1,50] */;
    uint32_t NN_P16tft_T1BMI[50 EXT_STRUCT_CTRL(NN_P16tft_T1BM)] /* [1,50] */;
    uint32_t NN_P16tft_T1BME[50 EXT_STRUCT_CTRL(NN_P16tft_T1BM)] /* [1,10000] */;
    uint32_t NN_P16tft_T1B /* [0,10000] */;
    uint32_t NN_P16tft_T1Bv[10000 EXT_STRUCT_CTRL(NN_P16tft_T1B)] /* [0,65535] */;
    uint32_t NN_P16tft_T2BM /* [1,50] */;
    uint32_t NN_P16tft_T2BMI[50 EXT_STRUCT_CTRL(NN_P16tft_T2BM)] /* [1,50] */;
    uint32_t NN_P16tft_T2BME[50 EXT_STRUCT_CTRL(NN_P16tft_T2BM)] /* [1,10000] */;
    uint32_t NN_P16tft_T2B /* [0,10000] */;
    uint32_t NN_P16tft_T2Bv[10000 EXT_STRUCT_CTRL(NN_P16tft_T2B)] /* [0,65535] */;
    uint32_t NN_P17tcl_T1BM /* [1,50] */;
    uint32_t NN_P17tcl_T1BMI[50 EXT_STRUCT_CTRL(NN_P17tcl_T1BM)] /* [1,50] */;
    uint32_t NN_P17tcl_T1BME[50 EXT_STRUCT_CTRL(NN_P17tcl_T1BM)] /* [1,10000] */;
    uint32_t NN_P17tcl_T1B /* [0,10000] */;
    uint32_t NN_P17tcl_T1Bv[10000 EXT_STRUCT_CTRL(NN_P17tcl_T1B)] /* [0,65535] */;
    uint32_t NN_P17tcl_T2BM /* [1,50] */;
    uint32_t NN_P17tcl_T2BMI[50 EXT_STRUCT_CTRL(NN_P17tcl_T2BM)] /* [1,50] */;
    uint32_t NN_P17tcl_T2BME[50 EXT_STRUCT_CTRL(NN_P17tcl_T2BM)] /* [1,10000] */;
    uint32_t NN_P17tcl_T2B /* [0,10000] */;
    uint32_t NN_P17tcl_T2Bv[10000 EXT_STRUCT_CTRL(NN_P17tcl_T2B)] /* [0,65535] */;
    uint32_t NN_P17tfl_T1BM /* [1,50] */;
    uint32_t NN_P17tfl_T1BMI[50 EXT_STRUCT_CTRL(NN_P17tfl_T1BM)] /* [1,50] */;
    uint32_t NN_P17tfl_T1BME[50 EXT_STRUCT_CTRL(NN_P17tfl_T1BM)] /* [1,10000] */;
    uint32_t NN_P17tfl_T1B /* [0,10000] */;
    uint32_t NN_P17tfl_T1Bv[10000 EXT_STRUCT_CTRL(NN_P17tfl_T1B)] /* [0,65535] */;
    uint32_t NN_P17tfl_T2BM /* [1,50] */;
    uint32_t NN_P17tfl_T2BMI[50 EXT_STRUCT_CTRL(NN_P17tfl_T2BM)] /* [1,50] */;
    uint32_t NN_P17tfl_T2BME[50 EXT_STRUCT_CTRL(NN_P17tfl_T2BM)] /* [1,10000] */;
    uint32_t NN_P17tfl_T2B /* [0,10000] */;
    uint32_t NN_P17tfl_T2Bv[10000 EXT_STRUCT_CTRL(NN_P17tfl_T2B)] /* [0,65535] */;
    uint32_t NN_P17tct_T1BM /* [1,50] */;
    uint32_t NN_P17tct_T1BMI[50 EXT_STRUCT_CTRL(NN_P17tct_T1BM)] /* [1,50] */;
    uint32_t NN_P17tct_T1BME[50 EXT_STRUCT_CTRL(NN_P17tct_T1BM)] /* [1,10000] */;
    uint32_t NN_P17tct_T1B /* [0,10000] */;
    uint32_t NN_P17tct_T1Bv[10000 EXT_STRUCT_CTRL(NN_P17tct_T1B)] /* [0,65535] */;
    uint32_t NN_P17tct_T2BM /* [1,50] */;
    uint32_t NN_P17tct_T2BMI[50 EXT_STRUCT_CTRL(NN_P17tct_T2BM)] /* [1,50] */;
    uint32_t NN_P17tct_T2BME[50 EXT_STRUCT_CTRL(NN_P17tct_T2BM)] /* [1,10000] */;
    uint32_t NN_P17tct_T2B /* [0,10000] */;
    uint32_t NN_P17tct_T2Bv[10000 EXT_STRUCT_CTRL(NN_P17tct_T2B)] /* [0,65535] */;
    uint32_t NN_P17tft_T1BM /* [1,50] */;
    uint32_t NN_P17tft_T1BMI[50 EXT_STRUCT_CTRL(NN_P17tft_T1BM)] /* [1,50] */;
    uint32_t NN_P17tft_T1BME[50 EXT_STRUCT_CTRL(NN_P17tft_T1BM)] /* [1,10000] */;
    uint32_t NN_P17tft_T1B /* [0,10000] */;
    uint32_t NN_P17tft_T1Bv[10000 EXT_STRUCT_CTRL(NN_P17tft_T1B)] /* [0,65535] */;
    uint32_t NN_P17tft_T2BM /* [1,50] */;
    uint32_t NN_P17tft_T2BMI[50 EXT_STRUCT_CTRL(NN_P17tft_T2BM)] /* [1,50] */;
    uint32_t NN_P17tft_T2BME[50 EXT_STRUCT_CTRL(NN_P17tft_T2BM)] /* [1,10000] */;
    uint32_t NN_P17tft_T2B /* [0,10000] */;
    uint32_t NN_P17tft_T2Bv[10000 EXT_STRUCT_CTRL(NN_P17tft_T2B)] /* [0,65535] */;
    uint32_t NN_P18tcl_T1BM /* [1,50] */;
    uint32_t NN_P18tcl_T1BMI[50 EXT_STRUCT_CTRL(NN_P18tcl_T1BM)] /* [1,50] */;
    uint32_t NN_P18tcl_T1BME[50 EXT_STRUCT_CTRL(NN_P18tcl_T1BM)] /* [1,10000] */;
    uint32_t NN_P18tcl_T1B /* [0,10000] */;
    uint32_t NN_P18tcl_T1Bv[10000 EXT_STRUCT_CTRL(NN_P18tcl_T1B)] /* [0,65535] */;
    uint32_t NN_P18tcl_T2BM /* [1,50] */;
    uint32_t NN_P18tcl_T2BMI[50 EXT_STRUCT_CTRL(NN_P18tcl_T2BM)] /* [1,50] */;
    uint32_t NN_P18tcl_T2BME[50 EXT_STRUCT_CTRL(NN_P18tcl_T2BM)] /* [1,10000] */;
    uint32_t NN_P18tcl_T2B /* [0,10000] */;
    uint32_t NN_P18tcl_T2Bv[10000 EXT_STRUCT_CTRL(NN_P18tcl_T2B)] /* [0,65535] */;
    uint32_t NN_P18tfl_T1BM /* [1,50] */;
    uint32_t NN_P18tfl_T1BMI[50 EXT_STRUCT_CTRL(NN_P18tfl_T1BM)] /* [1,50] */;
    uint32_t NN_P18tfl_T1BME[50 EXT_STRUCT_CTRL(NN_P18tfl_T1BM)] /* [1,10000] */;
    uint32_t NN_P18tfl_T1B /* [0,10000] */;
    uint32_t NN_P18tfl_T1Bv[10000 EXT_STRUCT_CTRL(NN_P18tfl_T1B)] /* [0,65535] */;
    uint32_t NN_P18tfl_T2BM /* [1,50] */;
    uint32_t NN_P18tfl_T2BMI[50 EXT_STRUCT_CTRL(NN_P18tfl_T2BM)] /* [1,50] */;
    uint32_t NN_P18tfl_T2BME[50 EXT_STRUCT_CTRL(NN_P18tfl_T2BM)] /* [1,10000] */;
    uint32_t NN_P18tfl_T2B /* [0,10000] */;
    uint32_t NN_P18tfl_T2Bv[10000 EXT_STRUCT_CTRL(NN_P18tfl_T2B)] /* [0,65535] */;
    uint32_t NN_P18tct_T1BM /* [1,50] */;
    uint32_t NN_P18tct_T1BMI[50 EXT_STRUCT_CTRL(NN_P18tct_T1BM)] /* [1,50] */;
    uint32_t NN_P18tct_T1BME[50 EXT_STRUCT_CTRL(NN_P18tct_T1BM)] /* [1,10000] */;
    uint32_t NN_P18tct_T1B /* [0,10000] */;
    uint32_t NN_P18tct_T1Bv[10000 EXT_STRUCT_CTRL(NN_P18tct_T1B)] /* [0,65535] */;
    uint32_t NN_P18tct_T2BM /* [1,50] */;
    uint32_t NN_P18tct_T2BMI[50 EXT_STRUCT_CTRL(NN_P18tct_T2BM)] /* [1,50] */;
    uint32_t NN_P18tct_T2BME[50 EXT_STRUCT_CTRL(NN_P18tct_T2BM)] /* [1,10000] */;
    uint32_t NN_P18tct_T2B /* [0,10000] */;
    uint32_t NN_P18tct_T2Bv[10000 EXT_STRUCT_CTRL(NN_P18tct_T2B)] /* [0,65535] */;
    uint32_t NN_P18tft_T1BM /* [1,50] */;
    uint32_t NN_P18tft_T1BMI[50 EXT_STRUCT_CTRL(NN_P18tft_T1BM)] /* [1,50] */;
    uint32_t NN_P18tft_T1BME[50 EXT_STRUCT_CTRL(NN_P18tft_T1BM)] /* [1,10000] */;
    uint32_t NN_P18tft_T1B /* [0,10000] */;
    uint32_t NN_P18tft_T1Bv[10000 EXT_STRUCT_CTRL(NN_P18tft_T1B)] /* [0,65535] */;
    uint32_t NN_P18tft_T2BM /* [1,50] */;
    uint32_t NN_P18tft_T2BMI[50 EXT_STRUCT_CTRL(NN_P18tft_T2BM)] /* [1,50] */;
    uint32_t NN_P18tft_T2BME[50 EXT_STRUCT_CTRL(NN_P18tft_T2BM)] /* [1,10000] */;
    uint32_t NN_P18tft_T2B /* [0,10000] */;
    uint32_t NN_P18tft_T2Bv[10000 EXT_STRUCT_CTRL(NN_P18tft_T2B)] /* [0,65535] */;
    uint32_t NN_P19tcl_T1BM /* [1,50] */;
    uint32_t NN_P19tcl_T1BMI[50 EXT_STRUCT_CTRL(NN_P19tcl_T1BM)] /* [1,50] */;
    uint32_t NN_P19tcl_T1BME[50 EXT_STRUCT_CTRL(NN_P19tcl_T1BM)] /* [1,10000] */;
    uint32_t NN_P19tcl_T1B /* [0,10000] */;
    uint32_t NN_P19tcl_T1Bv[10000 EXT_STRUCT_CTRL(NN_P19tcl_T1B)] /* [0,65535] */;
    uint32_t NN_P19tcl_T2BM /* [1,50] */;
    uint32_t NN_P19tcl_T2BMI[50 EXT_STRUCT_CTRL(NN_P19tcl_T2BM)] /* [1,50] */;
    uint32_t NN_P19tcl_T2BME[50 EXT_STRUCT_CTRL(NN_P19tcl_T2BM)] /* [1,10000] */;
    uint32_t NN_P19tcl_T2B /* [0,10000] */;
    uint32_t NN_P19tcl_T2Bv[10000 EXT_STRUCT_CTRL(NN_P19tcl_T2B)] /* [0,65535] */;
    uint32_t NN_P19tfl_T1BM /* [1,50] */;
    uint32_t NN_P19tfl_T1BMI[50 EXT_STRUCT_CTRL(NN_P19tfl_T1BM)] /* [1,50] */;
    uint32_t NN_P19tfl_T1BME[50 EXT_STRUCT_CTRL(NN_P19tfl_T1BM)] /* [1,10000] */;
    uint32_t NN_P19tfl_T1B /* [0,10000] */;
    uint32_t NN_P19tfl_T1Bv[10000 EXT_STRUCT_CTRL(NN_P19tfl_T1B)] /* [0,65535] */;
    uint32_t NN_P19tfl_T2BM /* [1,50] */;
    uint32_t NN_P19tfl_T2BMI[50 EXT_STRUCT_CTRL(NN_P19tfl_T2BM)] /* [1,50] */;
    uint32_t NN_P19tfl_T2BME[50 EXT_STRUCT_CTRL(NN_P19tfl_T2BM)] /* [1,10000] */;
    uint32_t NN_P19tfl_T2B /* [0,10000] */;
    uint32_t NN_P19tfl_T2Bv[10000 EXT_STRUCT_CTRL(NN_P19tfl_T2B)] /* [0,65535] */;
    uint32_t NN_P19tct_T1BM /* [1,50] */;
    uint32_t NN_P19tct_T1BMI[50 EXT_STRUCT_CTRL(NN_P19tct_T1BM)] /* [1,50] */;
    uint32_t NN_P19tct_T1BME[50 EXT_STRUCT_CTRL(NN_P19tct_T1BM)] /* [1,10000] */;
    uint32_t NN_P19tct_T1B /* [0,10000] */;
    uint32_t NN_P19tct_T1Bv[10000 EXT_STRUCT_CTRL(NN_P19tct_T1B)] /* [0,65535] */;
    uint32_t NN_P19tct_T2BM /* [1,50] */;
    uint32_t NN_P19tct_T2BMI[50 EXT_STRUCT_CTRL(NN_P19tct_T2BM)] /* [1,50] */;
    uint32_t NN_P19tct_T2BME[50 EXT_STRUCT_CTRL(NN_P19tct_T2BM)] /* [1,10000] */;
    uint32_t NN_P19tct_T2B /* [0,10000] */;
    uint32_t NN_P19tct_T2Bv[10000 EXT_STRUCT_CTRL(NN_P19tct_T2B)] /* [0,65535] */;
    uint32_t NN_P19tft_T1BM /* [1,50] */;
    uint32_t NN_P19tft_T1BMI[50 EXT_STRUCT_CTRL(NN_P19tft_T1BM)] /* [1,50] */;
    uint32_t NN_P19tft_T1BME[50 EXT_STRUCT_CTRL(NN_P19tft_T1BM)] /* [1,10000] */;
    uint32_t NN_P19tft_T1B /* [0,10000] */;
    uint32_t NN_P19tft_T1Bv[10000 EXT_STRUCT_CTRL(NN_P19tft_T1B)] /* [0,65535] */;
    uint32_t NN_P19tft_T2BM /* [1,50] */;
    uint32_t NN_P19tft_T2BMI[50 EXT_STRUCT_CTRL(NN_P19tft_T2BM)] /* [1,50] */;
    uint32_t NN_P19tft_T2BME[50 EXT_STRUCT_CTRL(NN_P19tft_T2BM)] /* [1,10000] */;
    uint32_t NN_P19tft_T2B /* [0,10000] */;
    uint32_t NN_P19tft_T2Bv[10000 EXT_STRUCT_CTRL(NN_P19tft_T2B)] /* [0,65535] */;
    uint32_t NN_P20tcl_T1BM /* [1,50] */;
    uint32_t NN_P20tcl_T1BMI[50 EXT_STRUCT_CTRL(NN_P20tcl_T1BM)] /* [1,50] */;
    uint32_t NN_P20tcl_T1BME[50 EXT_STRUCT_CTRL(NN_P20tcl_T1BM)] /* [1,10000] */;
    uint32_t NN_P20tcl_T1B /* [0,10000] */;
    uint32_t NN_P20tcl_T1Bv[10000 EXT_STRUCT_CTRL(NN_P20tcl_T1B)] /* [0,65535] */;
    uint32_t NN_P20tcl_T2BM /* [1,50] */;
    uint32_t NN_P20tcl_T2BMI[50 EXT_STRUCT_CTRL(NN_P20tcl_T2BM)] /* [1,50] */;
    uint32_t NN_P20tcl_T2BME[50 EXT_STRUCT_CTRL(NN_P20tcl_T2BM)] /* [1,10000] */;
    uint32_t NN_P20tcl_T2B /* [0,10000] */;
    uint32_t NN_P20tcl_T2Bv[10000 EXT_STRUCT_CTRL(NN_P20tcl_T2B)] /* [0,65535] */;
    uint32_t NN_P20tfl_T1BM /* [1,50] */;
    uint32_t NN_P20tfl_T1BMI[50 EXT_STRUCT_CTRL(NN_P20tfl_T1BM)] /* [1,50] */;
    uint32_t NN_P20tfl_T1BME[50 EXT_STRUCT_CTRL(NN_P20tfl_T1BM)] /* [1,10000] */;
    uint32_t NN_P20tfl_T1B /* [0,10000] */;
    uint32_t NN_P20tfl_T1Bv[10000 EXT_STRUCT_CTRL(NN_P20tfl_T1B)] /* [0,65535] */;
    uint32_t NN_P20tfl_T2BM /* [1,50] */;
    uint32_t NN_P20tfl_T2BMI[50 EXT_STRUCT_CTRL(NN_P20tfl_T2BM)] /* [1,50] */;
    uint32_t NN_P20tfl_T2BME[50 EXT_STRUCT_CTRL(NN_P20tfl_T2BM)] /* [1,10000] */;
    uint32_t NN_P20tfl_T2B /* [0,10000] */;
    uint32_t NN_P20tfl_T2Bv[10000 EXT_STRUCT_CTRL(NN_P20tfl_T2B)] /* [0,65535] */;
    uint32_t NN_P20tct_T1BM /* [1,50] */;
    uint32_t NN_P20tct_T1BMI[50 EXT_STRUCT_CTRL(NN_P20tct_T1BM)] /* [1,50] */;
    uint32_t NN_P20tct_T1BME[50 EXT_STRUCT_CTRL(NN_P20tct_T1BM)] /* [1,10000] */;
    uint32_t NN_P20tct_T1B /* [0,10000] */;
    uint32_t NN_P20tct_T1Bv[10000 EXT_STRUCT_CTRL(NN_P20tct_T1B)] /* [0,65535] */;
    uint32_t NN_P20tct_T2BM /* [1,50] */;
    uint32_t NN_P20tct_T2BMI[50 EXT_STRUCT_CTRL(NN_P20tct_T2BM)] /* [1,50] */;
    uint32_t NN_P20tct_T2BME[50 EXT_STRUCT_CTRL(NN_P20tct_T2BM)] /* [1,10000] */;
    uint32_t NN_P20tct_T2B /* [0,10000] */;
    uint32_t NN_P20tct_T2Bv[10000 EXT_STRUCT_CTRL(NN_P20tct_T2B)] /* [0,65535] */;
    uint32_t NN_P20tft_T1BM /* [1,50] */;
    uint32_t NN_P20tft_T1BMI[50 EXT_STRUCT_CTRL(NN_P20tft_T1BM)] /* [1,50] */;
    uint32_t NN_P20tft_T1BME[50 EXT_STRUCT_CTRL(NN_P20tft_T1BM)] /* [1,10000] */;
    uint32_t NN_P20tft_T1B /* [0,10000] */;
    uint32_t NN_P20tft_T1Bv[10000 EXT_STRUCT_CTRL(NN_P20tft_T1B)] /* [0,65535] */;
    uint32_t NN_P20tft_T2BM /* [1,50] */;
    uint32_t NN_P20tft_T2BMI[50 EXT_STRUCT_CTRL(NN_P20tft_T2BM)] /* [1,50] */;
    uint32_t NN_P20tft_T2BME[50 EXT_STRUCT_CTRL(NN_P20tft_T2BM)] /* [1,10000] */;
    uint32_t NN_P20tft_T2B /* [0,10000] */;
    uint32_t NN_P20tft_T2Bv[10000 EXT_STRUCT_CTRL(NN_P20tft_T2B)] /* [0,65535] */;
    uint32_t NN_P21tcl_T1BM /* [1,50] */;
    uint32_t NN_P21tcl_T1BMI[50 EXT_STRUCT_CTRL(NN_P21tcl_T1BM)] /* [1,50] */;
    uint32_t NN_P21tcl_T1BME[50 EXT_STRUCT_CTRL(NN_P21tcl_T1BM)] /* [1,10000] */;
    uint32_t NN_P21tcl_T1B /* [0,10000] */;
    uint32_t NN_P21tcl_T1Bv[10000 EXT_STRUCT_CTRL(NN_P21tcl_T1B)] /* [0,65535] */;
    uint32_t NN_P21tcl_T2BM /* [1,50] */;
    uint32_t NN_P21tcl_T2BMI[50 EXT_STRUCT_CTRL(NN_P21tcl_T2BM)] /* [1,50] */;
    uint32_t NN_P21tcl_T2BME[50 EXT_STRUCT_CTRL(NN_P21tcl_T2BM)] /* [1,10000] */;
    uint32_t NN_P21tcl_T2B /* [0,10000] */;
    uint32_t NN_P21tcl_T2Bv[10000 EXT_STRUCT_CTRL(NN_P21tcl_T2B)] /* [0,65535] */;
    uint32_t NN_P21tfl_T1BM /* [1,50] */;
    uint32_t NN_P21tfl_T1BMI[50 EXT_STRUCT_CTRL(NN_P21tfl_T1BM)] /* [1,50] */;
    uint32_t NN_P21tfl_T1BME[50 EXT_STRUCT_CTRL(NN_P21tfl_T1BM)] /* [1,10000] */;
    uint32_t NN_P21tfl_T1B /* [0,10000] */;
    uint32_t NN_P21tfl_T1Bv[10000 EXT_STRUCT_CTRL(NN_P21tfl_T1B)] /* [0,65535] */;
    uint32_t NN_P21tfl_T2BM /* [1,50] */;
    uint32_t NN_P21tfl_T2BMI[50 EXT_STRUCT_CTRL(NN_P21tfl_T2BM)] /* [1,50] */;
    uint32_t NN_P21tfl_T2BME[50 EXT_STRUCT_CTRL(NN_P21tfl_T2BM)] /* [1,10000] */;
    uint32_t NN_P21tfl_T2B /* [0,10000] */;
    uint32_t NN_P21tfl_T2Bv[10000 EXT_STRUCT_CTRL(NN_P21tfl_T2B)] /* [0,65535] */;
    uint32_t NN_P21tct_T1BM /* [1,50] */;
    uint32_t NN_P21tct_T1BMI[50 EXT_STRUCT_CTRL(NN_P21tct_T1BM)] /* [1,50] */;
    uint32_t NN_P21tct_T1BME[50 EXT_STRUCT_CTRL(NN_P21tct_T1BM)] /* [1,10000] */;
    uint32_t NN_P21tct_T1B /* [0,10000] */;
    uint32_t NN_P21tct_T1Bv[10000 EXT_STRUCT_CTRL(NN_P21tct_T1B)] /* [0,65535] */;
    uint32_t NN_P21tct_T2BM /* [1,50] */;
    uint32_t NN_P21tct_T2BMI[50 EXT_STRUCT_CTRL(NN_P21tct_T2BM)] /* [1,50] */;
    uint32_t NN_P21tct_T2BME[50 EXT_STRUCT_CTRL(NN_P21tct_T2BM)] /* [1,10000] */;
    uint32_t NN_P21tct_T2B /* [0,10000] */;
    uint32_t NN_P21tct_T2Bv[10000 EXT_STRUCT_CTRL(NN_P21tct_T2B)] /* [0,65535] */;
    uint32_t NN_P21tft_T1BM /* [1,50] */;
    uint32_t NN_P21tft_T1BMI[50 EXT_STRUCT_CTRL(NN_P21tft_T1BM)] /* [1,50] */;
    uint32_t NN_P21tft_T1BME[50 EXT_STRUCT_CTRL(NN_P21tft_T1BM)] /* [1,10000] */;
    uint32_t NN_P21tft_T1B /* [0,10000] */;
    uint32_t NN_P21tft_T1Bv[10000 EXT_STRUCT_CTRL(NN_P21tft_T1B)] /* [0,65535] */;
    uint32_t NN_P21tft_T2BM /* [1,50] */;
    uint32_t NN_P21tft_T2BMI[50 EXT_STRUCT_CTRL(NN_P21tft_T2BM)] /* [1,50] */;
    uint32_t NN_P21tft_T2BME[50 EXT_STRUCT_CTRL(NN_P21tft_T2BM)] /* [1,10000] */;
    uint32_t NN_P21tft_T2B /* [0,10000] */;
    uint32_t NN_P21tft_T2Bv[10000 EXT_STRUCT_CTRL(NN_P21tft_T2B)] /* [0,65535] */;
    uint32_t NN_P22tcl_T1BM /* [1,50] */;
    uint32_t NN_P22tcl_T1BMI[50 EXT_STRUCT_CTRL(NN_P22tcl_T1BM)] /* [1,50] */;
    uint32_t NN_P22tcl_T1BME[50 EXT_STRUCT_CTRL(NN_P22tcl_T1BM)] /* [1,10000] */;
    uint32_t NN_P22tcl_T1B /* [0,10000] */;
    uint32_t NN_P22tcl_T1Bv[10000 EXT_STRUCT_CTRL(NN_P22tcl_T1B)] /* [0,65535] */;
    uint32_t NN_P22tcl_T2BM /* [1,50] */;
    uint32_t NN_P22tcl_T2BMI[50 EXT_STRUCT_CTRL(NN_P22tcl_T2BM)] /* [1,50] */;
    uint32_t NN_P22tcl_T2BME[50 EXT_STRUCT_CTRL(NN_P22tcl_T2BM)] /* [1,10000] */;
    uint32_t NN_P22tcl_T2B /* [0,10000] */;
    uint32_t NN_P22tcl_T2Bv[10000 EXT_STRUCT_CTRL(NN_P22tcl_T2B)] /* [0,65535] */;
    uint32_t NN_P22tfl_T1BM /* [1,50] */;
    uint32_t NN_P22tfl_T1BMI[50 EXT_STRUCT_CTRL(NN_P22tfl_T1BM)] /* [1,50] */;
    uint32_t NN_P22tfl_T1BME[50 EXT_STRUCT_CTRL(NN_P22tfl_T1BM)] /* [1,10000] */;
    uint32_t NN_P22tfl_T1B /* [0,10000] */;
    uint32_t NN_P22tfl_T1Bv[10000 EXT_STRUCT_CTRL(NN_P22tfl_T1B)] /* [0,65535] */;
    uint32_t NN_P22tfl_T2BM /* [1,50] */;
    uint32_t NN_P22tfl_T2BMI[50 EXT_STRUCT_CTRL(NN_P22tfl_T2BM)] /* [1,50] */;
    uint32_t NN_P22tfl_T2BME[50 EXT_STRUCT_CTRL(NN_P22tfl_T2BM)] /* [1,10000] */;
    uint32_t NN_P22tfl_T2B /* [0,10000] */;
    uint32_t NN_P22tfl_T2Bv[10000 EXT_STRUCT_CTRL(NN_P22tfl_T2B)] /* [0,65535] */;
    uint32_t NN_P22tct_T1BM /* [1,50] */;
    uint32_t NN_P22tct_T1BMI[50 EXT_STRUCT_CTRL(NN_P22tct_T1BM)] /* [1,50] */;
    uint32_t NN_P22tct_T1BME[50 EXT_STRUCT_CTRL(NN_P22tct_T1BM)] /* [1,10000] */;
    uint32_t NN_P22tct_T1B /* [0,10000] */;
    uint32_t NN_P22tct_T1Bv[10000 EXT_STRUCT_CTRL(NN_P22tct_T1B)] /* [0,65535] */;
    uint32_t NN_P22tct_T2BM /* [1,50] */;
    uint32_t NN_P22tct_T2BMI[50 EXT_STRUCT_CTRL(NN_P22tct_T2BM)] /* [1,50] */;
    uint32_t NN_P22tct_T2BME[50 EXT_STRUCT_CTRL(NN_P22tct_T2BM)] /* [1,10000] */;
    uint32_t NN_P22tct_T2B /* [0,10000] */;
    uint32_t NN_P22tct_T2Bv[10000 EXT_STRUCT_CTRL(NN_P22tct_T2B)] /* [0,65535] */;
    uint32_t NN_P22tft_T1BM /* [1,50] */;
    uint32_t NN_P22tft_T1BMI[50 EXT_STRUCT_CTRL(NN_P22tft_T1BM)] /* [1,50] */;
    uint32_t NN_P22tft_T1BME[50 EXT_STRUCT_CTRL(NN_P22tft_T1BM)] /* [1,10000] */;
    uint32_t NN_P22tft_T1B /* [0,10000] */;
    uint32_t NN_P22tft_T1Bv[10000 EXT_STRUCT_CTRL(NN_P22tft_T1B)] /* [0,65535] */;
    uint32_t NN_P22tft_T2BM /* [1,50] */;
    uint32_t NN_P22tft_T2BMI[50 EXT_STRUCT_CTRL(NN_P22tft_T2BM)] /* [1,50] */;
    uint32_t NN_P22tft_T2BME[50 EXT_STRUCT_CTRL(NN_P22tft_T2BM)] /* [1,10000] */;
    uint32_t NN_P22tft_T2B /* [0,10000] */;
    uint32_t NN_P22tft_T2Bv[10000 EXT_STRUCT_CTRL(NN_P22tft_T2B)] /* [0,65535] */;
    uint32_t NN_P23tcl_T1BM /* [1,50] */;
    uint32_t NN_P23tcl_T1BMI[50 EXT_STRUCT_CTRL(NN_P23tcl_T1BM)] /* [1,50] */;
    uint32_t NN_P23tcl_T1BME[50 EXT_STRUCT_CTRL(NN_P23tcl_T1BM)] /* [1,10000] */;
    uint32_t NN_P23tcl_T1B /* [0,10000] */;
    uint32_t NN_P23tcl_T1Bv[10000 EXT_STRUCT_CTRL(NN_P23tcl_T1B)] /* [0,65535] */;
    uint32_t NN_P23tcl_T2BM /* [1,50] */;
    uint32_t NN_P23tcl_T2BMI[50 EXT_STRUCT_CTRL(NN_P23tcl_T2BM)] /* [1,50] */;
    uint32_t NN_P23tcl_T2BME[50 EXT_STRUCT_CTRL(NN_P23tcl_T2BM)] /* [1,10000] */;
    uint32_t NN_P23tcl_T2B /* [0,10000] */;
    uint32_t NN_P23tcl_T2Bv[10000 EXT_STRUCT_CTRL(NN_P23tcl_T2B)] /* [0,65535] */;
    uint32_t NN_P23tfl_T1BM /* [1,50] */;
    uint32_t NN_P23tfl_T1BMI[50 EXT_STRUCT_CTRL(NN_P23tfl_T1BM)] /* [1,50] */;
    uint32_t NN_P23tfl_T1BME[50 EXT_STRUCT_CTRL(NN_P23tfl_T1BM)] /* [1,10000] */;
    uint32_t NN_P23tfl_T1B /* [0,10000] */;
    uint32_t NN_P23tfl_T1Bv[10000 EXT_STRUCT_CTRL(NN_P23tfl_T1B)] /* [0,65535] */;
    uint32_t NN_P23tfl_T2BM /* [1,50] */;
    uint32_t NN_P23tfl_T2BMI[50 EXT_STRUCT_CTRL(NN_P23tfl_T2BM)] /* [1,50] */;
    uint32_t NN_P23tfl_T2BME[50 EXT_STRUCT_CTRL(NN_P23tfl_T2BM)] /* [1,10000] */;
    uint32_t NN_P23tfl_T2B /* [0,10000] */;
    uint32_t NN_P23tfl_T2Bv[10000 EXT_STRUCT_CTRL(NN_P23tfl_T2B)] /* [0,65535] */;
    uint32_t NN_P23tct_T1BM /* [1,50] */;
    uint32_t NN_P23tct_T1BMI[50 EXT_STRUCT_CTRL(NN_P23tct_T1BM)] /* [1,50] */;
    uint32_t NN_P23tct_T1BME[50 EXT_STRUCT_CTRL(NN_P23tct_T1BM)] /* [1,10000] */;
    uint32_t NN_P23tct_T1B /* [0,10000] */;
    uint32_t NN_P23tct_T1Bv[10000 EXT_STRUCT_CTRL(NN_P23tct_T1B)] /* [0,65535] */;
    uint32_t NN_P23tct_T2BM /* [1,50] */;
    uint32_t NN_P23tct_T2BMI[50 EXT_STRUCT_CTRL(NN_P23tct_T2BM)] /* [1,50] */;
    uint32_t NN_P23tct_T2BME[50 EXT_STRUCT_CTRL(NN_P23tct_T2BM)] /* [1,10000] */;
    uint32_t NN_P23tct_T2B /* [0,10000] */;
    uint32_t NN_P23tct_T2Bv[10000 EXT_STRUCT_CTRL(NN_P23tct_T2B)] /* [0,65535] */;
    uint32_t NN_P23tft_T1BM /* [1,50] */;
    uint32_t NN_P23tft_T1BMI[50 EXT_STRUCT_CTRL(NN_P23tft_T1BM)] /* [1,50] */;
    uint32_t NN_P23tft_T1BME[50 EXT_STRUCT_CTRL(NN_P23tft_T1BM)] /* [1,10000] */;
    uint32_t NN_P23tft_T1B /* [0,10000] */;
    uint32_t NN_P23tft_T1Bv[10000 EXT_STRUCT_CTRL(NN_P23tft_T1B)] /* [0,65535] */;
    uint32_t NN_P23tft_T2BM /* [1,50] */;
    uint32_t NN_P23tft_T2BMI[50 EXT_STRUCT_CTRL(NN_P23tft_T2BM)] /* [1,50] */;
    uint32_t NN_P23tft_T2BME[50 EXT_STRUCT_CTRL(NN_P23tft_T2BM)] /* [1,10000] */;
    uint32_t NN_P23tft_T2B /* [0,10000] */;
    uint32_t NN_P23tft_T2Bv[10000 EXT_STRUCT_CTRL(NN_P23tft_T2B)] /* [0,65535] */;
    uint32_t NN_P24tcl_T1BM /* [1,50] */;
    uint32_t NN_P24tcl_T1BMI[50 EXT_STRUCT_CTRL(NN_P24tcl_T1BM)] /* [1,50] */;
    uint32_t NN_P24tcl_T1BME[50 EXT_STRUCT_CTRL(NN_P24tcl_T1BM)] /* [1,10000] */;
    uint32_t NN_P24tcl_T1B /* [0,10000] */;
    uint32_t NN_P24tcl_T1Bv[10000 EXT_STRUCT_CTRL(NN_P24tcl_T1B)] /* [0,65535] */;
    uint32_t NN_P24tcl_T2BM /* [1,50] */;
    uint32_t NN_P24tcl_T2BMI[50 EXT_STRUCT_CTRL(NN_P24tcl_T2BM)] /* [1,50] */;
    uint32_t NN_P24tcl_T2BME[50 EXT_STRUCT_CTRL(NN_P24tcl_T2BM)] /* [1,10000] */;
    uint32_t NN_P24tcl_T2B /* [0,10000] */;
    uint32_t NN_P24tcl_T2Bv[10000 EXT_STRUCT_CTRL(NN_P24tcl_T2B)] /* [0,65535] */;
    uint32_t NN_P24tfl_T1BM /* [1,50] */;
    uint32_t NN_P24tfl_T1BMI[50 EXT_STRUCT_CTRL(NN_P24tfl_T1BM)] /* [1,50] */;
    uint32_t NN_P24tfl_T1BME[50 EXT_STRUCT_CTRL(NN_P24tfl_T1BM)] /* [1,10000] */;
    uint32_t NN_P24tfl_T1B /* [0,10000] */;
    uint32_t NN_P24tfl_T1Bv[10000 EXT_STRUCT_CTRL(NN_P24tfl_T1B)] /* [0,65535] */;
    uint32_t NN_P24tfl_T2BM /* [1,50] */;
    uint32_t NN_P24tfl_T2BMI[50 EXT_STRUCT_CTRL(NN_P24tfl_T2BM)] /* [1,50] */;
    uint32_t NN_P24tfl_T2BME[50 EXT_STRUCT_CTRL(NN_P24tfl_T2BM)] /* [1,10000] */;
    uint32_t NN_P24tfl_T2B /* [0,10000] */;
    uint32_t NN_P24tfl_T2Bv[10000 EXT_STRUCT_CTRL(NN_P24tfl_T2B)] /* [0,65535] */;
    uint32_t NN_P24tct_T1BM /* [1,50] */;
    uint32_t NN_P24tct_T1BMI[50 EXT_STRUCT_CTRL(NN_P24tct_T1BM)] /* [1,50] */;
    uint32_t NN_P24tct_T1BME[50 EXT_STRUCT_CTRL(NN_P24tct_T1BM)] /* [1,10000] */;
    uint32_t NN_P24tct_T1B /* [0,10000] */;
    uint32_t NN_P24tct_T1Bv[10000 EXT_STRUCT_CTRL(NN_P24tct_T1B)] /* [0,65535] */;
    uint32_t NN_P24tct_T2BM /* [1,50] */;
    uint32_t NN_P24tct_T2BMI[50 EXT_STRUCT_CTRL(NN_P24tct_T2BM)] /* [1,50] */;
    uint32_t NN_P24tct_T2BME[50 EXT_STRUCT_CTRL(NN_P24tct_T2BM)] /* [1,10000] */;
    uint32_t NN_P24tct_T2B /* [0,10000] */;
    uint32_t NN_P24tct_T2Bv[10000 EXT_STRUCT_CTRL(NN_P24tct_T2B)] /* [0,65535] */;
    uint32_t NN_P24tft_T1BM /* [1,50] */;
    uint32_t NN_P24tft_T1BMI[50 EXT_STRUCT_CTRL(NN_P24tft_T1BM)] /* [1,50] */;
    uint32_t NN_P24tft_T1BME[50 EXT_STRUCT_CTRL(NN_P24tft_T1BM)] /* [1,10000] */;
    uint32_t NN_P24tft_T1B /* [0,10000] */;
    uint32_t NN_P24tft_T1Bv[10000 EXT_STRUCT_CTRL(NN_P24tft_T1B)] /* [0,65535] */;
    uint32_t NN_P24tft_T2BM /* [1,50] */;
    uint32_t NN_P24tft_T2BMI[50 EXT_STRUCT_CTRL(NN_P24tft_T2BM)] /* [1,50] */;
    uint32_t NN_P24tft_T2BME[50 EXT_STRUCT_CTRL(NN_P24tft_T2BM)] /* [1,10000] */;
    uint32_t NN_P24tft_T2B /* [0,10000] */;
    uint32_t NN_P24tft_T2Bv[10000 EXT_STRUCT_CTRL(NN_P24tft_T2B)] /* [0,65535] */;

} EXT_STR_h101_raw_nnp_tamex;

/********************************************************
 *
 * Structure with multiple levels of arrays (partially)
 * recovered (recommended):
 */

typedef struct EXT_STR_h101_raw_nnp_tamex_onion_t
{
    /* RAW */
    struct
    {
        struct
        {
            uint32_t BM;
            uint32_t BMI[50 /* BM */];
            uint32_t BME[50 /* BM */];
            uint32_t B;
            uint32_t Bv[10000 /* B */];
        } tcl_T[2];
        struct
        {
            uint32_t BM;
            uint32_t BMI[50 /* BM */];
            uint32_t BME[50 /* BM */];
            uint32_t B;
            uint32_t Bv[10000 /* B */];
        } tfl_T[2];
        struct
        {
            uint32_t BM;
            uint32_t BMI[50 /* BM */];
            uint32_t BME[50 /* BM */];
            uint32_t B;
            uint32_t Bv[10000 /* B */];
        } tct_T[2];
        struct
        {
            uint32_t BM;
            uint32_t BMI[50 /* BM */];
            uint32_t BME[50 /* BM */];
            uint32_t B;
            uint32_t Bv[10000 /* B */];
        } tft_T[2];
    } NN_P[24];

} EXT_STR_h101_raw_nnp_tamex_onion;

/*******************************************************/

#define EXT_STR_h101_raw_nnp_tamex_ITEMS_INFO(ok, si, offset, struct_t, printerr)                                 \
    do                                                                                                            \
    {                                                                                                             \
        ok = 1;                                                                                                   \
        /* RAW */                                                                                                 \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, NN_P1tcl_T1BM, UINT32, "NN_P1tcl_T1BM", 50);    \
        EXT_STR_ITEM_INFO_ZZP(                                                                                    \
            ok, si, offset, struct_t, printerr, NN_P1tcl_T1BMI, UINT32, "NN_P1tcl_T1BMI", "NN_P1tcl_T1BM");       \
        EXT_STR_ITEM_INFO_ZZP(                                                                                    \
            ok, si, offset, struct_t, printerr, NN_P1tcl_T1BME, UINT32, "NN_P1tcl_T1BME", "NN_P1tcl_T1BM");       \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, NN_P1tcl_T1B, UINT32, "NN_P1tcl_T1B", 10000);   \
        EXT_STR_ITEM_INFO_ZZP(                                                                                    \
            ok, si, offset, struct_t, printerr, NN_P1tcl_T1Bv, UINT32, "NN_P1tcl_T1Bv", "NN_P1tcl_T1B");          \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, NN_P1tcl_T2BM, UINT32, "NN_P1tcl_T2BM", 50);    \
        EXT_STR_ITEM_INFO_ZZP(                                                                                    \
            ok, si, offset, struct_t, printerr, NN_P1tcl_T2BMI, UINT32, "NN_P1tcl_T2BMI", "NN_P1tcl_T2BM");       \
        EXT_STR_ITEM_INFO_ZZP(                                                                                    \
            ok, si, offset, struct_t, printerr, NN_P1tcl_T2BME, UINT32, "NN_P1tcl_T2BME", "NN_P1tcl_T2BM");       \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, NN_P1tcl_T2B, UINT32, "NN_P1tcl_T2B", 10000);   \
        EXT_STR_ITEM_INFO_ZZP(                                                                                    \
            ok, si, offset, struct_t, printerr, NN_P1tcl_T2Bv, UINT32, "NN_P1tcl_T2Bv", "NN_P1tcl_T2B");          \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, NN_P1tfl_T1BM, UINT32, "NN_P1tfl_T1BM", 50);    \
        EXT_STR_ITEM_INFO_ZZP(                                                                                    \
            ok, si, offset, struct_t, printerr, NN_P1tfl_T1BMI, UINT32, "NN_P1tfl_T1BMI", "NN_P1tfl_T1BM");       \
        EXT_STR_ITEM_INFO_ZZP(                                                                                    \
            ok, si, offset, struct_t, printerr, NN_P1tfl_T1BME, UINT32, "NN_P1tfl_T1BME", "NN_P1tfl_T1BM");       \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, NN_P1tfl_T1B, UINT32, "NN_P1tfl_T1B", 10000);   \
        EXT_STR_ITEM_INFO_ZZP(                                                                                    \
            ok, si, offset, struct_t, printerr, NN_P1tfl_T1Bv, UINT32, "NN_P1tfl_T1Bv", "NN_P1tfl_T1B");          \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, NN_P1tfl_T2BM, UINT32, "NN_P1tfl_T2BM", 50);    \
        EXT_STR_ITEM_INFO_ZZP(                                                                                    \
            ok, si, offset, struct_t, printerr, NN_P1tfl_T2BMI, UINT32, "NN_P1tfl_T2BMI", "NN_P1tfl_T2BM");       \
        EXT_STR_ITEM_INFO_ZZP(                                                                                    \
            ok, si, offset, struct_t, printerr, NN_P1tfl_T2BME, UINT32, "NN_P1tfl_T2BME", "NN_P1tfl_T2BM");       \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, NN_P1tfl_T2B, UINT32, "NN_P1tfl_T2B", 10000);   \
        EXT_STR_ITEM_INFO_ZZP(                                                                                    \
            ok, si, offset, struct_t, printerr, NN_P1tfl_T2Bv, UINT32, "NN_P1tfl_T2Bv", "NN_P1tfl_T2B");          \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, NN_P1tct_T1BM, UINT32, "NN_P1tct_T1BM", 50);    \
        EXT_STR_ITEM_INFO_ZZP(                                                                                    \
            ok, si, offset, struct_t, printerr, NN_P1tct_T1BMI, UINT32, "NN_P1tct_T1BMI", "NN_P1tct_T1BM");       \
        EXT_STR_ITEM_INFO_ZZP(                                                                                    \
            ok, si, offset, struct_t, printerr, NN_P1tct_T1BME, UINT32, "NN_P1tct_T1BME", "NN_P1tct_T1BM");       \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, NN_P1tct_T1B, UINT32, "NN_P1tct_T1B", 10000);   \
        EXT_STR_ITEM_INFO_ZZP(                                                                                    \
            ok, si, offset, struct_t, printerr, NN_P1tct_T1Bv, UINT32, "NN_P1tct_T1Bv", "NN_P1tct_T1B");          \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, NN_P1tct_T2BM, UINT32, "NN_P1tct_T2BM", 50);    \
        EXT_STR_ITEM_INFO_ZZP(                                                                                    \
            ok, si, offset, struct_t, printerr, NN_P1tct_T2BMI, UINT32, "NN_P1tct_T2BMI", "NN_P1tct_T2BM");       \
        EXT_STR_ITEM_INFO_ZZP(                                                                                    \
            ok, si, offset, struct_t, printerr, NN_P1tct_T2BME, UINT32, "NN_P1tct_T2BME", "NN_P1tct_T2BM");       \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, NN_P1tct_T2B, UINT32, "NN_P1tct_T2B", 10000);   \
        EXT_STR_ITEM_INFO_ZZP(                                                                                    \
            ok, si, offset, struct_t, printerr, NN_P1tct_T2Bv, UINT32, "NN_P1tct_T2Bv", "NN_P1tct_T2B");          \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, NN_P1tft_T1BM, UINT32, "NN_P1tft_T1BM", 50);    \
        EXT_STR_ITEM_INFO_ZZP(                                                                                    \
            ok, si, offset, struct_t, printerr, NN_P1tft_T1BMI, UINT32, "NN_P1tft_T1BMI", "NN_P1tft_T1BM");       \
        EXT_STR_ITEM_INFO_ZZP(                                                                                    \
            ok, si, offset, struct_t, printerr, NN_P1tft_T1BME, UINT32, "NN_P1tft_T1BME", "NN_P1tft_T1BM");       \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, NN_P1tft_T1B, UINT32, "NN_P1tft_T1B", 10000);   \
        EXT_STR_ITEM_INFO_ZZP(                                                                                    \
            ok, si, offset, struct_t, printerr, NN_P1tft_T1Bv, UINT32, "NN_P1tft_T1Bv", "NN_P1tft_T1B");          \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, NN_P1tft_T2BM, UINT32, "NN_P1tft_T2BM", 50);    \
        EXT_STR_ITEM_INFO_ZZP(                                                                                    \
            ok, si, offset, struct_t, printerr, NN_P1tft_T2BMI, UINT32, "NN_P1tft_T2BMI", "NN_P1tft_T2BM");       \
        EXT_STR_ITEM_INFO_ZZP(                                                                                    \
            ok, si, offset, struct_t, printerr, NN_P1tft_T2BME, UINT32, "NN_P1tft_T2BME", "NN_P1tft_T2BM");       \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, NN_P1tft_T2B, UINT32, "NN_P1tft_T2B", 10000);   \
        EXT_STR_ITEM_INFO_ZZP(                                                                                    \
            ok, si, offset, struct_t, printerr, NN_P1tft_T2Bv, UINT32, "NN_P1tft_T2Bv", "NN_P1tft_T2B");          \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, NN_P2tcl_T1BM, UINT32, "NN_P2tcl_T1BM", 50);    \
        EXT_STR_ITEM_INFO_ZZP(                                                                                    \
            ok, si, offset, struct_t, printerr, NN_P2tcl_T1BMI, UINT32, "NN_P2tcl_T1BMI", "NN_P2tcl_T1BM");       \
        EXT_STR_ITEM_INFO_ZZP(                                                                                    \
            ok, si, offset, struct_t, printerr, NN_P2tcl_T1BME, UINT32, "NN_P2tcl_T1BME", "NN_P2tcl_T1BM");       \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, NN_P2tcl_T1B, UINT32, "NN_P2tcl_T1B", 10000);   \
        EXT_STR_ITEM_INFO_ZZP(                                                                                    \
            ok, si, offset, struct_t, printerr, NN_P2tcl_T1Bv, UINT32, "NN_P2tcl_T1Bv", "NN_P2tcl_T1B");          \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, NN_P2tcl_T2BM, UINT32, "NN_P2tcl_T2BM", 50);    \
        EXT_STR_ITEM_INFO_ZZP(                                                                                    \
            ok, si, offset, struct_t, printerr, NN_P2tcl_T2BMI, UINT32, "NN_P2tcl_T2BMI", "NN_P2tcl_T2BM");       \
        EXT_STR_ITEM_INFO_ZZP(                                                                                    \
            ok, si, offset, struct_t, printerr, NN_P2tcl_T2BME, UINT32, "NN_P2tcl_T2BME", "NN_P2tcl_T2BM");       \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, NN_P2tcl_T2B, UINT32, "NN_P2tcl_T2B", 10000);   \
        EXT_STR_ITEM_INFO_ZZP(                                                                                    \
            ok, si, offset, struct_t, printerr, NN_P2tcl_T2Bv, UINT32, "NN_P2tcl_T2Bv", "NN_P2tcl_T2B");          \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, NN_P2tfl_T1BM, UINT32, "NN_P2tfl_T1BM", 50);    \
        EXT_STR_ITEM_INFO_ZZP(                                                                                    \
            ok, si, offset, struct_t, printerr, NN_P2tfl_T1BMI, UINT32, "NN_P2tfl_T1BMI", "NN_P2tfl_T1BM");       \
        EXT_STR_ITEM_INFO_ZZP(                                                                                    \
            ok, si, offset, struct_t, printerr, NN_P2tfl_T1BME, UINT32, "NN_P2tfl_T1BME", "NN_P2tfl_T1BM");       \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, NN_P2tfl_T1B, UINT32, "NN_P2tfl_T1B", 10000);   \
        EXT_STR_ITEM_INFO_ZZP(                                                                                    \
            ok, si, offset, struct_t, printerr, NN_P2tfl_T1Bv, UINT32, "NN_P2tfl_T1Bv", "NN_P2tfl_T1B");          \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, NN_P2tfl_T2BM, UINT32, "NN_P2tfl_T2BM", 50);    \
        EXT_STR_ITEM_INFO_ZZP(                                                                                    \
            ok, si, offset, struct_t, printerr, NN_P2tfl_T2BMI, UINT32, "NN_P2tfl_T2BMI", "NN_P2tfl_T2BM");       \
        EXT_STR_ITEM_INFO_ZZP(                                                                                    \
            ok, si, offset, struct_t, printerr, NN_P2tfl_T2BME, UINT32, "NN_P2tfl_T2BME", "NN_P2tfl_T2BM");       \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, NN_P2tfl_T2B, UINT32, "NN_P2tfl_T2B", 10000);   \
        EXT_STR_ITEM_INFO_ZZP(                                                                                    \
            ok, si, offset, struct_t, printerr, NN_P2tfl_T2Bv, UINT32, "NN_P2tfl_T2Bv", "NN_P2tfl_T2B");          \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, NN_P2tct_T1BM, UINT32, "NN_P2tct_T1BM", 50);    \
        EXT_STR_ITEM_INFO_ZZP(                                                                                    \
            ok, si, offset, struct_t, printerr, NN_P2tct_T1BMI, UINT32, "NN_P2tct_T1BMI", "NN_P2tct_T1BM");       \
        EXT_STR_ITEM_INFO_ZZP(                                                                                    \
            ok, si, offset, struct_t, printerr, NN_P2tct_T1BME, UINT32, "NN_P2tct_T1BME", "NN_P2tct_T1BM");       \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, NN_P2tct_T1B, UINT32, "NN_P2tct_T1B", 10000);   \
        EXT_STR_ITEM_INFO_ZZP(                                                                                    \
            ok, si, offset, struct_t, printerr, NN_P2tct_T1Bv, UINT32, "NN_P2tct_T1Bv", "NN_P2tct_T1B");          \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, NN_P2tct_T2BM, UINT32, "NN_P2tct_T2BM", 50);    \
        EXT_STR_ITEM_INFO_ZZP(                                                                                    \
            ok, si, offset, struct_t, printerr, NN_P2tct_T2BMI, UINT32, "NN_P2tct_T2BMI", "NN_P2tct_T2BM");       \
        EXT_STR_ITEM_INFO_ZZP(                                                                                    \
            ok, si, offset, struct_t, printerr, NN_P2tct_T2BME, UINT32, "NN_P2tct_T2BME", "NN_P2tct_T2BM");       \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, NN_P2tct_T2B, UINT32, "NN_P2tct_T2B", 10000);   \
        EXT_STR_ITEM_INFO_ZZP(                                                                                    \
            ok, si, offset, struct_t, printerr, NN_P2tct_T2Bv, UINT32, "NN_P2tct_T2Bv", "NN_P2tct_T2B");          \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, NN_P2tft_T1BM, UINT32, "NN_P2tft_T1BM", 50);    \
        EXT_STR_ITEM_INFO_ZZP(                                                                                    \
            ok, si, offset, struct_t, printerr, NN_P2tft_T1BMI, UINT32, "NN_P2tft_T1BMI", "NN_P2tft_T1BM");       \
        EXT_STR_ITEM_INFO_ZZP(                                                                                    \
            ok, si, offset, struct_t, printerr, NN_P2tft_T1BME, UINT32, "NN_P2tft_T1BME", "NN_P2tft_T1BM");       \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, NN_P2tft_T1B, UINT32, "NN_P2tft_T1B", 10000);   \
        EXT_STR_ITEM_INFO_ZZP(                                                                                    \
            ok, si, offset, struct_t, printerr, NN_P2tft_T1Bv, UINT32, "NN_P2tft_T1Bv", "NN_P2tft_T1B");          \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, NN_P2tft_T2BM, UINT32, "NN_P2tft_T2BM", 50);    \
        EXT_STR_ITEM_INFO_ZZP(                                                                                    \
            ok, si, offset, struct_t, printerr, NN_P2tft_T2BMI, UINT32, "NN_P2tft_T2BMI", "NN_P2tft_T2BM");       \
        EXT_STR_ITEM_INFO_ZZP(                                                                                    \
            ok, si, offset, struct_t, printerr, NN_P2tft_T2BME, UINT32, "NN_P2tft_T2BME", "NN_P2tft_T2BM");       \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, NN_P2tft_T2B, UINT32, "NN_P2tft_T2B", 10000);   \
        EXT_STR_ITEM_INFO_ZZP(                                                                                    \
            ok, si, offset, struct_t, printerr, NN_P2tft_T2Bv, UINT32, "NN_P2tft_T2Bv", "NN_P2tft_T2B");          \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, NN_P3tcl_T1BM, UINT32, "NN_P3tcl_T1BM", 50);    \
        EXT_STR_ITEM_INFO_ZZP(                                                                                    \
            ok, si, offset, struct_t, printerr, NN_P3tcl_T1BMI, UINT32, "NN_P3tcl_T1BMI", "NN_P3tcl_T1BM");       \
        EXT_STR_ITEM_INFO_ZZP(                                                                                    \
            ok, si, offset, struct_t, printerr, NN_P3tcl_T1BME, UINT32, "NN_P3tcl_T1BME", "NN_P3tcl_T1BM");       \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, NN_P3tcl_T1B, UINT32, "NN_P3tcl_T1B", 10000);   \
        EXT_STR_ITEM_INFO_ZZP(                                                                                    \
            ok, si, offset, struct_t, printerr, NN_P3tcl_T1Bv, UINT32, "NN_P3tcl_T1Bv", "NN_P3tcl_T1B");          \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, NN_P3tcl_T2BM, UINT32, "NN_P3tcl_T2BM", 50);    \
        EXT_STR_ITEM_INFO_ZZP(                                                                                    \
            ok, si, offset, struct_t, printerr, NN_P3tcl_T2BMI, UINT32, "NN_P3tcl_T2BMI", "NN_P3tcl_T2BM");       \
        EXT_STR_ITEM_INFO_ZZP(                                                                                    \
            ok, si, offset, struct_t, printerr, NN_P3tcl_T2BME, UINT32, "NN_P3tcl_T2BME", "NN_P3tcl_T2BM");       \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, NN_P3tcl_T2B, UINT32, "NN_P3tcl_T2B", 10000);   \
        EXT_STR_ITEM_INFO_ZZP(                                                                                    \
            ok, si, offset, struct_t, printerr, NN_P3tcl_T2Bv, UINT32, "NN_P3tcl_T2Bv", "NN_P3tcl_T2B");          \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, NN_P3tfl_T1BM, UINT32, "NN_P3tfl_T1BM", 50);    \
        EXT_STR_ITEM_INFO_ZZP(                                                                                    \
            ok, si, offset, struct_t, printerr, NN_P3tfl_T1BMI, UINT32, "NN_P3tfl_T1BMI", "NN_P3tfl_T1BM");       \
        EXT_STR_ITEM_INFO_ZZP(                                                                                    \
            ok, si, offset, struct_t, printerr, NN_P3tfl_T1BME, UINT32, "NN_P3tfl_T1BME", "NN_P3tfl_T1BM");       \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, NN_P3tfl_T1B, UINT32, "NN_P3tfl_T1B", 10000);   \
        EXT_STR_ITEM_INFO_ZZP(                                                                                    \
            ok, si, offset, struct_t, printerr, NN_P3tfl_T1Bv, UINT32, "NN_P3tfl_T1Bv", "NN_P3tfl_T1B");          \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, NN_P3tfl_T2BM, UINT32, "NN_P3tfl_T2BM", 50);    \
        EXT_STR_ITEM_INFO_ZZP(                                                                                    \
            ok, si, offset, struct_t, printerr, NN_P3tfl_T2BMI, UINT32, "NN_P3tfl_T2BMI", "NN_P3tfl_T2BM");       \
        EXT_STR_ITEM_INFO_ZZP(                                                                                    \
            ok, si, offset, struct_t, printerr, NN_P3tfl_T2BME, UINT32, "NN_P3tfl_T2BME", "NN_P3tfl_T2BM");       \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, NN_P3tfl_T2B, UINT32, "NN_P3tfl_T2B", 10000);   \
        EXT_STR_ITEM_INFO_ZZP(                                                                                    \
            ok, si, offset, struct_t, printerr, NN_P3tfl_T2Bv, UINT32, "NN_P3tfl_T2Bv", "NN_P3tfl_T2B");          \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, NN_P3tct_T1BM, UINT32, "NN_P3tct_T1BM", 50);    \
        EXT_STR_ITEM_INFO_ZZP(                                                                                    \
            ok, si, offset, struct_t, printerr, NN_P3tct_T1BMI, UINT32, "NN_P3tct_T1BMI", "NN_P3tct_T1BM");       \
        EXT_STR_ITEM_INFO_ZZP(                                                                                    \
            ok, si, offset, struct_t, printerr, NN_P3tct_T1BME, UINT32, "NN_P3tct_T1BME", "NN_P3tct_T1BM");       \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, NN_P3tct_T1B, UINT32, "NN_P3tct_T1B", 10000);   \
        EXT_STR_ITEM_INFO_ZZP(                                                                                    \
            ok, si, offset, struct_t, printerr, NN_P3tct_T1Bv, UINT32, "NN_P3tct_T1Bv", "NN_P3tct_T1B");          \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, NN_P3tct_T2BM, UINT32, "NN_P3tct_T2BM", 50);    \
        EXT_STR_ITEM_INFO_ZZP(                                                                                    \
            ok, si, offset, struct_t, printerr, NN_P3tct_T2BMI, UINT32, "NN_P3tct_T2BMI", "NN_P3tct_T2BM");       \
        EXT_STR_ITEM_INFO_ZZP(                                                                                    \
            ok, si, offset, struct_t, printerr, NN_P3tct_T2BME, UINT32, "NN_P3tct_T2BME", "NN_P3tct_T2BM");       \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, NN_P3tct_T2B, UINT32, "NN_P3tct_T2B", 10000);   \
        EXT_STR_ITEM_INFO_ZZP(                                                                                    \
            ok, si, offset, struct_t, printerr, NN_P3tct_T2Bv, UINT32, "NN_P3tct_T2Bv", "NN_P3tct_T2B");          \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, NN_P3tft_T1BM, UINT32, "NN_P3tft_T1BM", 50);    \
        EXT_STR_ITEM_INFO_ZZP(                                                                                    \
            ok, si, offset, struct_t, printerr, NN_P3tft_T1BMI, UINT32, "NN_P3tft_T1BMI", "NN_P3tft_T1BM");       \
        EXT_STR_ITEM_INFO_ZZP(                                                                                    \
            ok, si, offset, struct_t, printerr, NN_P3tft_T1BME, UINT32, "NN_P3tft_T1BME", "NN_P3tft_T1BM");       \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, NN_P3tft_T1B, UINT32, "NN_P3tft_T1B", 10000);   \
        EXT_STR_ITEM_INFO_ZZP(                                                                                    \
            ok, si, offset, struct_t, printerr, NN_P3tft_T1Bv, UINT32, "NN_P3tft_T1Bv", "NN_P3tft_T1B");          \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, NN_P3tft_T2BM, UINT32, "NN_P3tft_T2BM", 50);    \
        EXT_STR_ITEM_INFO_ZZP(                                                                                    \
            ok, si, offset, struct_t, printerr, NN_P3tft_T2BMI, UINT32, "NN_P3tft_T2BMI", "NN_P3tft_T2BM");       \
        EXT_STR_ITEM_INFO_ZZP(                                                                                    \
            ok, si, offset, struct_t, printerr, NN_P3tft_T2BME, UINT32, "NN_P3tft_T2BME", "NN_P3tft_T2BM");       \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, NN_P3tft_T2B, UINT32, "NN_P3tft_T2B", 10000);   \
        EXT_STR_ITEM_INFO_ZZP(                                                                                    \
            ok, si, offset, struct_t, printerr, NN_P3tft_T2Bv, UINT32, "NN_P3tft_T2Bv", "NN_P3tft_T2B");          \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, NN_P4tcl_T1BM, UINT32, "NN_P4tcl_T1BM", 50);    \
        EXT_STR_ITEM_INFO_ZZP(                                                                                    \
            ok, si, offset, struct_t, printerr, NN_P4tcl_T1BMI, UINT32, "NN_P4tcl_T1BMI", "NN_P4tcl_T1BM");       \
        EXT_STR_ITEM_INFO_ZZP(                                                                                    \
            ok, si, offset, struct_t, printerr, NN_P4tcl_T1BME, UINT32, "NN_P4tcl_T1BME", "NN_P4tcl_T1BM");       \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, NN_P4tcl_T1B, UINT32, "NN_P4tcl_T1B", 10000);   \
        EXT_STR_ITEM_INFO_ZZP(                                                                                    \
            ok, si, offset, struct_t, printerr, NN_P4tcl_T1Bv, UINT32, "NN_P4tcl_T1Bv", "NN_P4tcl_T1B");          \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, NN_P4tcl_T2BM, UINT32, "NN_P4tcl_T2BM", 50);    \
        EXT_STR_ITEM_INFO_ZZP(                                                                                    \
            ok, si, offset, struct_t, printerr, NN_P4tcl_T2BMI, UINT32, "NN_P4tcl_T2BMI", "NN_P4tcl_T2BM");       \
        EXT_STR_ITEM_INFO_ZZP(                                                                                    \
            ok, si, offset, struct_t, printerr, NN_P4tcl_T2BME, UINT32, "NN_P4tcl_T2BME", "NN_P4tcl_T2BM");       \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, NN_P4tcl_T2B, UINT32, "NN_P4tcl_T2B", 10000);   \
        EXT_STR_ITEM_INFO_ZZP(                                                                                    \
            ok, si, offset, struct_t, printerr, NN_P4tcl_T2Bv, UINT32, "NN_P4tcl_T2Bv", "NN_P4tcl_T2B");          \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, NN_P4tfl_T1BM, UINT32, "NN_P4tfl_T1BM", 50);    \
        EXT_STR_ITEM_INFO_ZZP(                                                                                    \
            ok, si, offset, struct_t, printerr, NN_P4tfl_T1BMI, UINT32, "NN_P4tfl_T1BMI", "NN_P4tfl_T1BM");       \
        EXT_STR_ITEM_INFO_ZZP(                                                                                    \
            ok, si, offset, struct_t, printerr, NN_P4tfl_T1BME, UINT32, "NN_P4tfl_T1BME", "NN_P4tfl_T1BM");       \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, NN_P4tfl_T1B, UINT32, "NN_P4tfl_T1B", 10000);   \
        EXT_STR_ITEM_INFO_ZZP(                                                                                    \
            ok, si, offset, struct_t, printerr, NN_P4tfl_T1Bv, UINT32, "NN_P4tfl_T1Bv", "NN_P4tfl_T1B");          \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, NN_P4tfl_T2BM, UINT32, "NN_P4tfl_T2BM", 50);    \
        EXT_STR_ITEM_INFO_ZZP(                                                                                    \
            ok, si, offset, struct_t, printerr, NN_P4tfl_T2BMI, UINT32, "NN_P4tfl_T2BMI", "NN_P4tfl_T2BM");       \
        EXT_STR_ITEM_INFO_ZZP(                                                                                    \
            ok, si, offset, struct_t, printerr, NN_P4tfl_T2BME, UINT32, "NN_P4tfl_T2BME", "NN_P4tfl_T2BM");       \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, NN_P4tfl_T2B, UINT32, "NN_P4tfl_T2B", 10000);   \
        EXT_STR_ITEM_INFO_ZZP(                                                                                    \
            ok, si, offset, struct_t, printerr, NN_P4tfl_T2Bv, UINT32, "NN_P4tfl_T2Bv", "NN_P4tfl_T2B");          \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, NN_P4tct_T1BM, UINT32, "NN_P4tct_T1BM", 50);    \
        EXT_STR_ITEM_INFO_ZZP(                                                                                    \
            ok, si, offset, struct_t, printerr, NN_P4tct_T1BMI, UINT32, "NN_P4tct_T1BMI", "NN_P4tct_T1BM");       \
        EXT_STR_ITEM_INFO_ZZP(                                                                                    \
            ok, si, offset, struct_t, printerr, NN_P4tct_T1BME, UINT32, "NN_P4tct_T1BME", "NN_P4tct_T1BM");       \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, NN_P4tct_T1B, UINT32, "NN_P4tct_T1B", 10000);   \
        EXT_STR_ITEM_INFO_ZZP(                                                                                    \
            ok, si, offset, struct_t, printerr, NN_P4tct_T1Bv, UINT32, "NN_P4tct_T1Bv", "NN_P4tct_T1B");          \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, NN_P4tct_T2BM, UINT32, "NN_P4tct_T2BM", 50);    \
        EXT_STR_ITEM_INFO_ZZP(                                                                                    \
            ok, si, offset, struct_t, printerr, NN_P4tct_T2BMI, UINT32, "NN_P4tct_T2BMI", "NN_P4tct_T2BM");       \
        EXT_STR_ITEM_INFO_ZZP(                                                                                    \
            ok, si, offset, struct_t, printerr, NN_P4tct_T2BME, UINT32, "NN_P4tct_T2BME", "NN_P4tct_T2BM");       \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, NN_P4tct_T2B, UINT32, "NN_P4tct_T2B", 10000);   \
        EXT_STR_ITEM_INFO_ZZP(                                                                                    \
            ok, si, offset, struct_t, printerr, NN_P4tct_T2Bv, UINT32, "NN_P4tct_T2Bv", "NN_P4tct_T2B");          \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, NN_P4tft_T1BM, UINT32, "NN_P4tft_T1BM", 50);    \
        EXT_STR_ITEM_INFO_ZZP(                                                                                    \
            ok, si, offset, struct_t, printerr, NN_P4tft_T1BMI, UINT32, "NN_P4tft_T1BMI", "NN_P4tft_T1BM");       \
        EXT_STR_ITEM_INFO_ZZP(                                                                                    \
            ok, si, offset, struct_t, printerr, NN_P4tft_T1BME, UINT32, "NN_P4tft_T1BME", "NN_P4tft_T1BM");       \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, NN_P4tft_T1B, UINT32, "NN_P4tft_T1B", 10000);   \
        EXT_STR_ITEM_INFO_ZZP(                                                                                    \
            ok, si, offset, struct_t, printerr, NN_P4tft_T1Bv, UINT32, "NN_P4tft_T1Bv", "NN_P4tft_T1B");          \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, NN_P4tft_T2BM, UINT32, "NN_P4tft_T2BM", 50);    \
        EXT_STR_ITEM_INFO_ZZP(                                                                                    \
            ok, si, offset, struct_t, printerr, NN_P4tft_T2BMI, UINT32, "NN_P4tft_T2BMI", "NN_P4tft_T2BM");       \
        EXT_STR_ITEM_INFO_ZZP(                                                                                    \
            ok, si, offset, struct_t, printerr, NN_P4tft_T2BME, UINT32, "NN_P4tft_T2BME", "NN_P4tft_T2BM");       \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, NN_P4tft_T2B, UINT32, "NN_P4tft_T2B", 10000);   \
        EXT_STR_ITEM_INFO_ZZP(                                                                                    \
            ok, si, offset, struct_t, printerr, NN_P4tft_T2Bv, UINT32, "NN_P4tft_T2Bv", "NN_P4tft_T2B");          \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, NN_P5tcl_T1BM, UINT32, "NN_P5tcl_T1BM", 50);    \
        EXT_STR_ITEM_INFO_ZZP(                                                                                    \
            ok, si, offset, struct_t, printerr, NN_P5tcl_T1BMI, UINT32, "NN_P5tcl_T1BMI", "NN_P5tcl_T1BM");       \
        EXT_STR_ITEM_INFO_ZZP(                                                                                    \
            ok, si, offset, struct_t, printerr, NN_P5tcl_T1BME, UINT32, "NN_P5tcl_T1BME", "NN_P5tcl_T1BM");       \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, NN_P5tcl_T1B, UINT32, "NN_P5tcl_T1B", 10000);   \
        EXT_STR_ITEM_INFO_ZZP(                                                                                    \
            ok, si, offset, struct_t, printerr, NN_P5tcl_T1Bv, UINT32, "NN_P5tcl_T1Bv", "NN_P5tcl_T1B");          \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, NN_P5tcl_T2BM, UINT32, "NN_P5tcl_T2BM", 50);    \
        EXT_STR_ITEM_INFO_ZZP(                                                                                    \
            ok, si, offset, struct_t, printerr, NN_P5tcl_T2BMI, UINT32, "NN_P5tcl_T2BMI", "NN_P5tcl_T2BM");       \
        EXT_STR_ITEM_INFO_ZZP(                                                                                    \
            ok, si, offset, struct_t, printerr, NN_P5tcl_T2BME, UINT32, "NN_P5tcl_T2BME", "NN_P5tcl_T2BM");       \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, NN_P5tcl_T2B, UINT32, "NN_P5tcl_T2B", 10000);   \
        EXT_STR_ITEM_INFO_ZZP(                                                                                    \
            ok, si, offset, struct_t, printerr, NN_P5tcl_T2Bv, UINT32, "NN_P5tcl_T2Bv", "NN_P5tcl_T2B");          \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, NN_P5tfl_T1BM, UINT32, "NN_P5tfl_T1BM", 50);    \
        EXT_STR_ITEM_INFO_ZZP(                                                                                    \
            ok, si, offset, struct_t, printerr, NN_P5tfl_T1BMI, UINT32, "NN_P5tfl_T1BMI", "NN_P5tfl_T1BM");       \
        EXT_STR_ITEM_INFO_ZZP(                                                                                    \
            ok, si, offset, struct_t, printerr, NN_P5tfl_T1BME, UINT32, "NN_P5tfl_T1BME", "NN_P5tfl_T1BM");       \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, NN_P5tfl_T1B, UINT32, "NN_P5tfl_T1B", 10000);   \
        EXT_STR_ITEM_INFO_ZZP(                                                                                    \
            ok, si, offset, struct_t, printerr, NN_P5tfl_T1Bv, UINT32, "NN_P5tfl_T1Bv", "NN_P5tfl_T1B");          \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, NN_P5tfl_T2BM, UINT32, "NN_P5tfl_T2BM", 50);    \
        EXT_STR_ITEM_INFO_ZZP(                                                                                    \
            ok, si, offset, struct_t, printerr, NN_P5tfl_T2BMI, UINT32, "NN_P5tfl_T2BMI", "NN_P5tfl_T2BM");       \
        EXT_STR_ITEM_INFO_ZZP(                                                                                    \
            ok, si, offset, struct_t, printerr, NN_P5tfl_T2BME, UINT32, "NN_P5tfl_T2BME", "NN_P5tfl_T2BM");       \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, NN_P5tfl_T2B, UINT32, "NN_P5tfl_T2B", 10000);   \
        EXT_STR_ITEM_INFO_ZZP(                                                                                    \
            ok, si, offset, struct_t, printerr, NN_P5tfl_T2Bv, UINT32, "NN_P5tfl_T2Bv", "NN_P5tfl_T2B");          \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, NN_P5tct_T1BM, UINT32, "NN_P5tct_T1BM", 50);    \
        EXT_STR_ITEM_INFO_ZZP(                                                                                    \
            ok, si, offset, struct_t, printerr, NN_P5tct_T1BMI, UINT32, "NN_P5tct_T1BMI", "NN_P5tct_T1BM");       \
        EXT_STR_ITEM_INFO_ZZP(                                                                                    \
            ok, si, offset, struct_t, printerr, NN_P5tct_T1BME, UINT32, "NN_P5tct_T1BME", "NN_P5tct_T1BM");       \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, NN_P5tct_T1B, UINT32, "NN_P5tct_T1B", 10000);   \
        EXT_STR_ITEM_INFO_ZZP(                                                                                    \
            ok, si, offset, struct_t, printerr, NN_P5tct_T1Bv, UINT32, "NN_P5tct_T1Bv", "NN_P5tct_T1B");          \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, NN_P5tct_T2BM, UINT32, "NN_P5tct_T2BM", 50);    \
        EXT_STR_ITEM_INFO_ZZP(                                                                                    \
            ok, si, offset, struct_t, printerr, NN_P5tct_T2BMI, UINT32, "NN_P5tct_T2BMI", "NN_P5tct_T2BM");       \
        EXT_STR_ITEM_INFO_ZZP(                                                                                    \
            ok, si, offset, struct_t, printerr, NN_P5tct_T2BME, UINT32, "NN_P5tct_T2BME", "NN_P5tct_T2BM");       \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, NN_P5tct_T2B, UINT32, "NN_P5tct_T2B", 10000);   \
        EXT_STR_ITEM_INFO_ZZP(                                                                                    \
            ok, si, offset, struct_t, printerr, NN_P5tct_T2Bv, UINT32, "NN_P5tct_T2Bv", "NN_P5tct_T2B");          \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, NN_P5tft_T1BM, UINT32, "NN_P5tft_T1BM", 50);    \
        EXT_STR_ITEM_INFO_ZZP(                                                                                    \
            ok, si, offset, struct_t, printerr, NN_P5tft_T1BMI, UINT32, "NN_P5tft_T1BMI", "NN_P5tft_T1BM");       \
        EXT_STR_ITEM_INFO_ZZP(                                                                                    \
            ok, si, offset, struct_t, printerr, NN_P5tft_T1BME, UINT32, "NN_P5tft_T1BME", "NN_P5tft_T1BM");       \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, NN_P5tft_T1B, UINT32, "NN_P5tft_T1B", 10000);   \
        EXT_STR_ITEM_INFO_ZZP(                                                                                    \
            ok, si, offset, struct_t, printerr, NN_P5tft_T1Bv, UINT32, "NN_P5tft_T1Bv", "NN_P5tft_T1B");          \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, NN_P5tft_T2BM, UINT32, "NN_P5tft_T2BM", 50);    \
        EXT_STR_ITEM_INFO_ZZP(                                                                                    \
            ok, si, offset, struct_t, printerr, NN_P5tft_T2BMI, UINT32, "NN_P5tft_T2BMI", "NN_P5tft_T2BM");       \
        EXT_STR_ITEM_INFO_ZZP(                                                                                    \
            ok, si, offset, struct_t, printerr, NN_P5tft_T2BME, UINT32, "NN_P5tft_T2BME", "NN_P5tft_T2BM");       \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, NN_P5tft_T2B, UINT32, "NN_P5tft_T2B", 10000);   \
        EXT_STR_ITEM_INFO_ZZP(                                                                                    \
            ok, si, offset, struct_t, printerr, NN_P5tft_T2Bv, UINT32, "NN_P5tft_T2Bv", "NN_P5tft_T2B");          \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, NN_P6tcl_T1BM, UINT32, "NN_P6tcl_T1BM", 50);    \
        EXT_STR_ITEM_INFO_ZZP(                                                                                    \
            ok, si, offset, struct_t, printerr, NN_P6tcl_T1BMI, UINT32, "NN_P6tcl_T1BMI", "NN_P6tcl_T1BM");       \
        EXT_STR_ITEM_INFO_ZZP(                                                                                    \
            ok, si, offset, struct_t, printerr, NN_P6tcl_T1BME, UINT32, "NN_P6tcl_T1BME", "NN_P6tcl_T1BM");       \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, NN_P6tcl_T1B, UINT32, "NN_P6tcl_T1B", 10000);   \
        EXT_STR_ITEM_INFO_ZZP(                                                                                    \
            ok, si, offset, struct_t, printerr, NN_P6tcl_T1Bv, UINT32, "NN_P6tcl_T1Bv", "NN_P6tcl_T1B");          \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, NN_P6tcl_T2BM, UINT32, "NN_P6tcl_T2BM", 50);    \
        EXT_STR_ITEM_INFO_ZZP(                                                                                    \
            ok, si, offset, struct_t, printerr, NN_P6tcl_T2BMI, UINT32, "NN_P6tcl_T2BMI", "NN_P6tcl_T2BM");       \
        EXT_STR_ITEM_INFO_ZZP(                                                                                    \
            ok, si, offset, struct_t, printerr, NN_P6tcl_T2BME, UINT32, "NN_P6tcl_T2BME", "NN_P6tcl_T2BM");       \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, NN_P6tcl_T2B, UINT32, "NN_P6tcl_T2B", 10000);   \
        EXT_STR_ITEM_INFO_ZZP(                                                                                    \
            ok, si, offset, struct_t, printerr, NN_P6tcl_T2Bv, UINT32, "NN_P6tcl_T2Bv", "NN_P6tcl_T2B");          \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, NN_P6tfl_T1BM, UINT32, "NN_P6tfl_T1BM", 50);    \
        EXT_STR_ITEM_INFO_ZZP(                                                                                    \
            ok, si, offset, struct_t, printerr, NN_P6tfl_T1BMI, UINT32, "NN_P6tfl_T1BMI", "NN_P6tfl_T1BM");       \
        EXT_STR_ITEM_INFO_ZZP(                                                                                    \
            ok, si, offset, struct_t, printerr, NN_P6tfl_T1BME, UINT32, "NN_P6tfl_T1BME", "NN_P6tfl_T1BM");       \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, NN_P6tfl_T1B, UINT32, "NN_P6tfl_T1B", 10000);   \
        EXT_STR_ITEM_INFO_ZZP(                                                                                    \
            ok, si, offset, struct_t, printerr, NN_P6tfl_T1Bv, UINT32, "NN_P6tfl_T1Bv", "NN_P6tfl_T1B");          \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, NN_P6tfl_T2BM, UINT32, "NN_P6tfl_T2BM", 50);    \
        EXT_STR_ITEM_INFO_ZZP(                                                                                    \
            ok, si, offset, struct_t, printerr, NN_P6tfl_T2BMI, UINT32, "NN_P6tfl_T2BMI", "NN_P6tfl_T2BM");       \
        EXT_STR_ITEM_INFO_ZZP(                                                                                    \
            ok, si, offset, struct_t, printerr, NN_P6tfl_T2BME, UINT32, "NN_P6tfl_T2BME", "NN_P6tfl_T2BM");       \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, NN_P6tfl_T2B, UINT32, "NN_P6tfl_T2B", 10000);   \
        EXT_STR_ITEM_INFO_ZZP(                                                                                    \
            ok, si, offset, struct_t, printerr, NN_P6tfl_T2Bv, UINT32, "NN_P6tfl_T2Bv", "NN_P6tfl_T2B");          \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, NN_P6tct_T1BM, UINT32, "NN_P6tct_T1BM", 50);    \
        EXT_STR_ITEM_INFO_ZZP(                                                                                    \
            ok, si, offset, struct_t, printerr, NN_P6tct_T1BMI, UINT32, "NN_P6tct_T1BMI", "NN_P6tct_T1BM");       \
        EXT_STR_ITEM_INFO_ZZP(                                                                                    \
            ok, si, offset, struct_t, printerr, NN_P6tct_T1BME, UINT32, "NN_P6tct_T1BME", "NN_P6tct_T1BM");       \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, NN_P6tct_T1B, UINT32, "NN_P6tct_T1B", 10000);   \
        EXT_STR_ITEM_INFO_ZZP(                                                                                    \
            ok, si, offset, struct_t, printerr, NN_P6tct_T1Bv, UINT32, "NN_P6tct_T1Bv", "NN_P6tct_T1B");          \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, NN_P6tct_T2BM, UINT32, "NN_P6tct_T2BM", 50);    \
        EXT_STR_ITEM_INFO_ZZP(                                                                                    \
            ok, si, offset, struct_t, printerr, NN_P6tct_T2BMI, UINT32, "NN_P6tct_T2BMI", "NN_P6tct_T2BM");       \
        EXT_STR_ITEM_INFO_ZZP(                                                                                    \
            ok, si, offset, struct_t, printerr, NN_P6tct_T2BME, UINT32, "NN_P6tct_T2BME", "NN_P6tct_T2BM");       \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, NN_P6tct_T2B, UINT32, "NN_P6tct_T2B", 10000);   \
        EXT_STR_ITEM_INFO_ZZP(                                                                                    \
            ok, si, offset, struct_t, printerr, NN_P6tct_T2Bv, UINT32, "NN_P6tct_T2Bv", "NN_P6tct_T2B");          \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, NN_P6tft_T1BM, UINT32, "NN_P6tft_T1BM", 50);    \
        EXT_STR_ITEM_INFO_ZZP(                                                                                    \
            ok, si, offset, struct_t, printerr, NN_P6tft_T1BMI, UINT32, "NN_P6tft_T1BMI", "NN_P6tft_T1BM");       \
        EXT_STR_ITEM_INFO_ZZP(                                                                                    \
            ok, si, offset, struct_t, printerr, NN_P6tft_T1BME, UINT32, "NN_P6tft_T1BME", "NN_P6tft_T1BM");       \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, NN_P6tft_T1B, UINT32, "NN_P6tft_T1B", 10000);   \
        EXT_STR_ITEM_INFO_ZZP(                                                                                    \
            ok, si, offset, struct_t, printerr, NN_P6tft_T1Bv, UINT32, "NN_P6tft_T1Bv", "NN_P6tft_T1B");          \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, NN_P6tft_T2BM, UINT32, "NN_P6tft_T2BM", 50);    \
        EXT_STR_ITEM_INFO_ZZP(                                                                                    \
            ok, si, offset, struct_t, printerr, NN_P6tft_T2BMI, UINT32, "NN_P6tft_T2BMI", "NN_P6tft_T2BM");       \
        EXT_STR_ITEM_INFO_ZZP(                                                                                    \
            ok, si, offset, struct_t, printerr, NN_P6tft_T2BME, UINT32, "NN_P6tft_T2BME", "NN_P6tft_T2BM");       \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, NN_P6tft_T2B, UINT32, "NN_P6tft_T2B", 10000);   \
        EXT_STR_ITEM_INFO_ZZP(                                                                                    \
            ok, si, offset, struct_t, printerr, NN_P6tft_T2Bv, UINT32, "NN_P6tft_T2Bv", "NN_P6tft_T2B");          \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, NN_P7tcl_T1BM, UINT32, "NN_P7tcl_T1BM", 50);    \
        EXT_STR_ITEM_INFO_ZZP(                                                                                    \
            ok, si, offset, struct_t, printerr, NN_P7tcl_T1BMI, UINT32, "NN_P7tcl_T1BMI", "NN_P7tcl_T1BM");       \
        EXT_STR_ITEM_INFO_ZZP(                                                                                    \
            ok, si, offset, struct_t, printerr, NN_P7tcl_T1BME, UINT32, "NN_P7tcl_T1BME", "NN_P7tcl_T1BM");       \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, NN_P7tcl_T1B, UINT32, "NN_P7tcl_T1B", 10000);   \
        EXT_STR_ITEM_INFO_ZZP(                                                                                    \
            ok, si, offset, struct_t, printerr, NN_P7tcl_T1Bv, UINT32, "NN_P7tcl_T1Bv", "NN_P7tcl_T1B");          \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, NN_P7tcl_T2BM, UINT32, "NN_P7tcl_T2BM", 50);    \
        EXT_STR_ITEM_INFO_ZZP(                                                                                    \
            ok, si, offset, struct_t, printerr, NN_P7tcl_T2BMI, UINT32, "NN_P7tcl_T2BMI", "NN_P7tcl_T2BM");       \
        EXT_STR_ITEM_INFO_ZZP(                                                                                    \
            ok, si, offset, struct_t, printerr, NN_P7tcl_T2BME, UINT32, "NN_P7tcl_T2BME", "NN_P7tcl_T2BM");       \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, NN_P7tcl_T2B, UINT32, "NN_P7tcl_T2B", 10000);   \
        EXT_STR_ITEM_INFO_ZZP(                                                                                    \
            ok, si, offset, struct_t, printerr, NN_P7tcl_T2Bv, UINT32, "NN_P7tcl_T2Bv", "NN_P7tcl_T2B");          \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, NN_P7tfl_T1BM, UINT32, "NN_P7tfl_T1BM", 50);    \
        EXT_STR_ITEM_INFO_ZZP(                                                                                    \
            ok, si, offset, struct_t, printerr, NN_P7tfl_T1BMI, UINT32, "NN_P7tfl_T1BMI", "NN_P7tfl_T1BM");       \
        EXT_STR_ITEM_INFO_ZZP(                                                                                    \
            ok, si, offset, struct_t, printerr, NN_P7tfl_T1BME, UINT32, "NN_P7tfl_T1BME", "NN_P7tfl_T1BM");       \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, NN_P7tfl_T1B, UINT32, "NN_P7tfl_T1B", 10000);   \
        EXT_STR_ITEM_INFO_ZZP(                                                                                    \
            ok, si, offset, struct_t, printerr, NN_P7tfl_T1Bv, UINT32, "NN_P7tfl_T1Bv", "NN_P7tfl_T1B");          \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, NN_P7tfl_T2BM, UINT32, "NN_P7tfl_T2BM", 50);    \
        EXT_STR_ITEM_INFO_ZZP(                                                                                    \
            ok, si, offset, struct_t, printerr, NN_P7tfl_T2BMI, UINT32, "NN_P7tfl_T2BMI", "NN_P7tfl_T2BM");       \
        EXT_STR_ITEM_INFO_ZZP(                                                                                    \
            ok, si, offset, struct_t, printerr, NN_P7tfl_T2BME, UINT32, "NN_P7tfl_T2BME", "NN_P7tfl_T2BM");       \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, NN_P7tfl_T2B, UINT32, "NN_P7tfl_T2B", 10000);   \
        EXT_STR_ITEM_INFO_ZZP(                                                                                    \
            ok, si, offset, struct_t, printerr, NN_P7tfl_T2Bv, UINT32, "NN_P7tfl_T2Bv", "NN_P7tfl_T2B");          \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, NN_P7tct_T1BM, UINT32, "NN_P7tct_T1BM", 50);    \
        EXT_STR_ITEM_INFO_ZZP(                                                                                    \
            ok, si, offset, struct_t, printerr, NN_P7tct_T1BMI, UINT32, "NN_P7tct_T1BMI", "NN_P7tct_T1BM");       \
        EXT_STR_ITEM_INFO_ZZP(                                                                                    \
            ok, si, offset, struct_t, printerr, NN_P7tct_T1BME, UINT32, "NN_P7tct_T1BME", "NN_P7tct_T1BM");       \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, NN_P7tct_T1B, UINT32, "NN_P7tct_T1B", 10000);   \
        EXT_STR_ITEM_INFO_ZZP(                                                                                    \
            ok, si, offset, struct_t, printerr, NN_P7tct_T1Bv, UINT32, "NN_P7tct_T1Bv", "NN_P7tct_T1B");          \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, NN_P7tct_T2BM, UINT32, "NN_P7tct_T2BM", 50);    \
        EXT_STR_ITEM_INFO_ZZP(                                                                                    \
            ok, si, offset, struct_t, printerr, NN_P7tct_T2BMI, UINT32, "NN_P7tct_T2BMI", "NN_P7tct_T2BM");       \
        EXT_STR_ITEM_INFO_ZZP(                                                                                    \
            ok, si, offset, struct_t, printerr, NN_P7tct_T2BME, UINT32, "NN_P7tct_T2BME", "NN_P7tct_T2BM");       \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, NN_P7tct_T2B, UINT32, "NN_P7tct_T2B", 10000);   \
        EXT_STR_ITEM_INFO_ZZP(                                                                                    \
            ok, si, offset, struct_t, printerr, NN_P7tct_T2Bv, UINT32, "NN_P7tct_T2Bv", "NN_P7tct_T2B");          \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, NN_P7tft_T1BM, UINT32, "NN_P7tft_T1BM", 50);    \
        EXT_STR_ITEM_INFO_ZZP(                                                                                    \
            ok, si, offset, struct_t, printerr, NN_P7tft_T1BMI, UINT32, "NN_P7tft_T1BMI", "NN_P7tft_T1BM");       \
        EXT_STR_ITEM_INFO_ZZP(                                                                                    \
            ok, si, offset, struct_t, printerr, NN_P7tft_T1BME, UINT32, "NN_P7tft_T1BME", "NN_P7tft_T1BM");       \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, NN_P7tft_T1B, UINT32, "NN_P7tft_T1B", 10000);   \
        EXT_STR_ITEM_INFO_ZZP(                                                                                    \
            ok, si, offset, struct_t, printerr, NN_P7tft_T1Bv, UINT32, "NN_P7tft_T1Bv", "NN_P7tft_T1B");          \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, NN_P7tft_T2BM, UINT32, "NN_P7tft_T2BM", 50);    \
        EXT_STR_ITEM_INFO_ZZP(                                                                                    \
            ok, si, offset, struct_t, printerr, NN_P7tft_T2BMI, UINT32, "NN_P7tft_T2BMI", "NN_P7tft_T2BM");       \
        EXT_STR_ITEM_INFO_ZZP(                                                                                    \
            ok, si, offset, struct_t, printerr, NN_P7tft_T2BME, UINT32, "NN_P7tft_T2BME", "NN_P7tft_T2BM");       \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, NN_P7tft_T2B, UINT32, "NN_P7tft_T2B", 10000);   \
        EXT_STR_ITEM_INFO_ZZP(                                                                                    \
            ok, si, offset, struct_t, printerr, NN_P7tft_T2Bv, UINT32, "NN_P7tft_T2Bv", "NN_P7tft_T2B");          \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, NN_P8tcl_T1BM, UINT32, "NN_P8tcl_T1BM", 50);    \
        EXT_STR_ITEM_INFO_ZZP(                                                                                    \
            ok, si, offset, struct_t, printerr, NN_P8tcl_T1BMI, UINT32, "NN_P8tcl_T1BMI", "NN_P8tcl_T1BM");       \
        EXT_STR_ITEM_INFO_ZZP(                                                                                    \
            ok, si, offset, struct_t, printerr, NN_P8tcl_T1BME, UINT32, "NN_P8tcl_T1BME", "NN_P8tcl_T1BM");       \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, NN_P8tcl_T1B, UINT32, "NN_P8tcl_T1B", 10000);   \
        EXT_STR_ITEM_INFO_ZZP(                                                                                    \
            ok, si, offset, struct_t, printerr, NN_P8tcl_T1Bv, UINT32, "NN_P8tcl_T1Bv", "NN_P8tcl_T1B");          \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, NN_P8tcl_T2BM, UINT32, "NN_P8tcl_T2BM", 50);    \
        EXT_STR_ITEM_INFO_ZZP(                                                                                    \
            ok, si, offset, struct_t, printerr, NN_P8tcl_T2BMI, UINT32, "NN_P8tcl_T2BMI", "NN_P8tcl_T2BM");       \
        EXT_STR_ITEM_INFO_ZZP(                                                                                    \
            ok, si, offset, struct_t, printerr, NN_P8tcl_T2BME, UINT32, "NN_P8tcl_T2BME", "NN_P8tcl_T2BM");       \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, NN_P8tcl_T2B, UINT32, "NN_P8tcl_T2B", 10000);   \
        EXT_STR_ITEM_INFO_ZZP(                                                                                    \
            ok, si, offset, struct_t, printerr, NN_P8tcl_T2Bv, UINT32, "NN_P8tcl_T2Bv", "NN_P8tcl_T2B");          \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, NN_P8tfl_T1BM, UINT32, "NN_P8tfl_T1BM", 50);    \
        EXT_STR_ITEM_INFO_ZZP(                                                                                    \
            ok, si, offset, struct_t, printerr, NN_P8tfl_T1BMI, UINT32, "NN_P8tfl_T1BMI", "NN_P8tfl_T1BM");       \
        EXT_STR_ITEM_INFO_ZZP(                                                                                    \
            ok, si, offset, struct_t, printerr, NN_P8tfl_T1BME, UINT32, "NN_P8tfl_T1BME", "NN_P8tfl_T1BM");       \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, NN_P8tfl_T1B, UINT32, "NN_P8tfl_T1B", 10000);   \
        EXT_STR_ITEM_INFO_ZZP(                                                                                    \
            ok, si, offset, struct_t, printerr, NN_P8tfl_T1Bv, UINT32, "NN_P8tfl_T1Bv", "NN_P8tfl_T1B");          \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, NN_P8tfl_T2BM, UINT32, "NN_P8tfl_T2BM", 50);    \
        EXT_STR_ITEM_INFO_ZZP(                                                                                    \
            ok, si, offset, struct_t, printerr, NN_P8tfl_T2BMI, UINT32, "NN_P8tfl_T2BMI", "NN_P8tfl_T2BM");       \
        EXT_STR_ITEM_INFO_ZZP(                                                                                    \
            ok, si, offset, struct_t, printerr, NN_P8tfl_T2BME, UINT32, "NN_P8tfl_T2BME", "NN_P8tfl_T2BM");       \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, NN_P8tfl_T2B, UINT32, "NN_P8tfl_T2B", 10000);   \
        EXT_STR_ITEM_INFO_ZZP(                                                                                    \
            ok, si, offset, struct_t, printerr, NN_P8tfl_T2Bv, UINT32, "NN_P8tfl_T2Bv", "NN_P8tfl_T2B");          \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, NN_P8tct_T1BM, UINT32, "NN_P8tct_T1BM", 50);    \
        EXT_STR_ITEM_INFO_ZZP(                                                                                    \
            ok, si, offset, struct_t, printerr, NN_P8tct_T1BMI, UINT32, "NN_P8tct_T1BMI", "NN_P8tct_T1BM");       \
        EXT_STR_ITEM_INFO_ZZP(                                                                                    \
            ok, si, offset, struct_t, printerr, NN_P8tct_T1BME, UINT32, "NN_P8tct_T1BME", "NN_P8tct_T1BM");       \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, NN_P8tct_T1B, UINT32, "NN_P8tct_T1B", 10000);   \
        EXT_STR_ITEM_INFO_ZZP(                                                                                    \
            ok, si, offset, struct_t, printerr, NN_P8tct_T1Bv, UINT32, "NN_P8tct_T1Bv", "NN_P8tct_T1B");          \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, NN_P8tct_T2BM, UINT32, "NN_P8tct_T2BM", 50);    \
        EXT_STR_ITEM_INFO_ZZP(                                                                                    \
            ok, si, offset, struct_t, printerr, NN_P8tct_T2BMI, UINT32, "NN_P8tct_T2BMI", "NN_P8tct_T2BM");       \
        EXT_STR_ITEM_INFO_ZZP(                                                                                    \
            ok, si, offset, struct_t, printerr, NN_P8tct_T2BME, UINT32, "NN_P8tct_T2BME", "NN_P8tct_T2BM");       \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, NN_P8tct_T2B, UINT32, "NN_P8tct_T2B", 10000);   \
        EXT_STR_ITEM_INFO_ZZP(                                                                                    \
            ok, si, offset, struct_t, printerr, NN_P8tct_T2Bv, UINT32, "NN_P8tct_T2Bv", "NN_P8tct_T2B");          \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, NN_P8tft_T1BM, UINT32, "NN_P8tft_T1BM", 50);    \
        EXT_STR_ITEM_INFO_ZZP(                                                                                    \
            ok, si, offset, struct_t, printerr, NN_P8tft_T1BMI, UINT32, "NN_P8tft_T1BMI", "NN_P8tft_T1BM");       \
        EXT_STR_ITEM_INFO_ZZP(                                                                                    \
            ok, si, offset, struct_t, printerr, NN_P8tft_T1BME, UINT32, "NN_P8tft_T1BME", "NN_P8tft_T1BM");       \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, NN_P8tft_T1B, UINT32, "NN_P8tft_T1B", 10000);   \
        EXT_STR_ITEM_INFO_ZZP(                                                                                    \
            ok, si, offset, struct_t, printerr, NN_P8tft_T1Bv, UINT32, "NN_P8tft_T1Bv", "NN_P8tft_T1B");          \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, NN_P8tft_T2BM, UINT32, "NN_P8tft_T2BM", 50);    \
        EXT_STR_ITEM_INFO_ZZP(                                                                                    \
            ok, si, offset, struct_t, printerr, NN_P8tft_T2BMI, UINT32, "NN_P8tft_T2BMI", "NN_P8tft_T2BM");       \
        EXT_STR_ITEM_INFO_ZZP(                                                                                    \
            ok, si, offset, struct_t, printerr, NN_P8tft_T2BME, UINT32, "NN_P8tft_T2BME", "NN_P8tft_T2BM");       \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, NN_P8tft_T2B, UINT32, "NN_P8tft_T2B", 10000);   \
        EXT_STR_ITEM_INFO_ZZP(                                                                                    \
            ok, si, offset, struct_t, printerr, NN_P8tft_T2Bv, UINT32, "NN_P8tft_T2Bv", "NN_P8tft_T2B");          \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, NN_P9tcl_T1BM, UINT32, "NN_P9tcl_T1BM", 50);    \
        EXT_STR_ITEM_INFO_ZZP(                                                                                    \
            ok, si, offset, struct_t, printerr, NN_P9tcl_T1BMI, UINT32, "NN_P9tcl_T1BMI", "NN_P9tcl_T1BM");       \
        EXT_STR_ITEM_INFO_ZZP(                                                                                    \
            ok, si, offset, struct_t, printerr, NN_P9tcl_T1BME, UINT32, "NN_P9tcl_T1BME", "NN_P9tcl_T1BM");       \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, NN_P9tcl_T1B, UINT32, "NN_P9tcl_T1B", 10000);   \
        EXT_STR_ITEM_INFO_ZZP(                                                                                    \
            ok, si, offset, struct_t, printerr, NN_P9tcl_T1Bv, UINT32, "NN_P9tcl_T1Bv", "NN_P9tcl_T1B");          \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, NN_P9tcl_T2BM, UINT32, "NN_P9tcl_T2BM", 50);    \
        EXT_STR_ITEM_INFO_ZZP(                                                                                    \
            ok, si, offset, struct_t, printerr, NN_P9tcl_T2BMI, UINT32, "NN_P9tcl_T2BMI", "NN_P9tcl_T2BM");       \
        EXT_STR_ITEM_INFO_ZZP(                                                                                    \
            ok, si, offset, struct_t, printerr, NN_P9tcl_T2BME, UINT32, "NN_P9tcl_T2BME", "NN_P9tcl_T2BM");       \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, NN_P9tcl_T2B, UINT32, "NN_P9tcl_T2B", 10000);   \
        EXT_STR_ITEM_INFO_ZZP(                                                                                    \
            ok, si, offset, struct_t, printerr, NN_P9tcl_T2Bv, UINT32, "NN_P9tcl_T2Bv", "NN_P9tcl_T2B");          \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, NN_P9tfl_T1BM, UINT32, "NN_P9tfl_T1BM", 50);    \
        EXT_STR_ITEM_INFO_ZZP(                                                                                    \
            ok, si, offset, struct_t, printerr, NN_P9tfl_T1BMI, UINT32, "NN_P9tfl_T1BMI", "NN_P9tfl_T1BM");       \
        EXT_STR_ITEM_INFO_ZZP(                                                                                    \
            ok, si, offset, struct_t, printerr, NN_P9tfl_T1BME, UINT32, "NN_P9tfl_T1BME", "NN_P9tfl_T1BM");       \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, NN_P9tfl_T1B, UINT32, "NN_P9tfl_T1B", 10000);   \
        EXT_STR_ITEM_INFO_ZZP(                                                                                    \
            ok, si, offset, struct_t, printerr, NN_P9tfl_T1Bv, UINT32, "NN_P9tfl_T1Bv", "NN_P9tfl_T1B");          \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, NN_P9tfl_T2BM, UINT32, "NN_P9tfl_T2BM", 50);    \
        EXT_STR_ITEM_INFO_ZZP(                                                                                    \
            ok, si, offset, struct_t, printerr, NN_P9tfl_T2BMI, UINT32, "NN_P9tfl_T2BMI", "NN_P9tfl_T2BM");       \
        EXT_STR_ITEM_INFO_ZZP(                                                                                    \
            ok, si, offset, struct_t, printerr, NN_P9tfl_T2BME, UINT32, "NN_P9tfl_T2BME", "NN_P9tfl_T2BM");       \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, NN_P9tfl_T2B, UINT32, "NN_P9tfl_T2B", 10000);   \
        EXT_STR_ITEM_INFO_ZZP(                                                                                    \
            ok, si, offset, struct_t, printerr, NN_P9tfl_T2Bv, UINT32, "NN_P9tfl_T2Bv", "NN_P9tfl_T2B");          \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, NN_P9tct_T1BM, UINT32, "NN_P9tct_T1BM", 50);    \
        EXT_STR_ITEM_INFO_ZZP(                                                                                    \
            ok, si, offset, struct_t, printerr, NN_P9tct_T1BMI, UINT32, "NN_P9tct_T1BMI", "NN_P9tct_T1BM");       \
        EXT_STR_ITEM_INFO_ZZP(                                                                                    \
            ok, si, offset, struct_t, printerr, NN_P9tct_T1BME, UINT32, "NN_P9tct_T1BME", "NN_P9tct_T1BM");       \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, NN_P9tct_T1B, UINT32, "NN_P9tct_T1B", 10000);   \
        EXT_STR_ITEM_INFO_ZZP(                                                                                    \
            ok, si, offset, struct_t, printerr, NN_P9tct_T1Bv, UINT32, "NN_P9tct_T1Bv", "NN_P9tct_T1B");          \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, NN_P9tct_T2BM, UINT32, "NN_P9tct_T2BM", 50);    \
        EXT_STR_ITEM_INFO_ZZP(                                                                                    \
            ok, si, offset, struct_t, printerr, NN_P9tct_T2BMI, UINT32, "NN_P9tct_T2BMI", "NN_P9tct_T2BM");       \
        EXT_STR_ITEM_INFO_ZZP(                                                                                    \
            ok, si, offset, struct_t, printerr, NN_P9tct_T2BME, UINT32, "NN_P9tct_T2BME", "NN_P9tct_T2BM");       \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, NN_P9tct_T2B, UINT32, "NN_P9tct_T2B", 10000);   \
        EXT_STR_ITEM_INFO_ZZP(                                                                                    \
            ok, si, offset, struct_t, printerr, NN_P9tct_T2Bv, UINT32, "NN_P9tct_T2Bv", "NN_P9tct_T2B");          \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, NN_P9tft_T1BM, UINT32, "NN_P9tft_T1BM", 50);    \
        EXT_STR_ITEM_INFO_ZZP(                                                                                    \
            ok, si, offset, struct_t, printerr, NN_P9tft_T1BMI, UINT32, "NN_P9tft_T1BMI", "NN_P9tft_T1BM");       \
        EXT_STR_ITEM_INFO_ZZP(                                                                                    \
            ok, si, offset, struct_t, printerr, NN_P9tft_T1BME, UINT32, "NN_P9tft_T1BME", "NN_P9tft_T1BM");       \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, NN_P9tft_T1B, UINT32, "NN_P9tft_T1B", 10000);   \
        EXT_STR_ITEM_INFO_ZZP(                                                                                    \
            ok, si, offset, struct_t, printerr, NN_P9tft_T1Bv, UINT32, "NN_P9tft_T1Bv", "NN_P9tft_T1B");          \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, NN_P9tft_T2BM, UINT32, "NN_P9tft_T2BM", 50);    \
        EXT_STR_ITEM_INFO_ZZP(                                                                                    \
            ok, si, offset, struct_t, printerr, NN_P9tft_T2BMI, UINT32, "NN_P9tft_T2BMI", "NN_P9tft_T2BM");       \
        EXT_STR_ITEM_INFO_ZZP(                                                                                    \
            ok, si, offset, struct_t, printerr, NN_P9tft_T2BME, UINT32, "NN_P9tft_T2BME", "NN_P9tft_T2BM");       \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, NN_P9tft_T2B, UINT32, "NN_P9tft_T2B", 10000);   \
        EXT_STR_ITEM_INFO_ZZP(                                                                                    \
            ok, si, offset, struct_t, printerr, NN_P9tft_T2Bv, UINT32, "NN_P9tft_T2Bv", "NN_P9tft_T2B");          \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, NN_P10tcl_T1BM, UINT32, "NN_P10tcl_T1BM", 50);  \
        EXT_STR_ITEM_INFO_ZZP(                                                                                    \
            ok, si, offset, struct_t, printerr, NN_P10tcl_T1BMI, UINT32, "NN_P10tcl_T1BMI", "NN_P10tcl_T1BM");    \
        EXT_STR_ITEM_INFO_ZZP(                                                                                    \
            ok, si, offset, struct_t, printerr, NN_P10tcl_T1BME, UINT32, "NN_P10tcl_T1BME", "NN_P10tcl_T1BM");    \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, NN_P10tcl_T1B, UINT32, "NN_P10tcl_T1B", 10000); \
        EXT_STR_ITEM_INFO_ZZP(                                                                                    \
            ok, si, offset, struct_t, printerr, NN_P10tcl_T1Bv, UINT32, "NN_P10tcl_T1Bv", "NN_P10tcl_T1B");       \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, NN_P10tcl_T2BM, UINT32, "NN_P10tcl_T2BM", 50);  \
        EXT_STR_ITEM_INFO_ZZP(                                                                                    \
            ok, si, offset, struct_t, printerr, NN_P10tcl_T2BMI, UINT32, "NN_P10tcl_T2BMI", "NN_P10tcl_T2BM");    \
        EXT_STR_ITEM_INFO_ZZP(                                                                                    \
            ok, si, offset, struct_t, printerr, NN_P10tcl_T2BME, UINT32, "NN_P10tcl_T2BME", "NN_P10tcl_T2BM");    \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, NN_P10tcl_T2B, UINT32, "NN_P10tcl_T2B", 10000); \
        EXT_STR_ITEM_INFO_ZZP(                                                                                    \
            ok, si, offset, struct_t, printerr, NN_P10tcl_T2Bv, UINT32, "NN_P10tcl_T2Bv", "NN_P10tcl_T2B");       \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, NN_P10tfl_T1BM, UINT32, "NN_P10tfl_T1BM", 50);  \
        EXT_STR_ITEM_INFO_ZZP(                                                                                    \
            ok, si, offset, struct_t, printerr, NN_P10tfl_T1BMI, UINT32, "NN_P10tfl_T1BMI", "NN_P10tfl_T1BM");    \
        EXT_STR_ITEM_INFO_ZZP(                                                                                    \
            ok, si, offset, struct_t, printerr, NN_P10tfl_T1BME, UINT32, "NN_P10tfl_T1BME", "NN_P10tfl_T1BM");    \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, NN_P10tfl_T1B, UINT32, "NN_P10tfl_T1B", 10000); \
        EXT_STR_ITEM_INFO_ZZP(                                                                                    \
            ok, si, offset, struct_t, printerr, NN_P10tfl_T1Bv, UINT32, "NN_P10tfl_T1Bv", "NN_P10tfl_T1B");       \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, NN_P10tfl_T2BM, UINT32, "NN_P10tfl_T2BM", 50);  \
        EXT_STR_ITEM_INFO_ZZP(                                                                                    \
            ok, si, offset, struct_t, printerr, NN_P10tfl_T2BMI, UINT32, "NN_P10tfl_T2BMI", "NN_P10tfl_T2BM");    \
        EXT_STR_ITEM_INFO_ZZP(                                                                                    \
            ok, si, offset, struct_t, printerr, NN_P10tfl_T2BME, UINT32, "NN_P10tfl_T2BME", "NN_P10tfl_T2BM");    \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, NN_P10tfl_T2B, UINT32, "NN_P10tfl_T2B", 10000); \
        EXT_STR_ITEM_INFO_ZZP(                                                                                    \
            ok, si, offset, struct_t, printerr, NN_P10tfl_T2Bv, UINT32, "NN_P10tfl_T2Bv", "NN_P10tfl_T2B");       \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, NN_P10tct_T1BM, UINT32, "NN_P10tct_T1BM", 50);  \
        EXT_STR_ITEM_INFO_ZZP(                                                                                    \
            ok, si, offset, struct_t, printerr, NN_P10tct_T1BMI, UINT32, "NN_P10tct_T1BMI", "NN_P10tct_T1BM");    \
        EXT_STR_ITEM_INFO_ZZP(                                                                                    \
            ok, si, offset, struct_t, printerr, NN_P10tct_T1BME, UINT32, "NN_P10tct_T1BME", "NN_P10tct_T1BM");    \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, NN_P10tct_T1B, UINT32, "NN_P10tct_T1B", 10000); \
        EXT_STR_ITEM_INFO_ZZP(                                                                                    \
            ok, si, offset, struct_t, printerr, NN_P10tct_T1Bv, UINT32, "NN_P10tct_T1Bv", "NN_P10tct_T1B");       \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, NN_P10tct_T2BM, UINT32, "NN_P10tct_T2BM", 50);  \
        EXT_STR_ITEM_INFO_ZZP(                                                                                    \
            ok, si, offset, struct_t, printerr, NN_P10tct_T2BMI, UINT32, "NN_P10tct_T2BMI", "NN_P10tct_T2BM");    \
        EXT_STR_ITEM_INFO_ZZP(                                                                                    \
            ok, si, offset, struct_t, printerr, NN_P10tct_T2BME, UINT32, "NN_P10tct_T2BME", "NN_P10tct_T2BM");    \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, NN_P10tct_T2B, UINT32, "NN_P10tct_T2B", 10000); \
        EXT_STR_ITEM_INFO_ZZP(                                                                                    \
            ok, si, offset, struct_t, printerr, NN_P10tct_T2Bv, UINT32, "NN_P10tct_T2Bv", "NN_P10tct_T2B");       \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, NN_P10tft_T1BM, UINT32, "NN_P10tft_T1BM", 50);  \
        EXT_STR_ITEM_INFO_ZZP(                                                                                    \
            ok, si, offset, struct_t, printerr, NN_P10tft_T1BMI, UINT32, "NN_P10tft_T1BMI", "NN_P10tft_T1BM");    \
        EXT_STR_ITEM_INFO_ZZP(                                                                                    \
            ok, si, offset, struct_t, printerr, NN_P10tft_T1BME, UINT32, "NN_P10tft_T1BME", "NN_P10tft_T1BM");    \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, NN_P10tft_T1B, UINT32, "NN_P10tft_T1B", 10000); \
        EXT_STR_ITEM_INFO_ZZP(                                                                                    \
            ok, si, offset, struct_t, printerr, NN_P10tft_T1Bv, UINT32, "NN_P10tft_T1Bv", "NN_P10tft_T1B");       \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, NN_P10tft_T2BM, UINT32, "NN_P10tft_T2BM", 50);  \
        EXT_STR_ITEM_INFO_ZZP(                                                                                    \
            ok, si, offset, struct_t, printerr, NN_P10tft_T2BMI, UINT32, "NN_P10tft_T2BMI", "NN_P10tft_T2BM");    \
        EXT_STR_ITEM_INFO_ZZP(                                                                                    \
            ok, si, offset, struct_t, printerr, NN_P10tft_T2BME, UINT32, "NN_P10tft_T2BME", "NN_P10tft_T2BM");    \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, NN_P10tft_T2B, UINT32, "NN_P10tft_T2B", 10000); \
        EXT_STR_ITEM_INFO_ZZP(                                                                                    \
            ok, si, offset, struct_t, printerr, NN_P10tft_T2Bv, UINT32, "NN_P10tft_T2Bv", "NN_P10tft_T2B");       \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, NN_P11tcl_T1BM, UINT32, "NN_P11tcl_T1BM", 50);  \
        EXT_STR_ITEM_INFO_ZZP(                                                                                    \
            ok, si, offset, struct_t, printerr, NN_P11tcl_T1BMI, UINT32, "NN_P11tcl_T1BMI", "NN_P11tcl_T1BM");    \
        EXT_STR_ITEM_INFO_ZZP(                                                                                    \
            ok, si, offset, struct_t, printerr, NN_P11tcl_T1BME, UINT32, "NN_P11tcl_T1BME", "NN_P11tcl_T1BM");    \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, NN_P11tcl_T1B, UINT32, "NN_P11tcl_T1B", 10000); \
        EXT_STR_ITEM_INFO_ZZP(                                                                                    \
            ok, si, offset, struct_t, printerr, NN_P11tcl_T1Bv, UINT32, "NN_P11tcl_T1Bv", "NN_P11tcl_T1B");       \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, NN_P11tcl_T2BM, UINT32, "NN_P11tcl_T2BM", 50);  \
        EXT_STR_ITEM_INFO_ZZP(                                                                                    \
            ok, si, offset, struct_t, printerr, NN_P11tcl_T2BMI, UINT32, "NN_P11tcl_T2BMI", "NN_P11tcl_T2BM");    \
        EXT_STR_ITEM_INFO_ZZP(                                                                                    \
            ok, si, offset, struct_t, printerr, NN_P11tcl_T2BME, UINT32, "NN_P11tcl_T2BME", "NN_P11tcl_T2BM");    \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, NN_P11tcl_T2B, UINT32, "NN_P11tcl_T2B", 10000); \
        EXT_STR_ITEM_INFO_ZZP(                                                                                    \
            ok, si, offset, struct_t, printerr, NN_P11tcl_T2Bv, UINT32, "NN_P11tcl_T2Bv", "NN_P11tcl_T2B");       \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, NN_P11tfl_T1BM, UINT32, "NN_P11tfl_T1BM", 50);  \
        EXT_STR_ITEM_INFO_ZZP(                                                                                    \
            ok, si, offset, struct_t, printerr, NN_P11tfl_T1BMI, UINT32, "NN_P11tfl_T1BMI", "NN_P11tfl_T1BM");    \
        EXT_STR_ITEM_INFO_ZZP(                                                                                    \
            ok, si, offset, struct_t, printerr, NN_P11tfl_T1BME, UINT32, "NN_P11tfl_T1BME", "NN_P11tfl_T1BM");    \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, NN_P11tfl_T1B, UINT32, "NN_P11tfl_T1B", 10000); \
        EXT_STR_ITEM_INFO_ZZP(                                                                                    \
            ok, si, offset, struct_t, printerr, NN_P11tfl_T1Bv, UINT32, "NN_P11tfl_T1Bv", "NN_P11tfl_T1B");       \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, NN_P11tfl_T2BM, UINT32, "NN_P11tfl_T2BM", 50);  \
        EXT_STR_ITEM_INFO_ZZP(                                                                                    \
            ok, si, offset, struct_t, printerr, NN_P11tfl_T2BMI, UINT32, "NN_P11tfl_T2BMI", "NN_P11tfl_T2BM");    \
        EXT_STR_ITEM_INFO_ZZP(                                                                                    \
            ok, si, offset, struct_t, printerr, NN_P11tfl_T2BME, UINT32, "NN_P11tfl_T2BME", "NN_P11tfl_T2BM");    \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, NN_P11tfl_T2B, UINT32, "NN_P11tfl_T2B", 10000); \
        EXT_STR_ITEM_INFO_ZZP(                                                                                    \
            ok, si, offset, struct_t, printerr, NN_P11tfl_T2Bv, UINT32, "NN_P11tfl_T2Bv", "NN_P11tfl_T2B");       \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, NN_P11tct_T1BM, UINT32, "NN_P11tct_T1BM", 50);  \
        EXT_STR_ITEM_INFO_ZZP(                                                                                    \
            ok, si, offset, struct_t, printerr, NN_P11tct_T1BMI, UINT32, "NN_P11tct_T1BMI", "NN_P11tct_T1BM");    \
        EXT_STR_ITEM_INFO_ZZP(                                                                                    \
            ok, si, offset, struct_t, printerr, NN_P11tct_T1BME, UINT32, "NN_P11tct_T1BME", "NN_P11tct_T1BM");    \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, NN_P11tct_T1B, UINT32, "NN_P11tct_T1B", 10000); \
        EXT_STR_ITEM_INFO_ZZP(                                                                                    \
            ok, si, offset, struct_t, printerr, NN_P11tct_T1Bv, UINT32, "NN_P11tct_T1Bv", "NN_P11tct_T1B");       \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, NN_P11tct_T2BM, UINT32, "NN_P11tct_T2BM", 50);  \
        EXT_STR_ITEM_INFO_ZZP(                                                                                    \
            ok, si, offset, struct_t, printerr, NN_P11tct_T2BMI, UINT32, "NN_P11tct_T2BMI", "NN_P11tct_T2BM");    \
        EXT_STR_ITEM_INFO_ZZP(                                                                                    \
            ok, si, offset, struct_t, printerr, NN_P11tct_T2BME, UINT32, "NN_P11tct_T2BME", "NN_P11tct_T2BM");    \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, NN_P11tct_T2B, UINT32, "NN_P11tct_T2B", 10000); \
        EXT_STR_ITEM_INFO_ZZP(                                                                                    \
            ok, si, offset, struct_t, printerr, NN_P11tct_T2Bv, UINT32, "NN_P11tct_T2Bv", "NN_P11tct_T2B");       \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, NN_P11tft_T1BM, UINT32, "NN_P11tft_T1BM", 50);  \
        EXT_STR_ITEM_INFO_ZZP(                                                                                    \
            ok, si, offset, struct_t, printerr, NN_P11tft_T1BMI, UINT32, "NN_P11tft_T1BMI", "NN_P11tft_T1BM");    \
        EXT_STR_ITEM_INFO_ZZP(                                                                                    \
            ok, si, offset, struct_t, printerr, NN_P11tft_T1BME, UINT32, "NN_P11tft_T1BME", "NN_P11tft_T1BM");    \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, NN_P11tft_T1B, UINT32, "NN_P11tft_T1B", 10000); \
        EXT_STR_ITEM_INFO_ZZP(                                                                                    \
            ok, si, offset, struct_t, printerr, NN_P11tft_T1Bv, UINT32, "NN_P11tft_T1Bv", "NN_P11tft_T1B");       \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, NN_P11tft_T2BM, UINT32, "NN_P11tft_T2BM", 50);  \
        EXT_STR_ITEM_INFO_ZZP(                                                                                    \
            ok, si, offset, struct_t, printerr, NN_P11tft_T2BMI, UINT32, "NN_P11tft_T2BMI", "NN_P11tft_T2BM");    \
        EXT_STR_ITEM_INFO_ZZP(                                                                                    \
            ok, si, offset, struct_t, printerr, NN_P11tft_T2BME, UINT32, "NN_P11tft_T2BME", "NN_P11tft_T2BM");    \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, NN_P11tft_T2B, UINT32, "NN_P11tft_T2B", 10000); \
        EXT_STR_ITEM_INFO_ZZP(                                                                                    \
            ok, si, offset, struct_t, printerr, NN_P11tft_T2Bv, UINT32, "NN_P11tft_T2Bv", "NN_P11tft_T2B");       \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, NN_P12tcl_T1BM, UINT32, "NN_P12tcl_T1BM", 50);  \
        EXT_STR_ITEM_INFO_ZZP(                                                                                    \
            ok, si, offset, struct_t, printerr, NN_P12tcl_T1BMI, UINT32, "NN_P12tcl_T1BMI", "NN_P12tcl_T1BM");    \
        EXT_STR_ITEM_INFO_ZZP(                                                                                    \
            ok, si, offset, struct_t, printerr, NN_P12tcl_T1BME, UINT32, "NN_P12tcl_T1BME", "NN_P12tcl_T1BM");    \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, NN_P12tcl_T1B, UINT32, "NN_P12tcl_T1B", 10000); \
        EXT_STR_ITEM_INFO_ZZP(                                                                                    \
            ok, si, offset, struct_t, printerr, NN_P12tcl_T1Bv, UINT32, "NN_P12tcl_T1Bv", "NN_P12tcl_T1B");       \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, NN_P12tcl_T2BM, UINT32, "NN_P12tcl_T2BM", 50);  \
        EXT_STR_ITEM_INFO_ZZP(                                                                                    \
            ok, si, offset, struct_t, printerr, NN_P12tcl_T2BMI, UINT32, "NN_P12tcl_T2BMI", "NN_P12tcl_T2BM");    \
        EXT_STR_ITEM_INFO_ZZP(                                                                                    \
            ok, si, offset, struct_t, printerr, NN_P12tcl_T2BME, UINT32, "NN_P12tcl_T2BME", "NN_P12tcl_T2BM");    \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, NN_P12tcl_T2B, UINT32, "NN_P12tcl_T2B", 10000); \
        EXT_STR_ITEM_INFO_ZZP(                                                                                    \
            ok, si, offset, struct_t, printerr, NN_P12tcl_T2Bv, UINT32, "NN_P12tcl_T2Bv", "NN_P12tcl_T2B");       \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, NN_P12tfl_T1BM, UINT32, "NN_P12tfl_T1BM", 50);  \
        EXT_STR_ITEM_INFO_ZZP(                                                                                    \
            ok, si, offset, struct_t, printerr, NN_P12tfl_T1BMI, UINT32, "NN_P12tfl_T1BMI", "NN_P12tfl_T1BM");    \
        EXT_STR_ITEM_INFO_ZZP(                                                                                    \
            ok, si, offset, struct_t, printerr, NN_P12tfl_T1BME, UINT32, "NN_P12tfl_T1BME", "NN_P12tfl_T1BM");    \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, NN_P12tfl_T1B, UINT32, "NN_P12tfl_T1B", 10000); \
        EXT_STR_ITEM_INFO_ZZP(                                                                                    \
            ok, si, offset, struct_t, printerr, NN_P12tfl_T1Bv, UINT32, "NN_P12tfl_T1Bv", "NN_P12tfl_T1B");       \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, NN_P12tfl_T2BM, UINT32, "NN_P12tfl_T2BM", 50);  \
        EXT_STR_ITEM_INFO_ZZP(                                                                                    \
            ok, si, offset, struct_t, printerr, NN_P12tfl_T2BMI, UINT32, "NN_P12tfl_T2BMI", "NN_P12tfl_T2BM");    \
        EXT_STR_ITEM_INFO_ZZP(                                                                                    \
            ok, si, offset, struct_t, printerr, NN_P12tfl_T2BME, UINT32, "NN_P12tfl_T2BME", "NN_P12tfl_T2BM");    \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, NN_P12tfl_T2B, UINT32, "NN_P12tfl_T2B", 10000); \
        EXT_STR_ITEM_INFO_ZZP(                                                                                    \
            ok, si, offset, struct_t, printerr, NN_P12tfl_T2Bv, UINT32, "NN_P12tfl_T2Bv", "NN_P12tfl_T2B");       \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, NN_P12tct_T1BM, UINT32, "NN_P12tct_T1BM", 50);  \
        EXT_STR_ITEM_INFO_ZZP(                                                                                    \
            ok, si, offset, struct_t, printerr, NN_P12tct_T1BMI, UINT32, "NN_P12tct_T1BMI", "NN_P12tct_T1BM");    \
        EXT_STR_ITEM_INFO_ZZP(                                                                                    \
            ok, si, offset, struct_t, printerr, NN_P12tct_T1BME, UINT32, "NN_P12tct_T1BME", "NN_P12tct_T1BM");    \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, NN_P12tct_T1B, UINT32, "NN_P12tct_T1B", 10000); \
        EXT_STR_ITEM_INFO_ZZP(                                                                                    \
            ok, si, offset, struct_t, printerr, NN_P12tct_T1Bv, UINT32, "NN_P12tct_T1Bv", "NN_P12tct_T1B");       \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, NN_P12tct_T2BM, UINT32, "NN_P12tct_T2BM", 50);  \
        EXT_STR_ITEM_INFO_ZZP(                                                                                    \
            ok, si, offset, struct_t, printerr, NN_P12tct_T2BMI, UINT32, "NN_P12tct_T2BMI", "NN_P12tct_T2BM");    \
        EXT_STR_ITEM_INFO_ZZP(                                                                                    \
            ok, si, offset, struct_t, printerr, NN_P12tct_T2BME, UINT32, "NN_P12tct_T2BME", "NN_P12tct_T2BM");    \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, NN_P12tct_T2B, UINT32, "NN_P12tct_T2B", 10000); \
        EXT_STR_ITEM_INFO_ZZP(                                                                                    \
            ok, si, offset, struct_t, printerr, NN_P12tct_T2Bv, UINT32, "NN_P12tct_T2Bv", "NN_P12tct_T2B");       \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, NN_P12tft_T1BM, UINT32, "NN_P12tft_T1BM", 50);  \
        EXT_STR_ITEM_INFO_ZZP(                                                                                    \
            ok, si, offset, struct_t, printerr, NN_P12tft_T1BMI, UINT32, "NN_P12tft_T1BMI", "NN_P12tft_T1BM");    \
        EXT_STR_ITEM_INFO_ZZP(                                                                                    \
            ok, si, offset, struct_t, printerr, NN_P12tft_T1BME, UINT32, "NN_P12tft_T1BME", "NN_P12tft_T1BM");    \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, NN_P12tft_T1B, UINT32, "NN_P12tft_T1B", 10000); \
        EXT_STR_ITEM_INFO_ZZP(                                                                                    \
            ok, si, offset, struct_t, printerr, NN_P12tft_T1Bv, UINT32, "NN_P12tft_T1Bv", "NN_P12tft_T1B");       \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, NN_P12tft_T2BM, UINT32, "NN_P12tft_T2BM", 50);  \
        EXT_STR_ITEM_INFO_ZZP(                                                                                    \
            ok, si, offset, struct_t, printerr, NN_P12tft_T2BMI, UINT32, "NN_P12tft_T2BMI", "NN_P12tft_T2BM");    \
        EXT_STR_ITEM_INFO_ZZP(                                                                                    \
            ok, si, offset, struct_t, printerr, NN_P12tft_T2BME, UINT32, "NN_P12tft_T2BME", "NN_P12tft_T2BM");    \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, NN_P12tft_T2B, UINT32, "NN_P12tft_T2B", 10000); \
        EXT_STR_ITEM_INFO_ZZP(                                                                                    \
            ok, si, offset, struct_t, printerr, NN_P12tft_T2Bv, UINT32, "NN_P12tft_T2Bv", "NN_P12tft_T2B");       \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, NN_P13tcl_T1BM, UINT32, "NN_P13tcl_T1BM", 50);  \
        EXT_STR_ITEM_INFO_ZZP(                                                                                    \
            ok, si, offset, struct_t, printerr, NN_P13tcl_T1BMI, UINT32, "NN_P13tcl_T1BMI", "NN_P13tcl_T1BM");    \
        EXT_STR_ITEM_INFO_ZZP(                                                                                    \
            ok, si, offset, struct_t, printerr, NN_P13tcl_T1BME, UINT32, "NN_P13tcl_T1BME", "NN_P13tcl_T1BM");    \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, NN_P13tcl_T1B, UINT32, "NN_P13tcl_T1B", 10000); \
        EXT_STR_ITEM_INFO_ZZP(                                                                                    \
            ok, si, offset, struct_t, printerr, NN_P13tcl_T1Bv, UINT32, "NN_P13tcl_T1Bv", "NN_P13tcl_T1B");       \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, NN_P13tcl_T2BM, UINT32, "NN_P13tcl_T2BM", 50);  \
        EXT_STR_ITEM_INFO_ZZP(                                                                                    \
            ok, si, offset, struct_t, printerr, NN_P13tcl_T2BMI, UINT32, "NN_P13tcl_T2BMI", "NN_P13tcl_T2BM");    \
        EXT_STR_ITEM_INFO_ZZP(                                                                                    \
            ok, si, offset, struct_t, printerr, NN_P13tcl_T2BME, UINT32, "NN_P13tcl_T2BME", "NN_P13tcl_T2BM");    \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, NN_P13tcl_T2B, UINT32, "NN_P13tcl_T2B", 10000); \
        EXT_STR_ITEM_INFO_ZZP(                                                                                    \
            ok, si, offset, struct_t, printerr, NN_P13tcl_T2Bv, UINT32, "NN_P13tcl_T2Bv", "NN_P13tcl_T2B");       \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, NN_P13tfl_T1BM, UINT32, "NN_P13tfl_T1BM", 50);  \
        EXT_STR_ITEM_INFO_ZZP(                                                                                    \
            ok, si, offset, struct_t, printerr, NN_P13tfl_T1BMI, UINT32, "NN_P13tfl_T1BMI", "NN_P13tfl_T1BM");    \
        EXT_STR_ITEM_INFO_ZZP(                                                                                    \
            ok, si, offset, struct_t, printerr, NN_P13tfl_T1BME, UINT32, "NN_P13tfl_T1BME", "NN_P13tfl_T1BM");    \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, NN_P13tfl_T1B, UINT32, "NN_P13tfl_T1B", 10000); \
        EXT_STR_ITEM_INFO_ZZP(                                                                                    \
            ok, si, offset, struct_t, printerr, NN_P13tfl_T1Bv, UINT32, "NN_P13tfl_T1Bv", "NN_P13tfl_T1B");       \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, NN_P13tfl_T2BM, UINT32, "NN_P13tfl_T2BM", 50);  \
        EXT_STR_ITEM_INFO_ZZP(                                                                                    \
            ok, si, offset, struct_t, printerr, NN_P13tfl_T2BMI, UINT32, "NN_P13tfl_T2BMI", "NN_P13tfl_T2BM");    \
        EXT_STR_ITEM_INFO_ZZP(                                                                                    \
            ok, si, offset, struct_t, printerr, NN_P13tfl_T2BME, UINT32, "NN_P13tfl_T2BME", "NN_P13tfl_T2BM");    \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, NN_P13tfl_T2B, UINT32, "NN_P13tfl_T2B", 10000); \
        EXT_STR_ITEM_INFO_ZZP(                                                                                    \
            ok, si, offset, struct_t, printerr, NN_P13tfl_T2Bv, UINT32, "NN_P13tfl_T2Bv", "NN_P13tfl_T2B");       \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, NN_P13tct_T1BM, UINT32, "NN_P13tct_T1BM", 50);  \
        EXT_STR_ITEM_INFO_ZZP(                                                                                    \
            ok, si, offset, struct_t, printerr, NN_P13tct_T1BMI, UINT32, "NN_P13tct_T1BMI", "NN_P13tct_T1BM");    \
        EXT_STR_ITEM_INFO_ZZP(                                                                                    \
            ok, si, offset, struct_t, printerr, NN_P13tct_T1BME, UINT32, "NN_P13tct_T1BME", "NN_P13tct_T1BM");    \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, NN_P13tct_T1B, UINT32, "NN_P13tct_T1B", 10000); \
        EXT_STR_ITEM_INFO_ZZP(                                                                                    \
            ok, si, offset, struct_t, printerr, NN_P13tct_T1Bv, UINT32, "NN_P13tct_T1Bv", "NN_P13tct_T1B");       \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, NN_P13tct_T2BM, UINT32, "NN_P13tct_T2BM", 50);  \
        EXT_STR_ITEM_INFO_ZZP(                                                                                    \
            ok, si, offset, struct_t, printerr, NN_P13tct_T2BMI, UINT32, "NN_P13tct_T2BMI", "NN_P13tct_T2BM");    \
        EXT_STR_ITEM_INFO_ZZP(                                                                                    \
            ok, si, offset, struct_t, printerr, NN_P13tct_T2BME, UINT32, "NN_P13tct_T2BME", "NN_P13tct_T2BM");    \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, NN_P13tct_T2B, UINT32, "NN_P13tct_T2B", 10000); \
        EXT_STR_ITEM_INFO_ZZP(                                                                                    \
            ok, si, offset, struct_t, printerr, NN_P13tct_T2Bv, UINT32, "NN_P13tct_T2Bv", "NN_P13tct_T2B");       \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, NN_P13tft_T1BM, UINT32, "NN_P13tft_T1BM", 50);  \
        EXT_STR_ITEM_INFO_ZZP(                                                                                    \
            ok, si, offset, struct_t, printerr, NN_P13tft_T1BMI, UINT32, "NN_P13tft_T1BMI", "NN_P13tft_T1BM");    \
        EXT_STR_ITEM_INFO_ZZP(                                                                                    \
            ok, si, offset, struct_t, printerr, NN_P13tft_T1BME, UINT32, "NN_P13tft_T1BME", "NN_P13tft_T1BM");    \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, NN_P13tft_T1B, UINT32, "NN_P13tft_T1B", 10000); \
        EXT_STR_ITEM_INFO_ZZP(                                                                                    \
            ok, si, offset, struct_t, printerr, NN_P13tft_T1Bv, UINT32, "NN_P13tft_T1Bv", "NN_P13tft_T1B");       \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, NN_P13tft_T2BM, UINT32, "NN_P13tft_T2BM", 50);  \
        EXT_STR_ITEM_INFO_ZZP(                                                                                    \
            ok, si, offset, struct_t, printerr, NN_P13tft_T2BMI, UINT32, "NN_P13tft_T2BMI", "NN_P13tft_T2BM");    \
        EXT_STR_ITEM_INFO_ZZP(                                                                                    \
            ok, si, offset, struct_t, printerr, NN_P13tft_T2BME, UINT32, "NN_P13tft_T2BME", "NN_P13tft_T2BM");    \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, NN_P13tft_T2B, UINT32, "NN_P13tft_T2B", 10000); \
        EXT_STR_ITEM_INFO_ZZP(                                                                                    \
            ok, si, offset, struct_t, printerr, NN_P13tft_T2Bv, UINT32, "NN_P13tft_T2Bv", "NN_P13tft_T2B");       \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, NN_P14tcl_T1BM, UINT32, "NN_P14tcl_T1BM", 50);  \
        EXT_STR_ITEM_INFO_ZZP(                                                                                    \
            ok, si, offset, struct_t, printerr, NN_P14tcl_T1BMI, UINT32, "NN_P14tcl_T1BMI", "NN_P14tcl_T1BM");    \
        EXT_STR_ITEM_INFO_ZZP(                                                                                    \
            ok, si, offset, struct_t, printerr, NN_P14tcl_T1BME, UINT32, "NN_P14tcl_T1BME", "NN_P14tcl_T1BM");    \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, NN_P14tcl_T1B, UINT32, "NN_P14tcl_T1B", 10000); \
        EXT_STR_ITEM_INFO_ZZP(                                                                                    \
            ok, si, offset, struct_t, printerr, NN_P14tcl_T1Bv, UINT32, "NN_P14tcl_T1Bv", "NN_P14tcl_T1B");       \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, NN_P14tcl_T2BM, UINT32, "NN_P14tcl_T2BM", 50);  \
        EXT_STR_ITEM_INFO_ZZP(                                                                                    \
            ok, si, offset, struct_t, printerr, NN_P14tcl_T2BMI, UINT32, "NN_P14tcl_T2BMI", "NN_P14tcl_T2BM");    \
        EXT_STR_ITEM_INFO_ZZP(                                                                                    \
            ok, si, offset, struct_t, printerr, NN_P14tcl_T2BME, UINT32, "NN_P14tcl_T2BME", "NN_P14tcl_T2BM");    \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, NN_P14tcl_T2B, UINT32, "NN_P14tcl_T2B", 10000); \
        EXT_STR_ITEM_INFO_ZZP(                                                                                    \
            ok, si, offset, struct_t, printerr, NN_P14tcl_T2Bv, UINT32, "NN_P14tcl_T2Bv", "NN_P14tcl_T2B");       \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, NN_P14tfl_T1BM, UINT32, "NN_P14tfl_T1BM", 50);  \
        EXT_STR_ITEM_INFO_ZZP(                                                                                    \
            ok, si, offset, struct_t, printerr, NN_P14tfl_T1BMI, UINT32, "NN_P14tfl_T1BMI", "NN_P14tfl_T1BM");    \
        EXT_STR_ITEM_INFO_ZZP(                                                                                    \
            ok, si, offset, struct_t, printerr, NN_P14tfl_T1BME, UINT32, "NN_P14tfl_T1BME", "NN_P14tfl_T1BM");    \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, NN_P14tfl_T1B, UINT32, "NN_P14tfl_T1B", 10000); \
        EXT_STR_ITEM_INFO_ZZP(                                                                                    \
            ok, si, offset, struct_t, printerr, NN_P14tfl_T1Bv, UINT32, "NN_P14tfl_T1Bv", "NN_P14tfl_T1B");       \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, NN_P14tfl_T2BM, UINT32, "NN_P14tfl_T2BM", 50);  \
        EXT_STR_ITEM_INFO_ZZP(                                                                                    \
            ok, si, offset, struct_t, printerr, NN_P14tfl_T2BMI, UINT32, "NN_P14tfl_T2BMI", "NN_P14tfl_T2BM");    \
        EXT_STR_ITEM_INFO_ZZP(                                                                                    \
            ok, si, offset, struct_t, printerr, NN_P14tfl_T2BME, UINT32, "NN_P14tfl_T2BME", "NN_P14tfl_T2BM");    \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, NN_P14tfl_T2B, UINT32, "NN_P14tfl_T2B", 10000); \
        EXT_STR_ITEM_INFO_ZZP(                                                                                    \
            ok, si, offset, struct_t, printerr, NN_P14tfl_T2Bv, UINT32, "NN_P14tfl_T2Bv", "NN_P14tfl_T2B");       \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, NN_P14tct_T1BM, UINT32, "NN_P14tct_T1BM", 50);  \
        EXT_STR_ITEM_INFO_ZZP(                                                                                    \
            ok, si, offset, struct_t, printerr, NN_P14tct_T1BMI, UINT32, "NN_P14tct_T1BMI", "NN_P14tct_T1BM");    \
        EXT_STR_ITEM_INFO_ZZP(                                                                                    \
            ok, si, offset, struct_t, printerr, NN_P14tct_T1BME, UINT32, "NN_P14tct_T1BME", "NN_P14tct_T1BM");    \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, NN_P14tct_T1B, UINT32, "NN_P14tct_T1B", 10000); \
        EXT_STR_ITEM_INFO_ZZP(                                                                                    \
            ok, si, offset, struct_t, printerr, NN_P14tct_T1Bv, UINT32, "NN_P14tct_T1Bv", "NN_P14tct_T1B");       \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, NN_P14tct_T2BM, UINT32, "NN_P14tct_T2BM", 50);  \
        EXT_STR_ITEM_INFO_ZZP(                                                                                    \
            ok, si, offset, struct_t, printerr, NN_P14tct_T2BMI, UINT32, "NN_P14tct_T2BMI", "NN_P14tct_T2BM");    \
        EXT_STR_ITEM_INFO_ZZP(                                                                                    \
            ok, si, offset, struct_t, printerr, NN_P14tct_T2BME, UINT32, "NN_P14tct_T2BME", "NN_P14tct_T2BM");    \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, NN_P14tct_T2B, UINT32, "NN_P14tct_T2B", 10000); \
        EXT_STR_ITEM_INFO_ZZP(                                                                                    \
            ok, si, offset, struct_t, printerr, NN_P14tct_T2Bv, UINT32, "NN_P14tct_T2Bv", "NN_P14tct_T2B");       \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, NN_P14tft_T1BM, UINT32, "NN_P14tft_T1BM", 50);  \
        EXT_STR_ITEM_INFO_ZZP(                                                                                    \
            ok, si, offset, struct_t, printerr, NN_P14tft_T1BMI, UINT32, "NN_P14tft_T1BMI", "NN_P14tft_T1BM");    \
        EXT_STR_ITEM_INFO_ZZP(                                                                                    \
            ok, si, offset, struct_t, printerr, NN_P14tft_T1BME, UINT32, "NN_P14tft_T1BME", "NN_P14tft_T1BM");    \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, NN_P14tft_T1B, UINT32, "NN_P14tft_T1B", 10000); \
        EXT_STR_ITEM_INFO_ZZP(                                                                                    \
            ok, si, offset, struct_t, printerr, NN_P14tft_T1Bv, UINT32, "NN_P14tft_T1Bv", "NN_P14tft_T1B");       \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, NN_P14tft_T2BM, UINT32, "NN_P14tft_T2BM", 50);  \
        EXT_STR_ITEM_INFO_ZZP(                                                                                    \
            ok, si, offset, struct_t, printerr, NN_P14tft_T2BMI, UINT32, "NN_P14tft_T2BMI", "NN_P14tft_T2BM");    \
        EXT_STR_ITEM_INFO_ZZP(                                                                                    \
            ok, si, offset, struct_t, printerr, NN_P14tft_T2BME, UINT32, "NN_P14tft_T2BME", "NN_P14tft_T2BM");    \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, NN_P14tft_T2B, UINT32, "NN_P14tft_T2B", 10000); \
        EXT_STR_ITEM_INFO_ZZP(                                                                                    \
            ok, si, offset, struct_t, printerr, NN_P14tft_T2Bv, UINT32, "NN_P14tft_T2Bv", "NN_P14tft_T2B");       \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, NN_P15tcl_T1BM, UINT32, "NN_P15tcl_T1BM", 50);  \
        EXT_STR_ITEM_INFO_ZZP(                                                                                    \
            ok, si, offset, struct_t, printerr, NN_P15tcl_T1BMI, UINT32, "NN_P15tcl_T1BMI", "NN_P15tcl_T1BM");    \
        EXT_STR_ITEM_INFO_ZZP(                                                                                    \
            ok, si, offset, struct_t, printerr, NN_P15tcl_T1BME, UINT32, "NN_P15tcl_T1BME", "NN_P15tcl_T1BM");    \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, NN_P15tcl_T1B, UINT32, "NN_P15tcl_T1B", 10000); \
        EXT_STR_ITEM_INFO_ZZP(                                                                                    \
            ok, si, offset, struct_t, printerr, NN_P15tcl_T1Bv, UINT32, "NN_P15tcl_T1Bv", "NN_P15tcl_T1B");       \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, NN_P15tcl_T2BM, UINT32, "NN_P15tcl_T2BM", 50);  \
        EXT_STR_ITEM_INFO_ZZP(                                                                                    \
            ok, si, offset, struct_t, printerr, NN_P15tcl_T2BMI, UINT32, "NN_P15tcl_T2BMI", "NN_P15tcl_T2BM");    \
        EXT_STR_ITEM_INFO_ZZP(                                                                                    \
            ok, si, offset, struct_t, printerr, NN_P15tcl_T2BME, UINT32, "NN_P15tcl_T2BME", "NN_P15tcl_T2BM");    \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, NN_P15tcl_T2B, UINT32, "NN_P15tcl_T2B", 10000); \
        EXT_STR_ITEM_INFO_ZZP(                                                                                    \
            ok, si, offset, struct_t, printerr, NN_P15tcl_T2Bv, UINT32, "NN_P15tcl_T2Bv", "NN_P15tcl_T2B");       \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, NN_P15tfl_T1BM, UINT32, "NN_P15tfl_T1BM", 50);  \
        EXT_STR_ITEM_INFO_ZZP(                                                                                    \
            ok, si, offset, struct_t, printerr, NN_P15tfl_T1BMI, UINT32, "NN_P15tfl_T1BMI", "NN_P15tfl_T1BM");    \
        EXT_STR_ITEM_INFO_ZZP(                                                                                    \
            ok, si, offset, struct_t, printerr, NN_P15tfl_T1BME, UINT32, "NN_P15tfl_T1BME", "NN_P15tfl_T1BM");    \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, NN_P15tfl_T1B, UINT32, "NN_P15tfl_T1B", 10000); \
        EXT_STR_ITEM_INFO_ZZP(                                                                                    \
            ok, si, offset, struct_t, printerr, NN_P15tfl_T1Bv, UINT32, "NN_P15tfl_T1Bv", "NN_P15tfl_T1B");       \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, NN_P15tfl_T2BM, UINT32, "NN_P15tfl_T2BM", 50);  \
        EXT_STR_ITEM_INFO_ZZP(                                                                                    \
            ok, si, offset, struct_t, printerr, NN_P15tfl_T2BMI, UINT32, "NN_P15tfl_T2BMI", "NN_P15tfl_T2BM");    \
        EXT_STR_ITEM_INFO_ZZP(                                                                                    \
            ok, si, offset, struct_t, printerr, NN_P15tfl_T2BME, UINT32, "NN_P15tfl_T2BME", "NN_P15tfl_T2BM");    \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, NN_P15tfl_T2B, UINT32, "NN_P15tfl_T2B", 10000); \
        EXT_STR_ITEM_INFO_ZZP(                                                                                    \
            ok, si, offset, struct_t, printerr, NN_P15tfl_T2Bv, UINT32, "NN_P15tfl_T2Bv", "NN_P15tfl_T2B");       \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, NN_P15tct_T1BM, UINT32, "NN_P15tct_T1BM", 50);  \
        EXT_STR_ITEM_INFO_ZZP(                                                                                    \
            ok, si, offset, struct_t, printerr, NN_P15tct_T1BMI, UINT32, "NN_P15tct_T1BMI", "NN_P15tct_T1BM");    \
        EXT_STR_ITEM_INFO_ZZP(                                                                                    \
            ok, si, offset, struct_t, printerr, NN_P15tct_T1BME, UINT32, "NN_P15tct_T1BME", "NN_P15tct_T1BM");    \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, NN_P15tct_T1B, UINT32, "NN_P15tct_T1B", 10000); \
        EXT_STR_ITEM_INFO_ZZP(                                                                                    \
            ok, si, offset, struct_t, printerr, NN_P15tct_T1Bv, UINT32, "NN_P15tct_T1Bv", "NN_P15tct_T1B");       \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, NN_P15tct_T2BM, UINT32, "NN_P15tct_T2BM", 50);  \
        EXT_STR_ITEM_INFO_ZZP(                                                                                    \
            ok, si, offset, struct_t, printerr, NN_P15tct_T2BMI, UINT32, "NN_P15tct_T2BMI", "NN_P15tct_T2BM");    \
        EXT_STR_ITEM_INFO_ZZP(                                                                                    \
            ok, si, offset, struct_t, printerr, NN_P15tct_T2BME, UINT32, "NN_P15tct_T2BME", "NN_P15tct_T2BM");    \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, NN_P15tct_T2B, UINT32, "NN_P15tct_T2B", 10000); \
        EXT_STR_ITEM_INFO_ZZP(                                                                                    \
            ok, si, offset, struct_t, printerr, NN_P15tct_T2Bv, UINT32, "NN_P15tct_T2Bv", "NN_P15tct_T2B");       \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, NN_P15tft_T1BM, UINT32, "NN_P15tft_T1BM", 50);  \
        EXT_STR_ITEM_INFO_ZZP(                                                                                    \
            ok, si, offset, struct_t, printerr, NN_P15tft_T1BMI, UINT32, "NN_P15tft_T1BMI", "NN_P15tft_T1BM");    \
        EXT_STR_ITEM_INFO_ZZP(                                                                                    \
            ok, si, offset, struct_t, printerr, NN_P15tft_T1BME, UINT32, "NN_P15tft_T1BME", "NN_P15tft_T1BM");    \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, NN_P15tft_T1B, UINT32, "NN_P15tft_T1B", 10000); \
        EXT_STR_ITEM_INFO_ZZP(                                                                                    \
            ok, si, offset, struct_t, printerr, NN_P15tft_T1Bv, UINT32, "NN_P15tft_T1Bv", "NN_P15tft_T1B");       \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, NN_P15tft_T2BM, UINT32, "NN_P15tft_T2BM", 50);  \
        EXT_STR_ITEM_INFO_ZZP(                                                                                    \
            ok, si, offset, struct_t, printerr, NN_P15tft_T2BMI, UINT32, "NN_P15tft_T2BMI", "NN_P15tft_T2BM");    \
        EXT_STR_ITEM_INFO_ZZP(                                                                                    \
            ok, si, offset, struct_t, printerr, NN_P15tft_T2BME, UINT32, "NN_P15tft_T2BME", "NN_P15tft_T2BM");    \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, NN_P15tft_T2B, UINT32, "NN_P15tft_T2B", 10000); \
        EXT_STR_ITEM_INFO_ZZP(                                                                                    \
            ok, si, offset, struct_t, printerr, NN_P15tft_T2Bv, UINT32, "NN_P15tft_T2Bv", "NN_P15tft_T2B");       \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, NN_P16tcl_T1BM, UINT32, "NN_P16tcl_T1BM", 50);  \
        EXT_STR_ITEM_INFO_ZZP(                                                                                    \
            ok, si, offset, struct_t, printerr, NN_P16tcl_T1BMI, UINT32, "NN_P16tcl_T1BMI", "NN_P16tcl_T1BM");    \
        EXT_STR_ITEM_INFO_ZZP(                                                                                    \
            ok, si, offset, struct_t, printerr, NN_P16tcl_T1BME, UINT32, "NN_P16tcl_T1BME", "NN_P16tcl_T1BM");    \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, NN_P16tcl_T1B, UINT32, "NN_P16tcl_T1B", 10000); \
        EXT_STR_ITEM_INFO_ZZP(                                                                                    \
            ok, si, offset, struct_t, printerr, NN_P16tcl_T1Bv, UINT32, "NN_P16tcl_T1Bv", "NN_P16tcl_T1B");       \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, NN_P16tcl_T2BM, UINT32, "NN_P16tcl_T2BM", 50);  \
        EXT_STR_ITEM_INFO_ZZP(                                                                                    \
            ok, si, offset, struct_t, printerr, NN_P16tcl_T2BMI, UINT32, "NN_P16tcl_T2BMI", "NN_P16tcl_T2BM");    \
        EXT_STR_ITEM_INFO_ZZP(                                                                                    \
            ok, si, offset, struct_t, printerr, NN_P16tcl_T2BME, UINT32, "NN_P16tcl_T2BME", "NN_P16tcl_T2BM");    \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, NN_P16tcl_T2B, UINT32, "NN_P16tcl_T2B", 10000); \
        EXT_STR_ITEM_INFO_ZZP(                                                                                    \
            ok, si, offset, struct_t, printerr, NN_P16tcl_T2Bv, UINT32, "NN_P16tcl_T2Bv", "NN_P16tcl_T2B");       \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, NN_P16tfl_T1BM, UINT32, "NN_P16tfl_T1BM", 50);  \
        EXT_STR_ITEM_INFO_ZZP(                                                                                    \
            ok, si, offset, struct_t, printerr, NN_P16tfl_T1BMI, UINT32, "NN_P16tfl_T1BMI", "NN_P16tfl_T1BM");    \
        EXT_STR_ITEM_INFO_ZZP(                                                                                    \
            ok, si, offset, struct_t, printerr, NN_P16tfl_T1BME, UINT32, "NN_P16tfl_T1BME", "NN_P16tfl_T1BM");    \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, NN_P16tfl_T1B, UINT32, "NN_P16tfl_T1B", 10000); \
        EXT_STR_ITEM_INFO_ZZP(                                                                                    \
            ok, si, offset, struct_t, printerr, NN_P16tfl_T1Bv, UINT32, "NN_P16tfl_T1Bv", "NN_P16tfl_T1B");       \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, NN_P16tfl_T2BM, UINT32, "NN_P16tfl_T2BM", 50);  \
        EXT_STR_ITEM_INFO_ZZP(                                                                                    \
            ok, si, offset, struct_t, printerr, NN_P16tfl_T2BMI, UINT32, "NN_P16tfl_T2BMI", "NN_P16tfl_T2BM");    \
        EXT_STR_ITEM_INFO_ZZP(                                                                                    \
            ok, si, offset, struct_t, printerr, NN_P16tfl_T2BME, UINT32, "NN_P16tfl_T2BME", "NN_P16tfl_T2BM");    \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, NN_P16tfl_T2B, UINT32, "NN_P16tfl_T2B", 10000); \
        EXT_STR_ITEM_INFO_ZZP(                                                                                    \
            ok, si, offset, struct_t, printerr, NN_P16tfl_T2Bv, UINT32, "NN_P16tfl_T2Bv", "NN_P16tfl_T2B");       \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, NN_P16tct_T1BM, UINT32, "NN_P16tct_T1BM", 50);  \
        EXT_STR_ITEM_INFO_ZZP(                                                                                    \
            ok, si, offset, struct_t, printerr, NN_P16tct_T1BMI, UINT32, "NN_P16tct_T1BMI", "NN_P16tct_T1BM");    \
        EXT_STR_ITEM_INFO_ZZP(                                                                                    \
            ok, si, offset, struct_t, printerr, NN_P16tct_T1BME, UINT32, "NN_P16tct_T1BME", "NN_P16tct_T1BM");    \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, NN_P16tct_T1B, UINT32, "NN_P16tct_T1B", 10000); \
        EXT_STR_ITEM_INFO_ZZP(                                                                                    \
            ok, si, offset, struct_t, printerr, NN_P16tct_T1Bv, UINT32, "NN_P16tct_T1Bv", "NN_P16tct_T1B");       \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, NN_P16tct_T2BM, UINT32, "NN_P16tct_T2BM", 50);  \
        EXT_STR_ITEM_INFO_ZZP(                                                                                    \
            ok, si, offset, struct_t, printerr, NN_P16tct_T2BMI, UINT32, "NN_P16tct_T2BMI", "NN_P16tct_T2BM");    \
        EXT_STR_ITEM_INFO_ZZP(                                                                                    \
            ok, si, offset, struct_t, printerr, NN_P16tct_T2BME, UINT32, "NN_P16tct_T2BME", "NN_P16tct_T2BM");    \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, NN_P16tct_T2B, UINT32, "NN_P16tct_T2B", 10000); \
        EXT_STR_ITEM_INFO_ZZP(                                                                                    \
            ok, si, offset, struct_t, printerr, NN_P16tct_T2Bv, UINT32, "NN_P16tct_T2Bv", "NN_P16tct_T2B");       \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, NN_P16tft_T1BM, UINT32, "NN_P16tft_T1BM", 50);  \
        EXT_STR_ITEM_INFO_ZZP(                                                                                    \
            ok, si, offset, struct_t, printerr, NN_P16tft_T1BMI, UINT32, "NN_P16tft_T1BMI", "NN_P16tft_T1BM");    \
        EXT_STR_ITEM_INFO_ZZP(                                                                                    \
            ok, si, offset, struct_t, printerr, NN_P16tft_T1BME, UINT32, "NN_P16tft_T1BME", "NN_P16tft_T1BM");    \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, NN_P16tft_T1B, UINT32, "NN_P16tft_T1B", 10000); \
        EXT_STR_ITEM_INFO_ZZP(                                                                                    \
            ok, si, offset, struct_t, printerr, NN_P16tft_T1Bv, UINT32, "NN_P16tft_T1Bv", "NN_P16tft_T1B");       \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, NN_P16tft_T2BM, UINT32, "NN_P16tft_T2BM", 50);  \
        EXT_STR_ITEM_INFO_ZZP(                                                                                    \
            ok, si, offset, struct_t, printerr, NN_P16tft_T2BMI, UINT32, "NN_P16tft_T2BMI", "NN_P16tft_T2BM");    \
        EXT_STR_ITEM_INFO_ZZP(                                                                                    \
            ok, si, offset, struct_t, printerr, NN_P16tft_T2BME, UINT32, "NN_P16tft_T2BME", "NN_P16tft_T2BM");    \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, NN_P16tft_T2B, UINT32, "NN_P16tft_T2B", 10000); \
        EXT_STR_ITEM_INFO_ZZP(                                                                                    \
            ok, si, offset, struct_t, printerr, NN_P16tft_T2Bv, UINT32, "NN_P16tft_T2Bv", "NN_P16tft_T2B");       \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, NN_P17tcl_T1BM, UINT32, "NN_P17tcl_T1BM", 50);  \
        EXT_STR_ITEM_INFO_ZZP(                                                                                    \
            ok, si, offset, struct_t, printerr, NN_P17tcl_T1BMI, UINT32, "NN_P17tcl_T1BMI", "NN_P17tcl_T1BM");    \
        EXT_STR_ITEM_INFO_ZZP(                                                                                    \
            ok, si, offset, struct_t, printerr, NN_P17tcl_T1BME, UINT32, "NN_P17tcl_T1BME", "NN_P17tcl_T1BM");    \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, NN_P17tcl_T1B, UINT32, "NN_P17tcl_T1B", 10000); \
        EXT_STR_ITEM_INFO_ZZP(                                                                                    \
            ok, si, offset, struct_t, printerr, NN_P17tcl_T1Bv, UINT32, "NN_P17tcl_T1Bv", "NN_P17tcl_T1B");       \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, NN_P17tcl_T2BM, UINT32, "NN_P17tcl_T2BM", 50);  \
        EXT_STR_ITEM_INFO_ZZP(                                                                                    \
            ok, si, offset, struct_t, printerr, NN_P17tcl_T2BMI, UINT32, "NN_P17tcl_T2BMI", "NN_P17tcl_T2BM");    \
        EXT_STR_ITEM_INFO_ZZP(                                                                                    \
            ok, si, offset, struct_t, printerr, NN_P17tcl_T2BME, UINT32, "NN_P17tcl_T2BME", "NN_P17tcl_T2BM");    \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, NN_P17tcl_T2B, UINT32, "NN_P17tcl_T2B", 10000); \
        EXT_STR_ITEM_INFO_ZZP(                                                                                    \
            ok, si, offset, struct_t, printerr, NN_P17tcl_T2Bv, UINT32, "NN_P17tcl_T2Bv", "NN_P17tcl_T2B");       \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, NN_P17tfl_T1BM, UINT32, "NN_P17tfl_T1BM", 50);  \
        EXT_STR_ITEM_INFO_ZZP(                                                                                    \
            ok, si, offset, struct_t, printerr, NN_P17tfl_T1BMI, UINT32, "NN_P17tfl_T1BMI", "NN_P17tfl_T1BM");    \
        EXT_STR_ITEM_INFO_ZZP(                                                                                    \
            ok, si, offset, struct_t, printerr, NN_P17tfl_T1BME, UINT32, "NN_P17tfl_T1BME", "NN_P17tfl_T1BM");    \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, NN_P17tfl_T1B, UINT32, "NN_P17tfl_T1B", 10000); \
        EXT_STR_ITEM_INFO_ZZP(                                                                                    \
            ok, si, offset, struct_t, printerr, NN_P17tfl_T1Bv, UINT32, "NN_P17tfl_T1Bv", "NN_P17tfl_T1B");       \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, NN_P17tfl_T2BM, UINT32, "NN_P17tfl_T2BM", 50);  \
        EXT_STR_ITEM_INFO_ZZP(                                                                                    \
            ok, si, offset, struct_t, printerr, NN_P17tfl_T2BMI, UINT32, "NN_P17tfl_T2BMI", "NN_P17tfl_T2BM");    \
        EXT_STR_ITEM_INFO_ZZP(                                                                                    \
            ok, si, offset, struct_t, printerr, NN_P17tfl_T2BME, UINT32, "NN_P17tfl_T2BME", "NN_P17tfl_T2BM");    \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, NN_P17tfl_T2B, UINT32, "NN_P17tfl_T2B", 10000); \
        EXT_STR_ITEM_INFO_ZZP(                                                                                    \
            ok, si, offset, struct_t, printerr, NN_P17tfl_T2Bv, UINT32, "NN_P17tfl_T2Bv", "NN_P17tfl_T2B");       \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, NN_P17tct_T1BM, UINT32, "NN_P17tct_T1BM", 50);  \
        EXT_STR_ITEM_INFO_ZZP(                                                                                    \
            ok, si, offset, struct_t, printerr, NN_P17tct_T1BMI, UINT32, "NN_P17tct_T1BMI", "NN_P17tct_T1BM");    \
        EXT_STR_ITEM_INFO_ZZP(                                                                                    \
            ok, si, offset, struct_t, printerr, NN_P17tct_T1BME, UINT32, "NN_P17tct_T1BME", "NN_P17tct_T1BM");    \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, NN_P17tct_T1B, UINT32, "NN_P17tct_T1B", 10000); \
        EXT_STR_ITEM_INFO_ZZP(                                                                                    \
            ok, si, offset, struct_t, printerr, NN_P17tct_T1Bv, UINT32, "NN_P17tct_T1Bv", "NN_P17tct_T1B");       \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, NN_P17tct_T2BM, UINT32, "NN_P17tct_T2BM", 50);  \
        EXT_STR_ITEM_INFO_ZZP(                                                                                    \
            ok, si, offset, struct_t, printerr, NN_P17tct_T2BMI, UINT32, "NN_P17tct_T2BMI", "NN_P17tct_T2BM");    \
        EXT_STR_ITEM_INFO_ZZP(                                                                                    \
            ok, si, offset, struct_t, printerr, NN_P17tct_T2BME, UINT32, "NN_P17tct_T2BME", "NN_P17tct_T2BM");    \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, NN_P17tct_T2B, UINT32, "NN_P17tct_T2B", 10000); \
        EXT_STR_ITEM_INFO_ZZP(                                                                                    \
            ok, si, offset, struct_t, printerr, NN_P17tct_T2Bv, UINT32, "NN_P17tct_T2Bv", "NN_P17tct_T2B");       \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, NN_P17tft_T1BM, UINT32, "NN_P17tft_T1BM", 50);  \
        EXT_STR_ITEM_INFO_ZZP(                                                                                    \
            ok, si, offset, struct_t, printerr, NN_P17tft_T1BMI, UINT32, "NN_P17tft_T1BMI", "NN_P17tft_T1BM");    \
        EXT_STR_ITEM_INFO_ZZP(                                                                                    \
            ok, si, offset, struct_t, printerr, NN_P17tft_T1BME, UINT32, "NN_P17tft_T1BME", "NN_P17tft_T1BM");    \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, NN_P17tft_T1B, UINT32, "NN_P17tft_T1B", 10000); \
        EXT_STR_ITEM_INFO_ZZP(                                                                                    \
            ok, si, offset, struct_t, printerr, NN_P17tft_T1Bv, UINT32, "NN_P17tft_T1Bv", "NN_P17tft_T1B");       \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, NN_P17tft_T2BM, UINT32, "NN_P17tft_T2BM", 50);  \
        EXT_STR_ITEM_INFO_ZZP(                                                                                    \
            ok, si, offset, struct_t, printerr, NN_P17tft_T2BMI, UINT32, "NN_P17tft_T2BMI", "NN_P17tft_T2BM");    \
        EXT_STR_ITEM_INFO_ZZP(                                                                                    \
            ok, si, offset, struct_t, printerr, NN_P17tft_T2BME, UINT32, "NN_P17tft_T2BME", "NN_P17tft_T2BM");    \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, NN_P17tft_T2B, UINT32, "NN_P17tft_T2B", 10000); \
        EXT_STR_ITEM_INFO_ZZP(                                                                                    \
            ok, si, offset, struct_t, printerr, NN_P17tft_T2Bv, UINT32, "NN_P17tft_T2Bv", "NN_P17tft_T2B");       \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, NN_P18tcl_T1BM, UINT32, "NN_P18tcl_T1BM", 50);  \
        EXT_STR_ITEM_INFO_ZZP(                                                                                    \
            ok, si, offset, struct_t, printerr, NN_P18tcl_T1BMI, UINT32, "NN_P18tcl_T1BMI", "NN_P18tcl_T1BM");    \
        EXT_STR_ITEM_INFO_ZZP(                                                                                    \
            ok, si, offset, struct_t, printerr, NN_P18tcl_T1BME, UINT32, "NN_P18tcl_T1BME", "NN_P18tcl_T1BM");    \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, NN_P18tcl_T1B, UINT32, "NN_P18tcl_T1B", 10000); \
        EXT_STR_ITEM_INFO_ZZP(                                                                                    \
            ok, si, offset, struct_t, printerr, NN_P18tcl_T1Bv, UINT32, "NN_P18tcl_T1Bv", "NN_P18tcl_T1B");       \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, NN_P18tcl_T2BM, UINT32, "NN_P18tcl_T2BM", 50);  \
        EXT_STR_ITEM_INFO_ZZP(                                                                                    \
            ok, si, offset, struct_t, printerr, NN_P18tcl_T2BMI, UINT32, "NN_P18tcl_T2BMI", "NN_P18tcl_T2BM");    \
        EXT_STR_ITEM_INFO_ZZP(                                                                                    \
            ok, si, offset, struct_t, printerr, NN_P18tcl_T2BME, UINT32, "NN_P18tcl_T2BME", "NN_P18tcl_T2BM");    \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, NN_P18tcl_T2B, UINT32, "NN_P18tcl_T2B", 10000); \
        EXT_STR_ITEM_INFO_ZZP(                                                                                    \
            ok, si, offset, struct_t, printerr, NN_P18tcl_T2Bv, UINT32, "NN_P18tcl_T2Bv", "NN_P18tcl_T2B");       \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, NN_P18tfl_T1BM, UINT32, "NN_P18tfl_T1BM", 50);  \
        EXT_STR_ITEM_INFO_ZZP(                                                                                    \
            ok, si, offset, struct_t, printerr, NN_P18tfl_T1BMI, UINT32, "NN_P18tfl_T1BMI", "NN_P18tfl_T1BM");    \
        EXT_STR_ITEM_INFO_ZZP(                                                                                    \
            ok, si, offset, struct_t, printerr, NN_P18tfl_T1BME, UINT32, "NN_P18tfl_T1BME", "NN_P18tfl_T1BM");    \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, NN_P18tfl_T1B, UINT32, "NN_P18tfl_T1B", 10000); \
        EXT_STR_ITEM_INFO_ZZP(                                                                                    \
            ok, si, offset, struct_t, printerr, NN_P18tfl_T1Bv, UINT32, "NN_P18tfl_T1Bv", "NN_P18tfl_T1B");       \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, NN_P18tfl_T2BM, UINT32, "NN_P18tfl_T2BM", 50);  \
        EXT_STR_ITEM_INFO_ZZP(                                                                                    \
            ok, si, offset, struct_t, printerr, NN_P18tfl_T2BMI, UINT32, "NN_P18tfl_T2BMI", "NN_P18tfl_T2BM");    \
        EXT_STR_ITEM_INFO_ZZP(                                                                                    \
            ok, si, offset, struct_t, printerr, NN_P18tfl_T2BME, UINT32, "NN_P18tfl_T2BME", "NN_P18tfl_T2BM");    \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, NN_P18tfl_T2B, UINT32, "NN_P18tfl_T2B", 10000); \
        EXT_STR_ITEM_INFO_ZZP(                                                                                    \
            ok, si, offset, struct_t, printerr, NN_P18tfl_T2Bv, UINT32, "NN_P18tfl_T2Bv", "NN_P18tfl_T2B");       \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, NN_P18tct_T1BM, UINT32, "NN_P18tct_T1BM", 50);  \
        EXT_STR_ITEM_INFO_ZZP(                                                                                    \
            ok, si, offset, struct_t, printerr, NN_P18tct_T1BMI, UINT32, "NN_P18tct_T1BMI", "NN_P18tct_T1BM");    \
        EXT_STR_ITEM_INFO_ZZP(                                                                                    \
            ok, si, offset, struct_t, printerr, NN_P18tct_T1BME, UINT32, "NN_P18tct_T1BME", "NN_P18tct_T1BM");    \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, NN_P18tct_T1B, UINT32, "NN_P18tct_T1B", 10000); \
        EXT_STR_ITEM_INFO_ZZP(                                                                                    \
            ok, si, offset, struct_t, printerr, NN_P18tct_T1Bv, UINT32, "NN_P18tct_T1Bv", "NN_P18tct_T1B");       \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, NN_P18tct_T2BM, UINT32, "NN_P18tct_T2BM", 50);  \
        EXT_STR_ITEM_INFO_ZZP(                                                                                    \
            ok, si, offset, struct_t, printerr, NN_P18tct_T2BMI, UINT32, "NN_P18tct_T2BMI", "NN_P18tct_T2BM");    \
        EXT_STR_ITEM_INFO_ZZP(                                                                                    \
            ok, si, offset, struct_t, printerr, NN_P18tct_T2BME, UINT32, "NN_P18tct_T2BME", "NN_P18tct_T2BM");    \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, NN_P18tct_T2B, UINT32, "NN_P18tct_T2B", 10000); \
        EXT_STR_ITEM_INFO_ZZP(                                                                                    \
            ok, si, offset, struct_t, printerr, NN_P18tct_T2Bv, UINT32, "NN_P18tct_T2Bv", "NN_P18tct_T2B");       \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, NN_P18tft_T1BM, UINT32, "NN_P18tft_T1BM", 50);  \
        EXT_STR_ITEM_INFO_ZZP(                                                                                    \
            ok, si, offset, struct_t, printerr, NN_P18tft_T1BMI, UINT32, "NN_P18tft_T1BMI", "NN_P18tft_T1BM");    \
        EXT_STR_ITEM_INFO_ZZP(                                                                                    \
            ok, si, offset, struct_t, printerr, NN_P18tft_T1BME, UINT32, "NN_P18tft_T1BME", "NN_P18tft_T1BM");    \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, NN_P18tft_T1B, UINT32, "NN_P18tft_T1B", 10000); \
        EXT_STR_ITEM_INFO_ZZP(                                                                                    \
            ok, si, offset, struct_t, printerr, NN_P18tft_T1Bv, UINT32, "NN_P18tft_T1Bv", "NN_P18tft_T1B");       \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, NN_P18tft_T2BM, UINT32, "NN_P18tft_T2BM", 50);  \
        EXT_STR_ITEM_INFO_ZZP(                                                                                    \
            ok, si, offset, struct_t, printerr, NN_P18tft_T2BMI, UINT32, "NN_P18tft_T2BMI", "NN_P18tft_T2BM");    \
        EXT_STR_ITEM_INFO_ZZP(                                                                                    \
            ok, si, offset, struct_t, printerr, NN_P18tft_T2BME, UINT32, "NN_P18tft_T2BME", "NN_P18tft_T2BM");    \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, NN_P18tft_T2B, UINT32, "NN_P18tft_T2B", 10000); \
        EXT_STR_ITEM_INFO_ZZP(                                                                                    \
            ok, si, offset, struct_t, printerr, NN_P18tft_T2Bv, UINT32, "NN_P18tft_T2Bv", "NN_P18tft_T2B");       \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, NN_P19tcl_T1BM, UINT32, "NN_P19tcl_T1BM", 50);  \
        EXT_STR_ITEM_INFO_ZZP(                                                                                    \
            ok, si, offset, struct_t, printerr, NN_P19tcl_T1BMI, UINT32, "NN_P19tcl_T1BMI", "NN_P19tcl_T1BM");    \
        EXT_STR_ITEM_INFO_ZZP(                                                                                    \
            ok, si, offset, struct_t, printerr, NN_P19tcl_T1BME, UINT32, "NN_P19tcl_T1BME", "NN_P19tcl_T1BM");    \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, NN_P19tcl_T1B, UINT32, "NN_P19tcl_T1B", 10000); \
        EXT_STR_ITEM_INFO_ZZP(                                                                                    \
            ok, si, offset, struct_t, printerr, NN_P19tcl_T1Bv, UINT32, "NN_P19tcl_T1Bv", "NN_P19tcl_T1B");       \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, NN_P19tcl_T2BM, UINT32, "NN_P19tcl_T2BM", 50);  \
        EXT_STR_ITEM_INFO_ZZP(                                                                                    \
            ok, si, offset, struct_t, printerr, NN_P19tcl_T2BMI, UINT32, "NN_P19tcl_T2BMI", "NN_P19tcl_T2BM");    \
        EXT_STR_ITEM_INFO_ZZP(                                                                                    \
            ok, si, offset, struct_t, printerr, NN_P19tcl_T2BME, UINT32, "NN_P19tcl_T2BME", "NN_P19tcl_T2BM");    \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, NN_P19tcl_T2B, UINT32, "NN_P19tcl_T2B", 10000); \
        EXT_STR_ITEM_INFO_ZZP(                                                                                    \
            ok, si, offset, struct_t, printerr, NN_P19tcl_T2Bv, UINT32, "NN_P19tcl_T2Bv", "NN_P19tcl_T2B");       \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, NN_P19tfl_T1BM, UINT32, "NN_P19tfl_T1BM", 50);  \
        EXT_STR_ITEM_INFO_ZZP(                                                                                    \
            ok, si, offset, struct_t, printerr, NN_P19tfl_T1BMI, UINT32, "NN_P19tfl_T1BMI", "NN_P19tfl_T1BM");    \
        EXT_STR_ITEM_INFO_ZZP(                                                                                    \
            ok, si, offset, struct_t, printerr, NN_P19tfl_T1BME, UINT32, "NN_P19tfl_T1BME", "NN_P19tfl_T1BM");    \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, NN_P19tfl_T1B, UINT32, "NN_P19tfl_T1B", 10000); \
        EXT_STR_ITEM_INFO_ZZP(                                                                                    \
            ok, si, offset, struct_t, printerr, NN_P19tfl_T1Bv, UINT32, "NN_P19tfl_T1Bv", "NN_P19tfl_T1B");       \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, NN_P19tfl_T2BM, UINT32, "NN_P19tfl_T2BM", 50);  \
        EXT_STR_ITEM_INFO_ZZP(                                                                                    \
            ok, si, offset, struct_t, printerr, NN_P19tfl_T2BMI, UINT32, "NN_P19tfl_T2BMI", "NN_P19tfl_T2BM");    \
        EXT_STR_ITEM_INFO_ZZP(                                                                                    \
            ok, si, offset, struct_t, printerr, NN_P19tfl_T2BME, UINT32, "NN_P19tfl_T2BME", "NN_P19tfl_T2BM");    \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, NN_P19tfl_T2B, UINT32, "NN_P19tfl_T2B", 10000); \
        EXT_STR_ITEM_INFO_ZZP(                                                                                    \
            ok, si, offset, struct_t, printerr, NN_P19tfl_T2Bv, UINT32, "NN_P19tfl_T2Bv", "NN_P19tfl_T2B");       \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, NN_P19tct_T1BM, UINT32, "NN_P19tct_T1BM", 50);  \
        EXT_STR_ITEM_INFO_ZZP(                                                                                    \
            ok, si, offset, struct_t, printerr, NN_P19tct_T1BMI, UINT32, "NN_P19tct_T1BMI", "NN_P19tct_T1BM");    \
        EXT_STR_ITEM_INFO_ZZP(                                                                                    \
            ok, si, offset, struct_t, printerr, NN_P19tct_T1BME, UINT32, "NN_P19tct_T1BME", "NN_P19tct_T1BM");    \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, NN_P19tct_T1B, UINT32, "NN_P19tct_T1B", 10000); \
        EXT_STR_ITEM_INFO_ZZP(                                                                                    \
            ok, si, offset, struct_t, printerr, NN_P19tct_T1Bv, UINT32, "NN_P19tct_T1Bv", "NN_P19tct_T1B");       \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, NN_P19tct_T2BM, UINT32, "NN_P19tct_T2BM", 50);  \
        EXT_STR_ITEM_INFO_ZZP(                                                                                    \
            ok, si, offset, struct_t, printerr, NN_P19tct_T2BMI, UINT32, "NN_P19tct_T2BMI", "NN_P19tct_T2BM");    \
        EXT_STR_ITEM_INFO_ZZP(                                                                                    \
            ok, si, offset, struct_t, printerr, NN_P19tct_T2BME, UINT32, "NN_P19tct_T2BME", "NN_P19tct_T2BM");    \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, NN_P19tct_T2B, UINT32, "NN_P19tct_T2B", 10000); \
        EXT_STR_ITEM_INFO_ZZP(                                                                                    \
            ok, si, offset, struct_t, printerr, NN_P19tct_T2Bv, UINT32, "NN_P19tct_T2Bv", "NN_P19tct_T2B");       \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, NN_P19tft_T1BM, UINT32, "NN_P19tft_T1BM", 50);  \
        EXT_STR_ITEM_INFO_ZZP(                                                                                    \
            ok, si, offset, struct_t, printerr, NN_P19tft_T1BMI, UINT32, "NN_P19tft_T1BMI", "NN_P19tft_T1BM");    \
        EXT_STR_ITEM_INFO_ZZP(                                                                                    \
            ok, si, offset, struct_t, printerr, NN_P19tft_T1BME, UINT32, "NN_P19tft_T1BME", "NN_P19tft_T1BM");    \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, NN_P19tft_T1B, UINT32, "NN_P19tft_T1B", 10000); \
        EXT_STR_ITEM_INFO_ZZP(                                                                                    \
            ok, si, offset, struct_t, printerr, NN_P19tft_T1Bv, UINT32, "NN_P19tft_T1Bv", "NN_P19tft_T1B");       \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, NN_P19tft_T2BM, UINT32, "NN_P19tft_T2BM", 50);  \
        EXT_STR_ITEM_INFO_ZZP(                                                                                    \
            ok, si, offset, struct_t, printerr, NN_P19tft_T2BMI, UINT32, "NN_P19tft_T2BMI", "NN_P19tft_T2BM");    \
        EXT_STR_ITEM_INFO_ZZP(                                                                                    \
            ok, si, offset, struct_t, printerr, NN_P19tft_T2BME, UINT32, "NN_P19tft_T2BME", "NN_P19tft_T2BM");    \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, NN_P19tft_T2B, UINT32, "NN_P19tft_T2B", 10000); \
        EXT_STR_ITEM_INFO_ZZP(                                                                                    \
            ok, si, offset, struct_t, printerr, NN_P19tft_T2Bv, UINT32, "NN_P19tft_T2Bv", "NN_P19tft_T2B");       \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, NN_P20tcl_T1BM, UINT32, "NN_P20tcl_T1BM", 50);  \
        EXT_STR_ITEM_INFO_ZZP(                                                                                    \
            ok, si, offset, struct_t, printerr, NN_P20tcl_T1BMI, UINT32, "NN_P20tcl_T1BMI", "NN_P20tcl_T1BM");    \
        EXT_STR_ITEM_INFO_ZZP(                                                                                    \
            ok, si, offset, struct_t, printerr, NN_P20tcl_T1BME, UINT32, "NN_P20tcl_T1BME", "NN_P20tcl_T1BM");    \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, NN_P20tcl_T1B, UINT32, "NN_P20tcl_T1B", 10000); \
        EXT_STR_ITEM_INFO_ZZP(                                                                                    \
            ok, si, offset, struct_t, printerr, NN_P20tcl_T1Bv, UINT32, "NN_P20tcl_T1Bv", "NN_P20tcl_T1B");       \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, NN_P20tcl_T2BM, UINT32, "NN_P20tcl_T2BM", 50);  \
        EXT_STR_ITEM_INFO_ZZP(                                                                                    \
            ok, si, offset, struct_t, printerr, NN_P20tcl_T2BMI, UINT32, "NN_P20tcl_T2BMI", "NN_P20tcl_T2BM");    \
        EXT_STR_ITEM_INFO_ZZP(                                                                                    \
            ok, si, offset, struct_t, printerr, NN_P20tcl_T2BME, UINT32, "NN_P20tcl_T2BME", "NN_P20tcl_T2BM");    \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, NN_P20tcl_T2B, UINT32, "NN_P20tcl_T2B", 10000); \
        EXT_STR_ITEM_INFO_ZZP(                                                                                    \
            ok, si, offset, struct_t, printerr, NN_P20tcl_T2Bv, UINT32, "NN_P20tcl_T2Bv", "NN_P20tcl_T2B");       \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, NN_P20tfl_T1BM, UINT32, "NN_P20tfl_T1BM", 50);  \
        EXT_STR_ITEM_INFO_ZZP(                                                                                    \
            ok, si, offset, struct_t, printerr, NN_P20tfl_T1BMI, UINT32, "NN_P20tfl_T1BMI", "NN_P20tfl_T1BM");    \
        EXT_STR_ITEM_INFO_ZZP(                                                                                    \
            ok, si, offset, struct_t, printerr, NN_P20tfl_T1BME, UINT32, "NN_P20tfl_T1BME", "NN_P20tfl_T1BM");    \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, NN_P20tfl_T1B, UINT32, "NN_P20tfl_T1B", 10000); \
        EXT_STR_ITEM_INFO_ZZP(                                                                                    \
            ok, si, offset, struct_t, printerr, NN_P20tfl_T1Bv, UINT32, "NN_P20tfl_T1Bv", "NN_P20tfl_T1B");       \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, NN_P20tfl_T2BM, UINT32, "NN_P20tfl_T2BM", 50);  \
        EXT_STR_ITEM_INFO_ZZP(                                                                                    \
            ok, si, offset, struct_t, printerr, NN_P20tfl_T2BMI, UINT32, "NN_P20tfl_T2BMI", "NN_P20tfl_T2BM");    \
        EXT_STR_ITEM_INFO_ZZP(                                                                                    \
            ok, si, offset, struct_t, printerr, NN_P20tfl_T2BME, UINT32, "NN_P20tfl_T2BME", "NN_P20tfl_T2BM");    \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, NN_P20tfl_T2B, UINT32, "NN_P20tfl_T2B", 10000); \
        EXT_STR_ITEM_INFO_ZZP(                                                                                    \
            ok, si, offset, struct_t, printerr, NN_P20tfl_T2Bv, UINT32, "NN_P20tfl_T2Bv", "NN_P20tfl_T2B");       \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, NN_P20tct_T1BM, UINT32, "NN_P20tct_T1BM", 50);  \
        EXT_STR_ITEM_INFO_ZZP(                                                                                    \
            ok, si, offset, struct_t, printerr, NN_P20tct_T1BMI, UINT32, "NN_P20tct_T1BMI", "NN_P20tct_T1BM");    \
        EXT_STR_ITEM_INFO_ZZP(                                                                                    \
            ok, si, offset, struct_t, printerr, NN_P20tct_T1BME, UINT32, "NN_P20tct_T1BME", "NN_P20tct_T1BM");    \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, NN_P20tct_T1B, UINT32, "NN_P20tct_T1B", 10000); \
        EXT_STR_ITEM_INFO_ZZP(                                                                                    \
            ok, si, offset, struct_t, printerr, NN_P20tct_T1Bv, UINT32, "NN_P20tct_T1Bv", "NN_P20tct_T1B");       \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, NN_P20tct_T2BM, UINT32, "NN_P20tct_T2BM", 50);  \
        EXT_STR_ITEM_INFO_ZZP(                                                                                    \
            ok, si, offset, struct_t, printerr, NN_P20tct_T2BMI, UINT32, "NN_P20tct_T2BMI", "NN_P20tct_T2BM");    \
        EXT_STR_ITEM_INFO_ZZP(                                                                                    \
            ok, si, offset, struct_t, printerr, NN_P20tct_T2BME, UINT32, "NN_P20tct_T2BME", "NN_P20tct_T2BM");    \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, NN_P20tct_T2B, UINT32, "NN_P20tct_T2B", 10000); \
        EXT_STR_ITEM_INFO_ZZP(                                                                                    \
            ok, si, offset, struct_t, printerr, NN_P20tct_T2Bv, UINT32, "NN_P20tct_T2Bv", "NN_P20tct_T2B");       \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, NN_P20tft_T1BM, UINT32, "NN_P20tft_T1BM", 50);  \
        EXT_STR_ITEM_INFO_ZZP(                                                                                    \
            ok, si, offset, struct_t, printerr, NN_P20tft_T1BMI, UINT32, "NN_P20tft_T1BMI", "NN_P20tft_T1BM");    \
        EXT_STR_ITEM_INFO_ZZP(                                                                                    \
            ok, si, offset, struct_t, printerr, NN_P20tft_T1BME, UINT32, "NN_P20tft_T1BME", "NN_P20tft_T1BM");    \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, NN_P20tft_T1B, UINT32, "NN_P20tft_T1B", 10000); \
        EXT_STR_ITEM_INFO_ZZP(                                                                                    \
            ok, si, offset, struct_t, printerr, NN_P20tft_T1Bv, UINT32, "NN_P20tft_T1Bv", "NN_P20tft_T1B");       \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, NN_P20tft_T2BM, UINT32, "NN_P20tft_T2BM", 50);  \
        EXT_STR_ITEM_INFO_ZZP(                                                                                    \
            ok, si, offset, struct_t, printerr, NN_P20tft_T2BMI, UINT32, "NN_P20tft_T2BMI", "NN_P20tft_T2BM");    \
        EXT_STR_ITEM_INFO_ZZP(                                                                                    \
            ok, si, offset, struct_t, printerr, NN_P20tft_T2BME, UINT32, "NN_P20tft_T2BME", "NN_P20tft_T2BM");    \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, NN_P20tft_T2B, UINT32, "NN_P20tft_T2B", 10000); \
        EXT_STR_ITEM_INFO_ZZP(                                                                                    \
            ok, si, offset, struct_t, printerr, NN_P20tft_T2Bv, UINT32, "NN_P20tft_T2Bv", "NN_P20tft_T2B");       \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, NN_P21tcl_T1BM, UINT32, "NN_P21tcl_T1BM", 50);  \
        EXT_STR_ITEM_INFO_ZZP(                                                                                    \
            ok, si, offset, struct_t, printerr, NN_P21tcl_T1BMI, UINT32, "NN_P21tcl_T1BMI", "NN_P21tcl_T1BM");    \
        EXT_STR_ITEM_INFO_ZZP(                                                                                    \
            ok, si, offset, struct_t, printerr, NN_P21tcl_T1BME, UINT32, "NN_P21tcl_T1BME", "NN_P21tcl_T1BM");    \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, NN_P21tcl_T1B, UINT32, "NN_P21tcl_T1B", 10000); \
        EXT_STR_ITEM_INFO_ZZP(                                                                                    \
            ok, si, offset, struct_t, printerr, NN_P21tcl_T1Bv, UINT32, "NN_P21tcl_T1Bv", "NN_P21tcl_T1B");       \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, NN_P21tcl_T2BM, UINT32, "NN_P21tcl_T2BM", 50);  \
        EXT_STR_ITEM_INFO_ZZP(                                                                                    \
            ok, si, offset, struct_t, printerr, NN_P21tcl_T2BMI, UINT32, "NN_P21tcl_T2BMI", "NN_P21tcl_T2BM");    \
        EXT_STR_ITEM_INFO_ZZP(                                                                                    \
            ok, si, offset, struct_t, printerr, NN_P21tcl_T2BME, UINT32, "NN_P21tcl_T2BME", "NN_P21tcl_T2BM");    \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, NN_P21tcl_T2B, UINT32, "NN_P21tcl_T2B", 10000); \
        EXT_STR_ITEM_INFO_ZZP(                                                                                    \
            ok, si, offset, struct_t, printerr, NN_P21tcl_T2Bv, UINT32, "NN_P21tcl_T2Bv", "NN_P21tcl_T2B");       \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, NN_P21tfl_T1BM, UINT32, "NN_P21tfl_T1BM", 50);  \
        EXT_STR_ITEM_INFO_ZZP(                                                                                    \
            ok, si, offset, struct_t, printerr, NN_P21tfl_T1BMI, UINT32, "NN_P21tfl_T1BMI", "NN_P21tfl_T1BM");    \
        EXT_STR_ITEM_INFO_ZZP(                                                                                    \
            ok, si, offset, struct_t, printerr, NN_P21tfl_T1BME, UINT32, "NN_P21tfl_T1BME", "NN_P21tfl_T1BM");    \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, NN_P21tfl_T1B, UINT32, "NN_P21tfl_T1B", 10000); \
        EXT_STR_ITEM_INFO_ZZP(                                                                                    \
            ok, si, offset, struct_t, printerr, NN_P21tfl_T1Bv, UINT32, "NN_P21tfl_T1Bv", "NN_P21tfl_T1B");       \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, NN_P21tfl_T2BM, UINT32, "NN_P21tfl_T2BM", 50);  \
        EXT_STR_ITEM_INFO_ZZP(                                                                                    \
            ok, si, offset, struct_t, printerr, NN_P21tfl_T2BMI, UINT32, "NN_P21tfl_T2BMI", "NN_P21tfl_T2BM");    \
        EXT_STR_ITEM_INFO_ZZP(                                                                                    \
            ok, si, offset, struct_t, printerr, NN_P21tfl_T2BME, UINT32, "NN_P21tfl_T2BME", "NN_P21tfl_T2BM");    \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, NN_P21tfl_T2B, UINT32, "NN_P21tfl_T2B", 10000); \
        EXT_STR_ITEM_INFO_ZZP(                                                                                    \
            ok, si, offset, struct_t, printerr, NN_P21tfl_T2Bv, UINT32, "NN_P21tfl_T2Bv", "NN_P21tfl_T2B");       \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, NN_P21tct_T1BM, UINT32, "NN_P21tct_T1BM", 50);  \
        EXT_STR_ITEM_INFO_ZZP(                                                                                    \
            ok, si, offset, struct_t, printerr, NN_P21tct_T1BMI, UINT32, "NN_P21tct_T1BMI", "NN_P21tct_T1BM");    \
        EXT_STR_ITEM_INFO_ZZP(                                                                                    \
            ok, si, offset, struct_t, printerr, NN_P21tct_T1BME, UINT32, "NN_P21tct_T1BME", "NN_P21tct_T1BM");    \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, NN_P21tct_T1B, UINT32, "NN_P21tct_T1B", 10000); \
        EXT_STR_ITEM_INFO_ZZP(                                                                                    \
            ok, si, offset, struct_t, printerr, NN_P21tct_T1Bv, UINT32, "NN_P21tct_T1Bv", "NN_P21tct_T1B");       \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, NN_P21tct_T2BM, UINT32, "NN_P21tct_T2BM", 50);  \
        EXT_STR_ITEM_INFO_ZZP(                                                                                    \
            ok, si, offset, struct_t, printerr, NN_P21tct_T2BMI, UINT32, "NN_P21tct_T2BMI", "NN_P21tct_T2BM");    \
        EXT_STR_ITEM_INFO_ZZP(                                                                                    \
            ok, si, offset, struct_t, printerr, NN_P21tct_T2BME, UINT32, "NN_P21tct_T2BME", "NN_P21tct_T2BM");    \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, NN_P21tct_T2B, UINT32, "NN_P21tct_T2B", 10000); \
        EXT_STR_ITEM_INFO_ZZP(                                                                                    \
            ok, si, offset, struct_t, printerr, NN_P21tct_T2Bv, UINT32, "NN_P21tct_T2Bv", "NN_P21tct_T2B");       \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, NN_P21tft_T1BM, UINT32, "NN_P21tft_T1BM", 50);  \
        EXT_STR_ITEM_INFO_ZZP(                                                                                    \
            ok, si, offset, struct_t, printerr, NN_P21tft_T1BMI, UINT32, "NN_P21tft_T1BMI", "NN_P21tft_T1BM");    \
        EXT_STR_ITEM_INFO_ZZP(                                                                                    \
            ok, si, offset, struct_t, printerr, NN_P21tft_T1BME, UINT32, "NN_P21tft_T1BME", "NN_P21tft_T1BM");    \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, NN_P21tft_T1B, UINT32, "NN_P21tft_T1B", 10000); \
        EXT_STR_ITEM_INFO_ZZP(                                                                                    \
            ok, si, offset, struct_t, printerr, NN_P21tft_T1Bv, UINT32, "NN_P21tft_T1Bv", "NN_P21tft_T1B");       \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, NN_P21tft_T2BM, UINT32, "NN_P21tft_T2BM", 50);  \
        EXT_STR_ITEM_INFO_ZZP(                                                                                    \
            ok, si, offset, struct_t, printerr, NN_P21tft_T2BMI, UINT32, "NN_P21tft_T2BMI", "NN_P21tft_T2BM");    \
        EXT_STR_ITEM_INFO_ZZP(                                                                                    \
            ok, si, offset, struct_t, printerr, NN_P21tft_T2BME, UINT32, "NN_P21tft_T2BME", "NN_P21tft_T2BM");    \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, NN_P21tft_T2B, UINT32, "NN_P21tft_T2B", 10000); \
        EXT_STR_ITEM_INFO_ZZP(                                                                                    \
            ok, si, offset, struct_t, printerr, NN_P21tft_T2Bv, UINT32, "NN_P21tft_T2Bv", "NN_P21tft_T2B");       \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, NN_P22tcl_T1BM, UINT32, "NN_P22tcl_T1BM", 50);  \
        EXT_STR_ITEM_INFO_ZZP(                                                                                    \
            ok, si, offset, struct_t, printerr, NN_P22tcl_T1BMI, UINT32, "NN_P22tcl_T1BMI", "NN_P22tcl_T1BM");    \
        EXT_STR_ITEM_INFO_ZZP(                                                                                    \
            ok, si, offset, struct_t, printerr, NN_P22tcl_T1BME, UINT32, "NN_P22tcl_T1BME", "NN_P22tcl_T1BM");    \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, NN_P22tcl_T1B, UINT32, "NN_P22tcl_T1B", 10000); \
        EXT_STR_ITEM_INFO_ZZP(                                                                                    \
            ok, si, offset, struct_t, printerr, NN_P22tcl_T1Bv, UINT32, "NN_P22tcl_T1Bv", "NN_P22tcl_T1B");       \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, NN_P22tcl_T2BM, UINT32, "NN_P22tcl_T2BM", 50);  \
        EXT_STR_ITEM_INFO_ZZP(                                                                                    \
            ok, si, offset, struct_t, printerr, NN_P22tcl_T2BMI, UINT32, "NN_P22tcl_T2BMI", "NN_P22tcl_T2BM");    \
        EXT_STR_ITEM_INFO_ZZP(                                                                                    \
            ok, si, offset, struct_t, printerr, NN_P22tcl_T2BME, UINT32, "NN_P22tcl_T2BME", "NN_P22tcl_T2BM");    \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, NN_P22tcl_T2B, UINT32, "NN_P22tcl_T2B", 10000); \
        EXT_STR_ITEM_INFO_ZZP(                                                                                    \
            ok, si, offset, struct_t, printerr, NN_P22tcl_T2Bv, UINT32, "NN_P22tcl_T2Bv", "NN_P22tcl_T2B");       \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, NN_P22tfl_T1BM, UINT32, "NN_P22tfl_T1BM", 50);  \
        EXT_STR_ITEM_INFO_ZZP(                                                                                    \
            ok, si, offset, struct_t, printerr, NN_P22tfl_T1BMI, UINT32, "NN_P22tfl_T1BMI", "NN_P22tfl_T1BM");    \
        EXT_STR_ITEM_INFO_ZZP(                                                                                    \
            ok, si, offset, struct_t, printerr, NN_P22tfl_T1BME, UINT32, "NN_P22tfl_T1BME", "NN_P22tfl_T1BM");    \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, NN_P22tfl_T1B, UINT32, "NN_P22tfl_T1B", 10000); \
        EXT_STR_ITEM_INFO_ZZP(                                                                                    \
            ok, si, offset, struct_t, printerr, NN_P22tfl_T1Bv, UINT32, "NN_P22tfl_T1Bv", "NN_P22tfl_T1B");       \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, NN_P22tfl_T2BM, UINT32, "NN_P22tfl_T2BM", 50);  \
        EXT_STR_ITEM_INFO_ZZP(                                                                                    \
            ok, si, offset, struct_t, printerr, NN_P22tfl_T2BMI, UINT32, "NN_P22tfl_T2BMI", "NN_P22tfl_T2BM");    \
        EXT_STR_ITEM_INFO_ZZP(                                                                                    \
            ok, si, offset, struct_t, printerr, NN_P22tfl_T2BME, UINT32, "NN_P22tfl_T2BME", "NN_P22tfl_T2BM");    \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, NN_P22tfl_T2B, UINT32, "NN_P22tfl_T2B", 10000); \
        EXT_STR_ITEM_INFO_ZZP(                                                                                    \
            ok, si, offset, struct_t, printerr, NN_P22tfl_T2Bv, UINT32, "NN_P22tfl_T2Bv", "NN_P22tfl_T2B");       \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, NN_P22tct_T1BM, UINT32, "NN_P22tct_T1BM", 50);  \
        EXT_STR_ITEM_INFO_ZZP(                                                                                    \
            ok, si, offset, struct_t, printerr, NN_P22tct_T1BMI, UINT32, "NN_P22tct_T1BMI", "NN_P22tct_T1BM");    \
        EXT_STR_ITEM_INFO_ZZP(                                                                                    \
            ok, si, offset, struct_t, printerr, NN_P22tct_T1BME, UINT32, "NN_P22tct_T1BME", "NN_P22tct_T1BM");    \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, NN_P22tct_T1B, UINT32, "NN_P22tct_T1B", 10000); \
        EXT_STR_ITEM_INFO_ZZP(                                                                                    \
            ok, si, offset, struct_t, printerr, NN_P22tct_T1Bv, UINT32, "NN_P22tct_T1Bv", "NN_P22tct_T1B");       \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, NN_P22tct_T2BM, UINT32, "NN_P22tct_T2BM", 50);  \
        EXT_STR_ITEM_INFO_ZZP(                                                                                    \
            ok, si, offset, struct_t, printerr, NN_P22tct_T2BMI, UINT32, "NN_P22tct_T2BMI", "NN_P22tct_T2BM");    \
        EXT_STR_ITEM_INFO_ZZP(                                                                                    \
            ok, si, offset, struct_t, printerr, NN_P22tct_T2BME, UINT32, "NN_P22tct_T2BME", "NN_P22tct_T2BM");    \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, NN_P22tct_T2B, UINT32, "NN_P22tct_T2B", 10000); \
        EXT_STR_ITEM_INFO_ZZP(                                                                                    \
            ok, si, offset, struct_t, printerr, NN_P22tct_T2Bv, UINT32, "NN_P22tct_T2Bv", "NN_P22tct_T2B");       \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, NN_P22tft_T1BM, UINT32, "NN_P22tft_T1BM", 50);  \
        EXT_STR_ITEM_INFO_ZZP(                                                                                    \
            ok, si, offset, struct_t, printerr, NN_P22tft_T1BMI, UINT32, "NN_P22tft_T1BMI", "NN_P22tft_T1BM");    \
        EXT_STR_ITEM_INFO_ZZP(                                                                                    \
            ok, si, offset, struct_t, printerr, NN_P22tft_T1BME, UINT32, "NN_P22tft_T1BME", "NN_P22tft_T1BM");    \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, NN_P22tft_T1B, UINT32, "NN_P22tft_T1B", 10000); \
        EXT_STR_ITEM_INFO_ZZP(                                                                                    \
            ok, si, offset, struct_t, printerr, NN_P22tft_T1Bv, UINT32, "NN_P22tft_T1Bv", "NN_P22tft_T1B");       \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, NN_P22tft_T2BM, UINT32, "NN_P22tft_T2BM", 50);  \
        EXT_STR_ITEM_INFO_ZZP(                                                                                    \
            ok, si, offset, struct_t, printerr, NN_P22tft_T2BMI, UINT32, "NN_P22tft_T2BMI", "NN_P22tft_T2BM");    \
        EXT_STR_ITEM_INFO_ZZP(                                                                                    \
            ok, si, offset, struct_t, printerr, NN_P22tft_T2BME, UINT32, "NN_P22tft_T2BME", "NN_P22tft_T2BM");    \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, NN_P22tft_T2B, UINT32, "NN_P22tft_T2B", 10000); \
        EXT_STR_ITEM_INFO_ZZP(                                                                                    \
            ok, si, offset, struct_t, printerr, NN_P22tft_T2Bv, UINT32, "NN_P22tft_T2Bv", "NN_P22tft_T2B");       \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, NN_P23tcl_T1BM, UINT32, "NN_P23tcl_T1BM", 50);  \
        EXT_STR_ITEM_INFO_ZZP(                                                                                    \
            ok, si, offset, struct_t, printerr, NN_P23tcl_T1BMI, UINT32, "NN_P23tcl_T1BMI", "NN_P23tcl_T1BM");    \
        EXT_STR_ITEM_INFO_ZZP(                                                                                    \
            ok, si, offset, struct_t, printerr, NN_P23tcl_T1BME, UINT32, "NN_P23tcl_T1BME", "NN_P23tcl_T1BM");    \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, NN_P23tcl_T1B, UINT32, "NN_P23tcl_T1B", 10000); \
        EXT_STR_ITEM_INFO_ZZP(                                                                                    \
            ok, si, offset, struct_t, printerr, NN_P23tcl_T1Bv, UINT32, "NN_P23tcl_T1Bv", "NN_P23tcl_T1B");       \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, NN_P23tcl_T2BM, UINT32, "NN_P23tcl_T2BM", 50);  \
        EXT_STR_ITEM_INFO_ZZP(                                                                                    \
            ok, si, offset, struct_t, printerr, NN_P23tcl_T2BMI, UINT32, "NN_P23tcl_T2BMI", "NN_P23tcl_T2BM");    \
        EXT_STR_ITEM_INFO_ZZP(                                                                                    \
            ok, si, offset, struct_t, printerr, NN_P23tcl_T2BME, UINT32, "NN_P23tcl_T2BME", "NN_P23tcl_T2BM");    \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, NN_P23tcl_T2B, UINT32, "NN_P23tcl_T2B", 10000); \
        EXT_STR_ITEM_INFO_ZZP(                                                                                    \
            ok, si, offset, struct_t, printerr, NN_P23tcl_T2Bv, UINT32, "NN_P23tcl_T2Bv", "NN_P23tcl_T2B");       \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, NN_P23tfl_T1BM, UINT32, "NN_P23tfl_T1BM", 50);  \
        EXT_STR_ITEM_INFO_ZZP(                                                                                    \
            ok, si, offset, struct_t, printerr, NN_P23tfl_T1BMI, UINT32, "NN_P23tfl_T1BMI", "NN_P23tfl_T1BM");    \
        EXT_STR_ITEM_INFO_ZZP(                                                                                    \
            ok, si, offset, struct_t, printerr, NN_P23tfl_T1BME, UINT32, "NN_P23tfl_T1BME", "NN_P23tfl_T1BM");    \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, NN_P23tfl_T1B, UINT32, "NN_P23tfl_T1B", 10000); \
        EXT_STR_ITEM_INFO_ZZP(                                                                                    \
            ok, si, offset, struct_t, printerr, NN_P23tfl_T1Bv, UINT32, "NN_P23tfl_T1Bv", "NN_P23tfl_T1B");       \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, NN_P23tfl_T2BM, UINT32, "NN_P23tfl_T2BM", 50);  \
        EXT_STR_ITEM_INFO_ZZP(                                                                                    \
            ok, si, offset, struct_t, printerr, NN_P23tfl_T2BMI, UINT32, "NN_P23tfl_T2BMI", "NN_P23tfl_T2BM");    \
        EXT_STR_ITEM_INFO_ZZP(                                                                                    \
            ok, si, offset, struct_t, printerr, NN_P23tfl_T2BME, UINT32, "NN_P23tfl_T2BME", "NN_P23tfl_T2BM");    \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, NN_P23tfl_T2B, UINT32, "NN_P23tfl_T2B", 10000); \
        EXT_STR_ITEM_INFO_ZZP(                                                                                    \
            ok, si, offset, struct_t, printerr, NN_P23tfl_T2Bv, UINT32, "NN_P23tfl_T2Bv", "NN_P23tfl_T2B");       \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, NN_P23tct_T1BM, UINT32, "NN_P23tct_T1BM", 50);  \
        EXT_STR_ITEM_INFO_ZZP(                                                                                    \
            ok, si, offset, struct_t, printerr, NN_P23tct_T1BMI, UINT32, "NN_P23tct_T1BMI", "NN_P23tct_T1BM");    \
        EXT_STR_ITEM_INFO_ZZP(                                                                                    \
            ok, si, offset, struct_t, printerr, NN_P23tct_T1BME, UINT32, "NN_P23tct_T1BME", "NN_P23tct_T1BM");    \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, NN_P23tct_T1B, UINT32, "NN_P23tct_T1B", 10000); \
        EXT_STR_ITEM_INFO_ZZP(                                                                                    \
            ok, si, offset, struct_t, printerr, NN_P23tct_T1Bv, UINT32, "NN_P23tct_T1Bv", "NN_P23tct_T1B");       \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, NN_P23tct_T2BM, UINT32, "NN_P23tct_T2BM", 50);  \
        EXT_STR_ITEM_INFO_ZZP(                                                                                    \
            ok, si, offset, struct_t, printerr, NN_P23tct_T2BMI, UINT32, "NN_P23tct_T2BMI", "NN_P23tct_T2BM");    \
        EXT_STR_ITEM_INFO_ZZP(                                                                                    \
            ok, si, offset, struct_t, printerr, NN_P23tct_T2BME, UINT32, "NN_P23tct_T2BME", "NN_P23tct_T2BM");    \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, NN_P23tct_T2B, UINT32, "NN_P23tct_T2B", 10000); \
        EXT_STR_ITEM_INFO_ZZP(                                                                                    \
            ok, si, offset, struct_t, printerr, NN_P23tct_T2Bv, UINT32, "NN_P23tct_T2Bv", "NN_P23tct_T2B");       \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, NN_P23tft_T1BM, UINT32, "NN_P23tft_T1BM", 50);  \
        EXT_STR_ITEM_INFO_ZZP(                                                                                    \
            ok, si, offset, struct_t, printerr, NN_P23tft_T1BMI, UINT32, "NN_P23tft_T1BMI", "NN_P23tft_T1BM");    \
        EXT_STR_ITEM_INFO_ZZP(                                                                                    \
            ok, si, offset, struct_t, printerr, NN_P23tft_T1BME, UINT32, "NN_P23tft_T1BME", "NN_P23tft_T1BM");    \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, NN_P23tft_T1B, UINT32, "NN_P23tft_T1B", 10000); \
        EXT_STR_ITEM_INFO_ZZP(                                                                                    \
            ok, si, offset, struct_t, printerr, NN_P23tft_T1Bv, UINT32, "NN_P23tft_T1Bv", "NN_P23tft_T1B");       \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, NN_P23tft_T2BM, UINT32, "NN_P23tft_T2BM", 50);  \
        EXT_STR_ITEM_INFO_ZZP(                                                                                    \
            ok, si, offset, struct_t, printerr, NN_P23tft_T2BMI, UINT32, "NN_P23tft_T2BMI", "NN_P23tft_T2BM");    \
        EXT_STR_ITEM_INFO_ZZP(                                                                                    \
            ok, si, offset, struct_t, printerr, NN_P23tft_T2BME, UINT32, "NN_P23tft_T2BME", "NN_P23tft_T2BM");    \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, NN_P23tft_T2B, UINT32, "NN_P23tft_T2B", 10000); \
        EXT_STR_ITEM_INFO_ZZP(                                                                                    \
            ok, si, offset, struct_t, printerr, NN_P23tft_T2Bv, UINT32, "NN_P23tft_T2Bv", "NN_P23tft_T2B");       \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, NN_P24tcl_T1BM, UINT32, "NN_P24tcl_T1BM", 50);  \
        EXT_STR_ITEM_INFO_ZZP(                                                                                    \
            ok, si, offset, struct_t, printerr, NN_P24tcl_T1BMI, UINT32, "NN_P24tcl_T1BMI", "NN_P24tcl_T1BM");    \
        EXT_STR_ITEM_INFO_ZZP(                                                                                    \
            ok, si, offset, struct_t, printerr, NN_P24tcl_T1BME, UINT32, "NN_P24tcl_T1BME", "NN_P24tcl_T1BM");    \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, NN_P24tcl_T1B, UINT32, "NN_P24tcl_T1B", 10000); \
        EXT_STR_ITEM_INFO_ZZP(                                                                                    \
            ok, si, offset, struct_t, printerr, NN_P24tcl_T1Bv, UINT32, "NN_P24tcl_T1Bv", "NN_P24tcl_T1B");       \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, NN_P24tcl_T2BM, UINT32, "NN_P24tcl_T2BM", 50);  \
        EXT_STR_ITEM_INFO_ZZP(                                                                                    \
            ok, si, offset, struct_t, printerr, NN_P24tcl_T2BMI, UINT32, "NN_P24tcl_T2BMI", "NN_P24tcl_T2BM");    \
        EXT_STR_ITEM_INFO_ZZP(                                                                                    \
            ok, si, offset, struct_t, printerr, NN_P24tcl_T2BME, UINT32, "NN_P24tcl_T2BME", "NN_P24tcl_T2BM");    \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, NN_P24tcl_T2B, UINT32, "NN_P24tcl_T2B", 10000); \
        EXT_STR_ITEM_INFO_ZZP(                                                                                    \
            ok, si, offset, struct_t, printerr, NN_P24tcl_T2Bv, UINT32, "NN_P24tcl_T2Bv", "NN_P24tcl_T2B");       \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, NN_P24tfl_T1BM, UINT32, "NN_P24tfl_T1BM", 50);  \
        EXT_STR_ITEM_INFO_ZZP(                                                                                    \
            ok, si, offset, struct_t, printerr, NN_P24tfl_T1BMI, UINT32, "NN_P24tfl_T1BMI", "NN_P24tfl_T1BM");    \
        EXT_STR_ITEM_INFO_ZZP(                                                                                    \
            ok, si, offset, struct_t, printerr, NN_P24tfl_T1BME, UINT32, "NN_P24tfl_T1BME", "NN_P24tfl_T1BM");    \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, NN_P24tfl_T1B, UINT32, "NN_P24tfl_T1B", 10000); \
        EXT_STR_ITEM_INFO_ZZP(                                                                                    \
            ok, si, offset, struct_t, printerr, NN_P24tfl_T1Bv, UINT32, "NN_P24tfl_T1Bv", "NN_P24tfl_T1B");       \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, NN_P24tfl_T2BM, UINT32, "NN_P24tfl_T2BM", 50);  \
        EXT_STR_ITEM_INFO_ZZP(                                                                                    \
            ok, si, offset, struct_t, printerr, NN_P24tfl_T2BMI, UINT32, "NN_P24tfl_T2BMI", "NN_P24tfl_T2BM");    \
        EXT_STR_ITEM_INFO_ZZP(                                                                                    \
            ok, si, offset, struct_t, printerr, NN_P24tfl_T2BME, UINT32, "NN_P24tfl_T2BME", "NN_P24tfl_T2BM");    \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, NN_P24tfl_T2B, UINT32, "NN_P24tfl_T2B", 10000); \
        EXT_STR_ITEM_INFO_ZZP(                                                                                    \
            ok, si, offset, struct_t, printerr, NN_P24tfl_T2Bv, UINT32, "NN_P24tfl_T2Bv", "NN_P24tfl_T2B");       \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, NN_P24tct_T1BM, UINT32, "NN_P24tct_T1BM", 50);  \
        EXT_STR_ITEM_INFO_ZZP(                                                                                    \
            ok, si, offset, struct_t, printerr, NN_P24tct_T1BMI, UINT32, "NN_P24tct_T1BMI", "NN_P24tct_T1BM");    \
        EXT_STR_ITEM_INFO_ZZP(                                                                                    \
            ok, si, offset, struct_t, printerr, NN_P24tct_T1BME, UINT32, "NN_P24tct_T1BME", "NN_P24tct_T1BM");    \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, NN_P24tct_T1B, UINT32, "NN_P24tct_T1B", 10000); \
        EXT_STR_ITEM_INFO_ZZP(                                                                                    \
            ok, si, offset, struct_t, printerr, NN_P24tct_T1Bv, UINT32, "NN_P24tct_T1Bv", "NN_P24tct_T1B");       \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, NN_P24tct_T2BM, UINT32, "NN_P24tct_T2BM", 50);  \
        EXT_STR_ITEM_INFO_ZZP(                                                                                    \
            ok, si, offset, struct_t, printerr, NN_P24tct_T2BMI, UINT32, "NN_P24tct_T2BMI", "NN_P24tct_T2BM");    \
        EXT_STR_ITEM_INFO_ZZP(                                                                                    \
            ok, si, offset, struct_t, printerr, NN_P24tct_T2BME, UINT32, "NN_P24tct_T2BME", "NN_P24tct_T2BM");    \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, NN_P24tct_T2B, UINT32, "NN_P24tct_T2B", 10000); \
        EXT_STR_ITEM_INFO_ZZP(                                                                                    \
            ok, si, offset, struct_t, printerr, NN_P24tct_T2Bv, UINT32, "NN_P24tct_T2Bv", "NN_P24tct_T2B");       \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, NN_P24tft_T1BM, UINT32, "NN_P24tft_T1BM", 50);  \
        EXT_STR_ITEM_INFO_ZZP(                                                                                    \
            ok, si, offset, struct_t, printerr, NN_P24tft_T1BMI, UINT32, "NN_P24tft_T1BMI", "NN_P24tft_T1BM");    \
        EXT_STR_ITEM_INFO_ZZP(                                                                                    \
            ok, si, offset, struct_t, printerr, NN_P24tft_T1BME, UINT32, "NN_P24tft_T1BME", "NN_P24tft_T1BM");    \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, NN_P24tft_T1B, UINT32, "NN_P24tft_T1B", 10000); \
        EXT_STR_ITEM_INFO_ZZP(                                                                                    \
            ok, si, offset, struct_t, printerr, NN_P24tft_T1Bv, UINT32, "NN_P24tft_T1Bv", "NN_P24tft_T1B");       \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, NN_P24tft_T2BM, UINT32, "NN_P24tft_T2BM", 50);  \
        EXT_STR_ITEM_INFO_ZZP(                                                                                    \
            ok, si, offset, struct_t, printerr, NN_P24tft_T2BMI, UINT32, "NN_P24tft_T2BMI", "NN_P24tft_T2BM");    \
        EXT_STR_ITEM_INFO_ZZP(                                                                                    \
            ok, si, offset, struct_t, printerr, NN_P24tft_T2BME, UINT32, "NN_P24tft_T2BME", "NN_P24tft_T2BM");    \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, NN_P24tft_T2B, UINT32, "NN_P24tft_T2B", 10000); \
        EXT_STR_ITEM_INFO_ZZP(                                                                                    \
            ok, si, offset, struct_t, printerr, NN_P24tft_T2Bv, UINT32, "NN_P24tft_T2Bv", "NN_P24tft_T2B");       \
                                                                                                                  \
    } while (0);

#endif /*__GUARD_H101_RAW_NNP_TAMEX_EXT_H101_RAW_NNP_TAMEX_H__*/
