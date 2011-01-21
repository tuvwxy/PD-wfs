/******************************************************************************
 *  wfs~
 *
 *  Copyright (C) 2008  Toshiro Yamada
 *
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *
 ******************************************************************************/
#ifndef _WFS_
#define _WFS_

#ifdef _MAXMSP_
    #include "ext.h"
    #include "ext_obex.h"
    #include "z_dsp.h"
#else
    #define _PD_ 1
    #include "m_pd.h"
#endif

#define VERSION "0.10"
#define COPYRIGHT "Copyright(c) 2008, Toshiro Yamada"

#define MAXSPEAKER (24)         /* the maximum # of speakers in an array */
#define DELAYTABLESIZE (65536)  /* delay tablesize, must be power of 2 */
#define EXTRADELAY (2048)       /* extra delay buffer */
static t_class *wfs_tilde_class;    /* Global pointer to this class */

/************************ source *********************************************/
typedef struct _source {
    float sr_loc[3];            /* source location (x,y,z) (point src only) */
    float sr_back;              /* directionality (point source only) */
    float sr_angle;             /* source direction */
    float sr_amp;               /* amplitude of source (plane wave only */

    float *sr_px;               /* pointers to vectors from inlets */
    float *sr_py;
    float *sr_pz; 
    float *sr_pangle;
    float *sr_pback;
    float *sr_pamp;
} t_source;

/************************ loudspeaker ****************************************/
typedef struct _speaker {
    float  sp_loc[3];           /* speaker location (x,y,z)*/
    float  sp_vec[3];           /* vector from source to speaker */
    double sp_dist;             /* distance between each spkr & notional src */
    double sp_gain;             /* attenuation factor */
    double sp_delay;            /* delay time in samples */
    double sp_phi;              /* angle to the source */
    double sp_radVec;           /* radiation vector */
} t_speaker;

/************************* WFS~ **********************************************/
typedef struct _wfs_tilde {
    t_object x_obj;             /* pd routine */
    
    float x_afreq;              /* aliasing frequency */
    float x_zref;               /* reference line */
    float *x_pzref;             /* pointer to inlet stream */
    int x_stype;                /* source type; 0: point, 1: plane */
    int x_scard;                /* Omni or cardioid directionality */
    int x_swindow;              /* window on/off */    
    int x_sdim;                 /* dimension; 0: 2D, 1: 3D (with y-axis) */
    
    t_sample *x_in;             /* input buffer */
    t_source *x_src;            /* notional (point) source */
    
    /* loudspeakers */
    int sp_n;                   /* number of loudspeakers */
    float sp_dx;                /* spacing between the loudspeakers */
    t_speaker *sp[MAXSPEAKER];  /* pointers to loudspeakers in the array */
    
    int sp_windowCount;         /* number of speaker on the edge to weight */
    double sp_windowAlpha;      /* window param 1 */
    double sp_windowBeta;       /* window param 2 */
    double sp_window[MAXSPEAKER];/* window to minimize truncation effect */
    
    /* Delay Line */
    float *delay_vec;           /* pointer to delay buffer */
    int delay_size;             /* buffer size */
    int delay_ptr;              /* phase for the write pointer */
    
    t_sample *x_out[MAXSPEAKER];/* output buffers for each speaker */
    
    float x_f;
} t_wfs_tilde;

/************************ Prototypes *****************************************/
/************************ various calculations *******************************/
static double wfs_gain (t_wfs_tilde *x, int i);
static double wfs_distance (float *a, float *b, int length);
static double wfs_dist2pow (float *a, float *b, int length);
static double wfs_dotProduct (float *a, float *b, int length);
static void wfs_vecAtoB (float *ab, float *a, float *b, int length);
static void wfs_vecCopy (float *a, float *b, int length);
static double wfs_hermite (double frac_pos, double xm1, double x0, double x1, double x2);

//static void wfs_addmethod ();
/************************ set methods ****************************************/
static void wfs_setWindow (t_wfs_tilde *x);
static void wfs_setWindowOpt (t_wfs_tilde *x, float on);
static void wfs_setWindowParam (t_wfs_tilde *x, float halfCount, float alpha, float beta);

/************************ print methods **************************************/
static void wfs_printAll (t_wfs_tilde *x);
static void wfs_printCopyright (t_wfs_tilde *x);
static void wfs_printGain (t_wfs_tilde *x); 
static void wfs_printDelay (t_wfs_tilde *x);
static void wfs_printDistance (t_wfs_tilde *x);
static void wfs_printHelp (t_wfs_tilde *x);
static void wfs_printInit (t_wfs_tilde *x);
static void wfs_printRadVec (t_wfs_tilde *x);
static void wfs_printReference (t_wfs_tilde *x);
static void wfs_printSpeaker (t_wfs_tilde *x);
static void wfs_printSpeakerVec (t_wfs_tilde *x);
static void wfs_printSource (t_wfs_tilde *x);
static void wfs_printVersion (t_wfs_tilde *x);
static void wfs_printWindow (t_wfs_tilde *x);

#endif /* _WFS_ */

