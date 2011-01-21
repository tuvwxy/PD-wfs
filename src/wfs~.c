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
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "wfs~.h" 

/* Macro for Max/MSP and Pd Compatibility */
#ifdef _MAXMSP_
#define METHODCLASS
#define METHOD method
#define SYM(x) (x)
#define GETATOMFLOAT(a,b) ((a)+(b))->a_w.w_float
#define GETATOMLONG(a,b) ((a)+(b))->a_w.w_long
#else /* Pure Data */
#define METHODCLASS wfs_tilde_class,
#define METHOD t_method
#define SYM(x) gensym(x)
#define GETATOMFLOAT(a,b) atom_getfloatarg((b),argc,(a))
#define GETATOMLONG(a,b) atom_getfloatarg((b),argc,(a))
#endif

#define SQR(z) ((z)*(z))    /* calculate square */

#define ATTEN_SPATIAL 0
#define ATTEN_LOUDSPEAKER 0
#define ATTEN_COMPENSATION 0

static double const pi = 3.14159265359;
static double const two_pi = 6.28318530718;
static double const C = 343;            /* velocity of sound */
static double const CINV = 0.00294117;  /* inverse of velocity of sound */
static double const D90 = 1.57079632;   /* 90 degress pi/2 */

static int const POINTSOURCE = 0;   /* notional source is a point source */
static int const PLANEWAVE   = 1;   /* notional source is a plane wave */
static int const OMNI        = 0;   /* same as back = 1 */
static int const CARD        = 1;   /* use Moore's radiation vec algorithm */
static int const WINDOW_OFF  = 0;   /* don't apply spatial window */
static int const WINDOW_ON   = 1;   /* apply window */
static int const TWOD        = 0;   /* 2 Dimension: X and Z coordinates only */
static int const THREED      = 1;   /* 3 Dimension: X, Y and Z coordinates */
static int const X           = 0;
static int const Y           = 1;
static int const Z           = 2;
static char const USAGE[] = "wfs~: Options: [2d/3d] [-loc X (Y) Z] [-n #] [-dx #] "
    "[point/plane] [card/omni] [window/nowindow]";

static int initBuffer = 0;      /* needs to buffer EXTRADELAY for delay */

/* points to class_addclass and addmess functions */
#ifdef _MAXMSP_
    static void (*func_addmess)(method f, char *s, short type, ...) = &addmess;  
#else
    static void (*func_addmess)(t_class *c, t_method fn, t_symbol *sel,
                                t_atomtype arg1, ...) = &class_addmethod;
#endif

/* BUGS:
 *  - In plane wave mode, crashes when the angle value exceeds certain values */
/**************************** wfs_tilde_new *********************************/
static void *wfs_tilde_new(t_symbol *s, int argc, t_atom *argv) 
{
#ifdef _MAXMSP_
    t_wfs_tilde *x = (t_wfs_tilde *) newobject(wfs_tilde_class);
#else
    t_wfs_tilde *x = (t_wfs_tilde *) pd_new(wfs_tilde_class);
#endif
    t_symbol *firstarg;
    
    int i = 0; /* to iterate each loudspeaker */

    /*----------------------*
     * Initialize System
     *----------------------*/
    /* loudspeaker coordinates */
    float xloc = 0;
    float yloc = 0;
    float zloc = 0;
    /* Default settings */
    x->sp_n         = 8;            /* 8 loudspeakers */
    x->sp_dx        = 0.1016;       /* m - 4" speaker spacing */
    x->x_stype      = POINTSOURCE;  /* Point source */
    x->x_scard      = OMNI;         /* Omni direction source */
    x->x_swindow    = WINDOW_OFF;   /* No weidghting on speakers */
    x->x_sdim       = TWOD;         /* 2D horizontal array */
    x->x_afreq      = 0;            /* Spatial aliasing frequency */ 
    x->x_zref       = 3;            /* m - Reference line at z = 3m */

    x->x_f          = 0;
    
    /* Custom window setup */
    x->sp_windowCount = 4;          /* Window 4 speakers on each edge */
    x->sp_windowAlpha = 1;          /* Param for window: Max possible spacing */
    x->sp_windowBeta  = 0;          /* Param for window: No margin from edge */
    
    x->delay_size   = DELAYTABLESIZE;
    x->delay_vec    = (float *) calloc(x->delay_size * sizeof(float), 1);
    x->delay_ptr    = 0;
    
    
    /*----------------------*
     * Parse Arguments
     *----------------------*/
    while (argc > 0) {
#if _MAXMSP_
        firstarg = argv->a_w.w_sym;
#else
        firstarg = atom_getsymbolarg(0, argc, argv);
#endif

        if (!strcmp(firstarg->s_name, "-loc") 
                && x->x_sdim == TWOD && argc > 2) {            
            xloc = GETATOMFLOAT(argv,1);
            zloc = GETATOMFLOAT(argv,2);
            yloc = 0;
            argc -= 3;
            argv += 3;
        }
        else if (!strcmp(firstarg->s_name, "-loc")
                    && x->x_sdim == THREED && argc > 3) {
            xloc = GETATOMFLOAT(argv,1);
            yloc = GETATOMFLOAT(argv,2);
            zloc = GETATOMFLOAT(argv,3);
            argc -= 4;
            argv += 4;            
        }
        else if (!strcmp(firstarg->s_name, "-n") && argc > 1) {
            x->sp_n = GETATOMLONG(argv,1);
            /* check size limit */
            if (x->sp_n > MAXSPEAKER) {
                post("wfs~: error: %d loudspeakers are allowed "
                         "in a single array", MAXSPEAKER);
                return ((void *) NULL);
            }
            argc -= 2;
            argv += 2;
        }
        else if (!strcmp(firstarg->s_name, "-dx") && argc > 1) {
            x->sp_dx = GETATOMFLOAT(argv,1);
            argc -= 2;
            argv += 2;
        }
        else if (!strcmp(firstarg->s_name, "point")) {
            x->x_stype = POINTSOURCE;
            argc -= 1;
            argv += 1;
        }
        else if (!strcmp(firstarg->s_name, "plane")) {
            x->x_stype = PLANEWAVE;
            argc -= 1;
            argv += 1;
        }
        else if (!strcmp(firstarg->s_name, "omni")) {
            x->x_scard = OMNI;
            argc -= 1;
            argv += 1;
        }
        else if (!strcmp(firstarg->s_name, "card")) {
            x->x_scard = CARD;
            argc -= 1;
            argv += 1;
        }
        else if (!strcmp(firstarg->s_name, "window")) {
            x->x_swindow = WINDOW_ON;
            argc -= 1;
            argv += 1;	
        }
        else if (!strcmp(firstarg->s_name, "nowindow")) {
            x->x_swindow = WINDOW_OFF;
            argc -= 1;
            argv += 1;

        }
        else if (!strcmp(firstarg->s_name, "2D") 
                    || !strcmp(firstarg->s_name, "2d")) {
            x->x_sdim = TWOD;
            argc -= 1;
            argv += 1;
        }
        else if (!strcmp(firstarg->s_name, "3D")
                    || !strcmp(firstarg->s_name, "3d")) {
            x->x_sdim = THREED;
            argc -= 1;
            argv += 1;
        }
        else {
            post("wfs~: %s: unknown flag or argument missing",
                    firstarg->s_name);
            // post("wfs~: error: unknown flag or argument missing");
            post(USAGE);

            return ((void *) NULL);
            argc -= 1;
            argv += 1;
        }
    }

    if (x->x_scard && x->x_stype == PLANEWAVE)
        post("wfs~: warning: Plane waves do not have card option.");
    
    /*----------------------*
     * Initialize each loudspeakers 
     *----------------------*/
    for (i = 0; i < x->sp_n; i++) { 
        x->sp[i] = (t_speaker *) calloc(sizeof(t_speaker), 1);
        x->sp[i]->sp_loc[X] = (i - (x->sp_n - 1)/2.0) * x->sp_dx; 
        x->sp[i]->sp_loc[Y] = yloc;
        x->sp[i]->sp_loc[Z] = zloc;
                
        x->sp[i]->sp_gain  = 0;
        x->sp[i]->sp_dist  = 0;
        x->sp[i]->sp_delay = EXTRADELAY;

#ifdef _MAXMSP_
        outlet_new((t_object *)x, "signal");
#else
        outlet_new(&x->x_obj, &s_signal);
#endif
    }

    if (x->x_swindow) {
        wfs_setWindow(x);
    }

    /*----------------------*
     * Initialize source 
     *----------------------*/
    x->x_src = (t_source *) calloc(sizeof(t_source), 1);
    x->x_src->sr_loc[X] = 0;
    x->x_src->sr_loc[Y] = 0;
    x->x_src->sr_loc[Z] = 0;
    /* sr_angle: value between -pi to pi
     * 0 points negative z-axis (toward the lisnter) */
    x->x_src->sr_angle  = 0;
    /* sr_back: value between 0 to 1
     * 1 == omnidirectional */
    x->x_src->sr_back   = 1;
    x->x_src->sr_amp    = 1;
    
    x->x_afreq          = C /(2.0*x->sp_dx);   /* worse aliasing frequency */
 
    wfs_printInit(x);

    /*----------------------*
     * Create inlets 
     *----------------------*/
    /* Point source: (signal, X, Y, Z, angle, back) */
    if (x->x_stype == POINTSOURCE) { 
#ifdef _MAXMSP_
        if (x->x_sdim == THREED && x->x_scard)
            dsp_setup((t_pxobject *)x, 6);
        else if (x->x_scard)
            dsp_setup((t_pxobject *)x, 5);
        else if (x->x_sdim == THREED)
            dsp_setup((t_pxobject *)x, 4);
        else
            dsp_setup((t_pxobject *)x, 3);            
#else            
        /* X */
        inlet_new(&x->x_obj, &x->x_obj.ob_pd, &s_signal, &s_signal);
        /* Y */
        if (x->x_sdim == THREED)
            inlet_new(&x->x_obj, &x->x_obj.ob_pd, &s_signal, &s_signal);
        /* Z */
        inlet_new(&x->x_obj, &x->x_obj.ob_pd, &s_signal, &s_signal); 
        /* angle and back */
        if (x->x_scard) {
            inlet_new(&x->x_obj, &x->x_obj.ob_pd, &s_signal, &s_signal);
            inlet_new(&x->x_obj, &x->x_obj.ob_pd, &s_signal, &s_signal); 
        }
#endif
    }
    /* Plane wave: (signal, angle, amplitude)*/
    else  {
#ifdef _MAXMSP_
        dsp_setup((t_pxobject *)x, 3);
#else
        inlet_new(&x->x_obj, &x->x_obj.ob_pd, &s_signal, &s_signal);
        inlet_new(&x->x_obj, &x->x_obj.ob_pd, &s_signal, &s_signal);
#endif
    }

    return (void *)x;
}

/**************************** wfs_tilde_perform *******************************/
static t_int *wfs_tilde_perform (t_int *w) 
{
    t_wfs_tilde *x = (t_wfs_tilde *) (w[1]);
    int blocksize = (int)(w[2]);    

    int i = 0;         /* to i each loudspeaker */
    int n = 0;         /* to iterate every sample in a block */
    
    t_sample *in = x->x_in;
    t_sample *out[MAXSPEAKER];
    
    float sr_loc[] = {0,0,0};   /* notional source coordinate (x,y,z) */
    float sr_angle = 0;         /* direction/angle of notional source */
    float sr_back  = 0;         /* directionality of notinoal source */
    float sr_amp   = 1;         /* overall amplitude of plane wave */

    double phi = 0;             /* angle between the source and speaker */
    float dx = x->sp_dx;        /* loudspeaker spacing */

    float *vp = x->delay_vec;       /* pointer to the delay table */
    int tableSize = x->delay_size;  /* table size of the delay line */
    int tableMask = tableSize-1;    /* mask i to table boundary */
    int readHead;                   /* read i of delay table */
    int writeHead = x->delay_ptr;   /* write i of delay table */
    int iDelay = 0;                 /* integer value of delay time */
    double frac = 0;                /* fraction of delay time */

    double newDist;                 /* new distance */
    double delay[MAXSPEAKER];       /* delay time */
    double gain[MAXSPEAKER];        /* gain values for each speaker*/
    double radVec[MAXSPEAKER];      /* radiation vector for each speaker */
    double delta_delay[MAXSPEAKER]; /* delta change between blocks */
    double delta_gain[MAXSPEAKER];  /* delta change between blocks */
    double delta_radVec[MAXSPEAKER];/* delta change between blocks */


    /*---------------------*
     *  Values from Inlets
     *---------------------*/
    /* Isotropic Wave */
    if (x->x_stype == POINTSOURCE) {
        sr_loc[X] = *x->x_src->sr_px;
        sr_loc[Z] = *x->x_src->sr_pz;
        
        if (x->x_sdim == THREED)
            sr_loc[Y] = *x->x_src->sr_py;
        
        if (x->x_scard) {
            sr_angle = *x->x_src->sr_pangle;
            sr_back  = *x->x_src->sr_pback;
    
            /* Wrap angle to (0, 1] */
            while (sr_angle >= 1)
                sr_angle -= 1;
            while (sr_angle < 0)
                sr_angle += 1;
                
            /* bound back between 0 to 1 */
            if (sr_back > 1)
                sr_back = 1;
            if (sr_back < 0)
                sr_back = 0;
        }
    }
    
    /* TODO: Plane Wave */
    else {
        sr_angle = *x->x_src->sr_pangle;
        sr_amp   = *x->x_src->sr_pamp;
        
        sr_angle -= 0.5;
        
        /* Wrap angle to (0, 1) */
        while (sr_angle > 1)
            sr_angle -= 1;
        while (sr_angle < 0)
            sr_angle += 1;
            
        /* Set angle boundary between -0.25 to 0.25 */
        sr_angle = (sr_angle + 0.5) * 0.5;
        //sr_angle *= 0.5;
    }
        
    /*---------------------*
     *  WFS Setup Routine
     *---------------------*/
    for (i = 0; i < x->sp_n; i++) {
        out[i]   = x->x_out[i];
        gain[i]  = x->sp[i]->sp_gain;
        delay[i] = x->sp[i]->sp_delay;
        radVec[i]= x->sp[i]->sp_radVec;
        
        /*** Point Source ***/
        if (x->x_stype == POINTSOURCE) {            
            /* Calculate distance, gain, and delay */
            newDist = wfs_distance(x->sp[i]->sp_loc, sr_loc, 3);
            x->sp[i]->sp_dist = newDist;
            x->sp[i]->sp_gain = wfs_gain(x, i);
            x->sp[i]->sp_delay = x->sp[i]->sp_dist * CINV * sys_getsr();

            /*** EXPERIMENTAL: Calcuate Radiation Vector ****/
            if (x->x_sdim == TWOD && x->x_scard) {
                float vecA[3];
                float vecB[] = {0, 0, 1};
                
                if (sr_back == 1 && x->x_scard) {
                    /* Find a vector from speaker to source */
                    wfs_vecAtoB(vecA, x->sp[i]->sp_loc, sr_loc, 3);

                    /* Find the angle to the source
                     * cos(phi) = A.B / |A||B| 
                     *     phi should be between 0 to 2pi */
                    phi = acos( wfs_dotProduct(vecA, vecB, 3) /
                            sqrt(wfs_dist2pow(vecA, NULL, 3) * wfs_dist2pow(vecB, NULL, 3))
                        );

                    /* Set angle to betwee -pi and pi */
                    phi -= pi; 

                    /* Store angle and vector values */
                    x->sp[i]->sp_phi = phi;
                    wfs_vecCopy(vecA, x->sp[i]->sp_vec, 3);

                    x->sp[i]->sp_radVec = SQR(1 + (sr_back - 1) * fabs(sr_angle - phi) / pi);
                }
                else
                    x->sp[i]->sp_radVec = 1;
            }
        }
        
        /*** TODO: Plane Wave ***/
        else {  
            /* Calculate Distance */
            if (sr_angle < 0)
                newDist = dx * sin(-sr_angle*two_pi) * (x->sp_n-1 - i);
            else
                newDist = dx * sin(sr_angle*two_pi) * i;
            
            /* add 1 to distance to avoid near field amplitude stuff */
            newDist += 1.0; 
                            
            /* Calculate gain and delay */
            x->sp[i]->sp_dist = newDist;
            /*x->sp[i]->sp_gain = sr_amp /sqrt(1 + newDist); */
            x->sp[i]->sp_gain = sr_amp /(1 + newDist);
            x->sp[i]->sp_delay = newDist * CINV * sys_getsr();

#if 0
            /* Plane wave propagates left to right */
            if (sr_angle < 0) {         
                x->sp[i]->sp_dist = newDist;
                x->sp[i]->sp_gain = sr_amp /sqrt(1 + newDist);
                x->sp[i]->sp_delay = newDist * CINV * sys_getsr();
            }
            /* Plane wave propagates right to left */
            else if (sr_angle > 0) {    
                x->sp[x->sp_n-i]->sp_dist = newDist;
                x->sp[x->sp_n-i]->sp_gain = sr_amp /sqrt(1 + newDist);
                x->sp[x->sp_n-i]->sp_delay = newDist * CINV * sys_getsr();
            }
            /* Plane wave propagates straight ahead */
            else {
                x->sp[i]->sp_dist = 1;
                x->sp[i]->sp_gain = sr_amp;
                x->sp[i]->sp_delay = 0;
            }
#endif
        }

        /*** interpolation routine ***/
        delta_gain[i]  = (x->sp[i]->sp_gain - gain[i]) / blocksize;
        delta_delay[i] = (x->sp[i]->sp_delay - delay[i]) / blocksize;
        delta_radVec[i]= (x->sp[i]->sp_radVec - radVec[i]) / blocksize;
    }

    /*---------------------*
     * Block Processing
     *---------------------*/
    for (n = 0; n < blocksize; n++) {   /* for each sample... */
        
        vp[writeHead] = *in++;
        
        /*** Process for each speaker ***/
        for (i = 0; i < x->sp_n; i++) {
            /* Variables for polynomial interpolation */
            double x2mx1, x0, x1, x2, x3;

            gain[i]   += delta_gain[i];
            delay[i]  += delta_delay[i];
            radVec[i] += delta_radVec[i];

            /**** Calculate Delay ****/
            iDelay = (int) floor(delay[i]);
            frac = delay[i] - (double)iDelay;
            
            /* TODO: this should really check if the source is infront of 
            the array not simply by comparing the z axis... this only works
            for arrays along x-axis */
            if (x->x_src->sr_loc[Z] > x->sp[0]->sp_loc[Z])
                iDelay = -iDelay;

            readHead = writeHead - (iDelay + EXTRADELAY);

#if 0 /* Hermite Spline interpolation - doesn't sound as good as Miller's */
            double xm1, x0, x1, x2; /* 4 samples around the delay sample */
            xm1 = vp[(readHead-1) & tableMask];
            x0  = vp[(readHead-0) & tableMask];
            x1  = vp[(readHead+1) & tableMask];
            x2  = vp[(readHead+2) & tableMask];
            
            *out[i]++ = x->sp_window[i] * gain[i]
                         * wfs_hermite(frac, xm1, x0, x1, x2);
#endif

#if 1 /* Miller Puckette's interpolation funciont in vd~ object */
            x3 = vp[(readHead-3) & tableMask];
            x2 = vp[(readHead-2) & tableMask];
            x1 = vp[(readHead-1) & tableMask];
            x0 = vp[(readHead-0) & tableMask];
            x2mx1 = x2-x1;
            *out[i] = ( x1 + frac * (
               x2mx1 - 0.1666667f * (1.-frac) * (
                 (x3 - x0 - 3.0f*x2mx1) * frac + (x3 + 2.0f*x0 - 3.0f*x1)
               )
            ));
#endif
            /* Apply Radiation Vector */
            if (x->x_scard && x->x_sdim == TWOD)
                *out[i] *= radVec[i];
            /* Apply window (Weighting) function */
            if (x->x_swindow)
                *out[i] *= x->sp_window[i];
            /* Apply gain function */
            *out[i] *= gain[i];
            
            *out[i]++;
        }
        
        writeHead++;
        
        /* wrap write pointer within the table boundary */
        while (writeHead >= tableSize)
            writeHead -= tableSize;
    }

    x->delay_ptr = writeHead; 

    /* Save delay for the block processing */
    for (i = 0; i < x->sp_n; i++) {
        x->sp[i]->sp_gain   = gain[i];
        x->sp[i]->sp_delay  = delay[i];
        x->sp[i]->sp_radVec = radVec[i];
    }
    
    /* Save source information */
    if (x->x_stype == POINTSOURCE) { 
        x->x_src->sr_loc[X] = sr_loc[X];
        if (x->x_sdim == THREED)
            x->x_src->sr_loc[Y] = sr_loc[Y];
        x->x_src->sr_loc[Z] = sr_loc[Z];
        if (x->x_scard) {
            x->x_src->sr_angle = sr_angle;
            x->x_src->sr_back = sr_back; 
        }
    }
    else {  /* Plane Wave */
        x->x_src->sr_angle = sr_angle;
        x->x_src->sr_amp = sr_amp;
    }

    return (w+3);
}

/**************************** wfs_tilde_dsp *********************************/
#ifdef _MAXMSP_
void wfs_tilde_dsp(t_wfs_tilde *x, t_signal **sp, short *count)
#else
void wfs_tilde_dsp(t_wfs_tilde *x, t_signal **sp)
#endif
{
    int i, j = 0;
    
    /* signal inlet */
    x->x_in = (t_sample *) sp[j++]->s_vec;
    
    /* source informations */
    if (x->x_stype == POINTSOURCE) {
        x->x_src->sr_px = (t_sample *) sp[j++]->s_vec;
        x->x_src->sr_pz = (t_sample *) sp[j++]->s_vec;
        if (x->x_sdim == THREED) 
            x->x_src->sr_py = (float *) sp[j++]->s_vec;
        if (x->x_scard) {
            x->x_src->sr_pangle = (float *) sp[j++]->s_vec;
            x->x_src->sr_pback = (float *) sp[j++]->s_vec;
        }
    }
    else { /* plane wave */
        x->x_src->sr_pangle = (float *) sp[j++]->s_vec;
        x->x_src->sr_pamp = (float *) sp[j++]->s_vec;
    }

    /* outlets */
    for (i = 0; i < x->sp_n; i++) { 
        x->x_out[i] = (t_sample *) sp[j++]->s_vec;
    }
    
    dsp_add(wfs_tilde_perform, 2, x, sp[0]->s_n);
}

/**************************** wfs_tilde_free ********************************/
static void wfs_tilde_free (t_wfs_tilde *x)
{
    int i; /* free memory for each loudspeaker */

    for (i = 0; i < x->sp_n; i++) {
        freebytes(&x->sp[i], sizeof(t_speaker));
        x->sp[i] = NULL;
    }
    
    freebytes(&x->delay_vec, x->delay_size * sizeof(float));
    x->delay_vec = NULL;

    freebytes(&x->x_src, sizeof(t_source));
    x->x_src = NULL;
}

/**************************** wfs_tilde_setup *******************************/
#ifdef _MAXMSP_
void main(void)
#else
void wfs_tilde_setup(void) 
#endif 
{    
#ifdef _MAXMSP_
    setup((t_messlist **)&wfs_tilde_class, 
            (method)wfs_tilde_new, (method)wfs_tilde_free, 
            (short)sizeof(t_wfs_tilde), 0L, A_GIMME, 0);
#else
    wfs_tilde_class = class_new(gensym("wfs~"), 
            (t_newmethod)wfs_tilde_new, (t_method)wfs_tilde_free,
            sizeof(t_wfs_tilde), 0, A_GIMME, 0);

    CLASS_MAINSIGNALIN(wfs_tilde_class, t_wfs_tilde, x_f);
#endif

    /*
    addmethod((METHODCLASS METHOD)wfs_tilde_dsp, SYM("dsp"), 0);
    addmethod((METHODCLASS METHOD)wfs_setWindowParam, 
               SYM("setWindow"), A_FLOAT, A_FLOAT, A_FLOAT, 0);
    addmethod((METHODCLASS METHOD)wfs_setWindowOpt, 
               SYM("window"), A_FLOAT, 0);
    addmethod((METHODCLASS METHOD)wfs_printInit, SYM("print"), 0); */
    /*
    class_addmethod(wfs_tilde_class, (METHOD)wfs_tilde_dsp,
        SYM("dsp"), 0); */
    func_addmess(METHODCLASS (METHOD)wfs_tilde_dsp,
            SYM("dsp"), 0);
    /* set methods  */
    func_addmess(METHODCLASS (METHOD)wfs_setWindowOpt,
            SYM("window"), A_FLOAT, 0); 
    func_addmess(METHODCLASS (METHOD)wfs_setWindowParam,
            SYM("setWindow"), A_FLOAT, A_FLOAT, A_FLOAT, 0);
    /* print methods */
    func_addmess(METHODCLASS (METHOD)wfs_printInit,
            SYM("print"), 0);
    func_addmess(METHODCLASS (METHOD)wfs_printAll,
            SYM("printAll"), 0);
    func_addmess(METHODCLASS (METHOD)wfs_printGain,
            SYM("printGain"), 0);
    func_addmess(METHODCLASS (METHOD)wfs_printDelay,
            SYM("printDelay"), 0);
    func_addmess(METHODCLASS (METHOD)wfs_printDistance,
            SYM("printDistance"), 0);
    func_addmess(METHODCLASS (METHOD)wfs_printHelp,
            SYM("printHelp"), 0);
    func_addmess(METHODCLASS (METHOD)wfs_printRadVec,
            SYM("printRadVec"), 0);
    func_addmess(METHODCLASS (METHOD)wfs_printReference,
            SYM("printReference"), 0);
    func_addmess(METHODCLASS (METHOD)wfs_printSpeaker,
            SYM("printSpeaker"), 0);
    func_addmess(METHODCLASS (METHOD)wfs_printSpeakerVec,
            SYM("printSpeakerVec"), 0);
    func_addmess(METHODCLASS (METHOD)wfs_printSource,
            SYM("printSource"), 0);
    func_addmess(METHODCLASS (METHOD)wfs_printVersion,
            SYM("printVersion"), 0);
    func_addmess(METHODCLASS (METHOD)wfs_printWindow,
            SYM("printWindow"), 0);
            
#ifdef _MAXMSP_
    /* rescopy ('STR#', ResourceID);   Copy the assistance messages resource into
                                       Max's temp file */
    finder_addclass("All Objects","wfs~");     /* add class to the New object list */
#endif
}

/************************ wfs_gain ******************************************/
static double wfs_gain (t_wfs_tilde *x, int i) 
{
    double A = 1.0;
    double d = fabs(x->sp[i]->sp_dist);
    double r0 = fabs(x->x_zref - x->x_src->sr_loc[Z]);
    double r1 = fabs(x->x_zref - x->sp[i]->sp_loc[Z]);

    /* TODO:
    - It seems like the first part of the loudspeaker characteristic factor
      changes the dynamic too much.  When the speakers are spaced, the
      dynamic fluctuates. The gain at far distance sounds more
      realistic with this, though. */
    
#if 1 /* Inverse distance law of sound pressure = 1/d */
    if (x->x_stype == POINTSOURCE) {
        if (d <= 1) {
            A *= 1;
        }
        else if (x->x_src->sr_loc[Z] > x->sp[0]->sp_loc[Z]) {/* source in front of array */
            A *= (2 - 1/d);
        }
        else {  /* source behind the array */
            A /= d;
        }
    }
#endif

#if 0 /* Compensation factor = sqrt(|z - z1|/|z - z0|) */
    if (r0 != 0 && r1 != 0)
        A *= sqrt(r1/r0);
#endif

#if 0
    /* loudspeaker characteristic = cos(angle)/sqrt(2pi)
     *                 cos(angle) = |z0 - z1|/|rm - rn| */
    if (d != 0) {
        A *= (fabs(x->x_src->sr_loc[2] - x->sp[i]->sp_loc[2])/d) * pow(two_pi, -0.5); 
        A *= pow(two_pi, -0.5);
    }
#endif

    return A;
}

/************************ wfs_distance **************************************/
static double wfs_distance (float *a, float *b, int length)
{
    return sqrt(wfs_dist2pow(a, b, length));
}

/************************ wfs_dist2pow **************************************/
static double wfs_dist2pow (float *a, float *b, int length)
{
    int i;
    double retVal = 0;
    
	if (a == NULL) {
        float newA[length];
        a = newA;
        for (i = 0; i < length; i++)
            a[i] = 0;
    }
     
    if (b == NULL) {
        float newB[length];
        b = newB;
        for (i = 0; i < length; i++)
            b[i] = 0;
    }

    for (i = 0; i < length; i++)
        retVal += SQR(a[i] - b[i]);

    return retVal;
}

/************************ wfs_dotProduct ************************************/
static double wfs_dotProduct (float *a, float *b, int length)
{
    int i;
    double retVal = 0;

    if (a == NULL || b == NULL)
        return 0.0;

    for (i = 0; i < length; i++)
        retVal += a[i] * b[i];

    return retVal;
}

/************************ wfs_distance **************************************/
static void wfs_vecAtoB (float *ab, float *a, float *b, int length)
{
    int i;

    for (i = 0; i < length; i++)
       ab[i] = a[i] - b[i];
}

/************************ wfs_vecCopy ***************************************/
static void wfs_vecCopy (float *a, float *b, int length)
{
    int i;
    
    for (i = 0; i < length; i++)
        b[i] = a[i];
}

/************************ wfs_hermite ***************************************/
/* from http://www.musicdsp.org/archive.php?classid=5#93
 * by laurent de soras */
static double wfs_hermite
    (double frac_pos, double xm1, double x0, double x1, double x2) 
{
   const double    c     = (x1 - xm1) * 0.5f;
   const double    v     = x0 - x1;
   const double    w     = c + v;
   const double    a     = w + v + (x2 - x0) * 0.5f;
   const double    b_neg = w + a;

   return ((((a * frac_pos) - b_neg) * frac_pos + c) * frac_pos + x0);
}

/************************ set methods ***************************************/
/************************ wfs_setWindow *************************************/
void wfs_setWindow(t_wfs_tilde *x)
{
    int i = 0;      /* iterate loop */
    
    int N = x->sp_windowCount;
    int L = x->sp_n;
    double a, b;
    
    if (N > 0) {
        /* spacing between the speakers */
        a = x->sp_windowAlpha * pi/(double)(N*(N-1)) + pi/(double)N;
        /* extra margin on the edge of window */
        b = x->sp_windowBeta * (pi - a*(N-1));
    
        for (i = 0; i < L - N; i++) {
            if (i < N) {    /* Edge of array */
                x->sp_window[i] = x->sp_window[L-i-1] = 0.53836 - 0.4164 * cos(a*i + b);
            }
            else {          /* Middle */
                x->sp_window[i] = 1;
            }
        }
    }
}

/************************ wfs_setWindowOpt **********************************/
static void wfs_setWindowOpt (t_wfs_tilde *x, float on)
{
    if (on)
        x->x_swindow = WINDOW_ON;
    else
        x->x_swindow = WINDOW_OFF;
}

/************************ wfs_setWindowParam ********************************/
static void wfs_setWindowParam (t_wfs_tilde *x, float halfCount, float alpha, float beta)
{
    if (halfCount > 0) {
        
        x->sp_windowCount = (int)halfCount;
        
        /* Fix alpha and beta params to 0 to 1 */
        if (alpha > 1) {
            post("wfs~: setWindowParam: Warning - alpha fixed to 1 "
                 "(should be between 0 to 1)");
            alpha = 1;
        }
        else if (alpha < 0) {
            post("wfs~: setWindowParam: Warning - alpha fixed to 0 "
                 "(should be between 0 to 1)");
            alpha = 0;
        }
        if (beta > 1) {
            post("wfs~: setWindowParam: Warning - beta fixed to 1 "
                 "(should be between 0 to 1)");
            beta = 1;
        }
        else if (beta < 0) {
            post("wfs~: setWindowParam: Warning - beta fixed to 0 "
                 "(should be between 0 to 1)");
            beta = 0;
        }
        
        x->sp_windowAlpha = alpha;
        x->sp_windowBeta = beta;
        wfs_setWindow(x);
    }
    else {
        post("wfs~: setWindowParam: Warning - Window size must be a positive integer");
    }
}

/************************ print methods *************************************/
static void wfs_printAll (t_wfs_tilde *x)
{
    wfs_printVersion(x);
    wfs_printCopyright(x);
    wfs_printReference(x);
    wfs_printSpeaker(x);
    wfs_printSpeakerVec(x);
    wfs_printRadVec(x);
    wfs_printSource(x);
    wfs_printDistance(x);
    wfs_printGain(x);
    wfs_printDelay(x);
    post("");
}

static void wfs_printCopyright (t_wfs_tilde *x)
{
    post("wfs~: %s", COPYRIGHT);
}
 
static void wfs_printGain (t_wfs_tilde *x) 
{
    int i = 0;
    for (i = 0; i < x->sp_n; i++)
        post("wfs~: Gain %d: %f", i+1, x->sp[i]->sp_gain);
}
static void wfs_printDelay (t_wfs_tilde *x) 
{
    int i = 0;
    for (i = 0; i < x->sp_n; i++)
        post("wfs~: Delay %d: %f", i+1, x->sp[i]->sp_delay);
}

static void wfs_printDistance (t_wfs_tilde *x) 
{
    int i = 0;
    for (i = 0; i < x->sp_n; i++)
        post("wfs~: Distance %d: %f", i+1, x->sp[i]->sp_dist);
}

static void wfs_printHelp (t_wfs_tilde *x)
{
    wfs_printVersion(x);
    wfs_printCopyright(x);
    post(USAGE);
}

static void wfs_printInit (t_wfs_tilde *x) 
{
    wfs_printVersion(x);
    /* wfs_printCopyright(x); */
    wfs_printReference(x);
    wfs_printSpeaker(x);
    wfs_printSource(x);
    post("");
}

static void wfs_printRadVec (t_wfs_tilde *x)
{
    int i;
    for (i = 0; i < x->sp_n; i++) {
        post("wfs~: Radiation Vector %d: %f", i+1, x->sp[i]->sp_radVec);
    }
}

static void wfs_printReference (t_wfs_tilde *x) 
{
    post("wfs~: Reference Z-axis: %2.2f", x->x_zref);
}

static void wfs_printSpeaker (t_wfs_tilde *x) 
{
    int i = 0;
    post("wfs~: Num. of Speakers: %d", x->sp_n);
    post("wfs~: Speaker Spacing: %0.4f m", x->sp_dx);

    for (i = 0; i < x->sp_n; i++) {
        if (x->x_sdim == TWOD)
            post("wfs~: Speaker %d: (%2.2f, %2.2f)", i+1,
                    x->sp[i]->sp_loc[X], x->sp[i]->sp_loc[Z]);
        else
            post("wfs~: Speaker %d: (%2.2f, %2.2f, %2.2f)", i+1,
                x->sp[i]->sp_loc[X], x->sp[i]->sp_loc[Y], x->sp[i]->sp_loc[Z]);
    }

    post("wfs~: Window: %s", x->x_swindow ? "ON":"OFF");
    post("wfs~: Worst Aliasing Frequency: %5.2f Hz", x->x_afreq);   
}

static void wfs_printSpeakerVec (t_wfs_tilde *x)
{
    int i = 0;
    for (i = 0; i < x->sp_n; i++) {
        if (x->x_sdim == TWOD)
            post("wfs~: Vector %d: (%2.2f, %2.2f)", i+1,
                    x->sp[i]->sp_vec[X], x->sp[i]->sp_vec[Z]);
        else
            post("wfs~: Vector %d: (%2.2f, %2.2f, %2.2f)", i+1,
                x->sp[i]->sp_vec[X], x->sp[i]->sp_vec[Y], x->sp[i]->sp_vec[Z]);
    }
}

static void wfs_printSource (t_wfs_tilde *x) 
{
    int i;

    post("wfs~: Source Type: %s", x->x_stype ? "Plane Wave":"Point");

    if (x->x_stype == PLANEWAVE)
        post("wfs~: Source Angle: %0.2f", x->x_src->sr_angle * 360);

    if (x->x_stype == POINTSOURCE) {
        post("wfs~: Source Patter: %s", x->x_scard ? "Cardioid":"Omni");

        if (x->x_scard)
            post("wfs~: Source Angle: %0.2f", x->x_src->sr_angle * 360);

        for (i = 0; i < x->sp_n; i++) {
            post("wfs~: Source Angle to Speaker %d: %0.2f", 
                i, x->sp[i]->sp_phi * 360);
        }
    }

    if (x->x_sdim == TWOD)
        post("wfs~: Source Location: (%2.2f, %2.2f)",
            x->x_src->sr_loc[X], x->x_src->sr_loc[Z]);
    else    
        post("wfs~: Source Location: (%2.2f, %2.2f, %2.2f)",
            x->x_src->sr_loc[X], x->x_src->sr_loc[Y], x->x_src->sr_loc[Z]);
}

static void wfs_printVersion (t_wfs_tilde *x) 
{
    post("wfs~: Version:  %s", VERSION);
}

static void wfs_printWindow (t_wfs_tilde *x)
{
    int i;
    
    if (x->x_swindow) {
        post("wfs~: Window weighting %d speakers on each edge", x->sp_windowCount);
        post("wfs~: Alpha: %0.3f, Beta: %0.3f", x->sp_windowAlpha, x->sp_windowBeta);
        for (i = 0; i < x->sp_n; i++)
            post("wfs~: Amplitude Weighting for %d: %0.3f", i, x->sp_window[i]);
    }
    else {
        post("wfs~: No Window is applied");
    }
}

