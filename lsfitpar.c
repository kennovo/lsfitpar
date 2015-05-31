/* $Id: lsfitpar.c,v 1.85 2015/05/16 23:45:40 kenno Exp $ */
/**********************************************************************\
*  Copyright (C) 2011-2015 Kenno Vanommeslaeghe                        *
*                                                                      *
* This program is free software: you can redistribute it and/or modify *
* it under the terms of the GNU Affero General Public License as       *
* published by the Free Software Foundation, either version 3 of the   *
* License, or (at your option) any later version. Please review the    *
* file COPYING that you should have received along with this program;  *
* it contain a preamble with important information.                    *
*                                                                      *
* This program is distributed in the hope that it will be useful,      *
* but WITHOUT ANY WARRANTY; without even the implied warranty of       *
* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the         *
* GNU Affero General Public License for more details.                  *
\**********************************************************************/

/* ACHTUNG: Do not forget to change the version string!!! */
const char fullname[] = "ParamChem least-squares fitting";
const char version[] = "0.9.1 beta";
const char reldate[] = "the 31st of May 2015";

#include <stdio.h>   /* duh! */
#include <stdint.h>  /* int32_t , int64_t (not interested in pre-C99 compatibility) */
#include <stdlib.h>  /* for malloc(), realloc(), free(), qsort(), bsearch() */
#include <string.h>  /* for strlen(), strncmp(), strncpy(), memcpy() */
#include <math.h>    /* for sqrt() and cosf()/cos() */

#include "gopt.h"  /* Public domain library for parsing command-line arguments. */

double alpha=1.0;
double beta=0.0;
/* int64_t is the default because we're almost exclusively running on little-endian
 * machines, so that a pointer to a memory location containing a 64-bit integer can be
 * treated as a pointer to a 32-bit integer without adverse effects. The opposite is not
 * necessarily true, and the scheme breaks down for big-endian (which is becoming
 * increasingly rare). */
#ifdef FINT32
#ifdef FINT64
#error Preprocessor macros FINT32 and FINT64 are mutually exclusive
#endif  /* FINT64 */
typedef int32_t fint;
#else   /* FINT32 */
typedef int64_t fint;  /* default */
#endif  /* FINT32 */
fint one=1;
char trans='T';
char notrans='N';
char lo='L';
/* dgemv_ from BLAS; dgels_ and dsfrk_ from LAPACK. Any compliant implementation will
 * do; we just have to point to it in the Makefile. */
void dgemv_(char *trans, fint *m, fint *n, double *alpha, double *a, fint *lda,
		double *x, fint *incx, double *beta, double *y, fint *incy);
/* we use DGELS because our null bias guarantees that the system is not underdetermined
 * (and the expert least-square routines are difficult to understand) */
void dgels_(char *trans, fint *m, fint *n, fint *nrhs, double *a, fint *lda,
		double *b, fint *ldb, double *work, fint *lwork, fint *info);
/* this is for generating the <f|g> matrix for determining the null bias */
void dsfrk_(char *transr, char *uplo, char *trans, fint *n, fint *k, double *alpha,
		double *a, fint *lda, double *beta, double *c);

/* We don't need a fully C99-compliant implementation for booleans and
 * inlining to work in the context of this program */
#ifdef __GNUC__
#if (__GNUC__ == 2 && __GNUC_MINOR__ >= 95) || __GNUC__ > 2
#define RUDIMENTARYC99 1
#endif    /* ... __GNUC__ > 2 */
#elif __STDC_VERSION__ >= 199901L    /* __GNUC__ */
#define RUDIMENTARYC99 1
#endif    /* ifdef __GNUC__  elif __STDC_VERSION__ >= 199901L */

#ifndef _STDBOOL_H  /* there are booleans in global data structures! */
#ifdef RUDIMENTARYC99
#include <stdbool.h>
#elif defined __cplusplus    /* RUDIMENTARYC99 */
#define _Bool   bool
#else    /* __cplusplus */
#define false   0
#define true    1
#define bool    _Bool
typedef int     _Bool;
#endif    /* ifdef RUDIMENTARYC99 elif __cplusplus */
#define _STDBOOL_H
#define __bool_true_false_are_defined   1
#endif    /* stdbool.h */

#define DEFBIAS 0.001
/* In CGenFF, the average dihedral K is ~1.5 . We could argue that for the current
 * calculation, we shoud use 3 because people will be instructed not to fit the aliphatic
 * hydrogens that pull the average down, but so far, we've obtained better results with
 * higher factors. The average angle and improper Ks are both ~60, so the factor for
 * dividing the biasing restraints is 60 / 1.5 = */
/* TODO: these factors are rather empirical. While adjusting them is not trivial and can
 * be dangerous, the user should be given at least some way to change them. */
#define ANGIMPF 40
/* The average for the bonds is ~300 */
#define BONDF 200

/* The following quantities are allocated dynamically;
 * we just want to keep them within sensible limits */
#define MAXRGRP     4096 /* maximum number of RMSD groups */
#define MAXDIH      992  /* maximum number of dihedrals. This automatically fixes the
			    maximum number of parameters to 992*6 = 5952 */

/* The following quantities are allocated statically */
#define MAXLINE    1024  /* full path to any file may not be longer than this (UDF) */
#define MAXGREQ    512  /* maximum number of (identical) strictly positive integers
			 that can be put in a 1024-byte space-separated list */
/* maximum size of target/trial vector including restraints; 4096 is 32k but Xiao told me
 * he sometimes fits his whole data set at once, comprising of 9000~10000 data points.
 * Note that higher values than 16384 require changing at least some data types
 * because (nt=ns+np) */
#define OPT_MAXT 16384
#define DIHSTR " measurement 1234"
#define DIHOFF 13  /* offset */
#define DIHF "%-4i"
#define HASHSZ 8
#define ATTF "8s"        /* Format of CHARMM atom type */
#define ATTL HASHSZ + 1       /* 1 for \0 */
#define DHASH HASHSZ * 4 + 4  /* 3 ' ' + 1 '\0' */
#define LASTATT HASHSZ * 3 + 2  /* index where the '\0' is (option base 0) */
#define ANGATT HASHSZ * 2 + 1  /* index where the '\0' is (option base 0) */
#define FNCHUNK 31744  /* Fits 31 maximum-size dihedrals (31*1024).
			* 2 of these chunks easily fit in 64k. */
#define MAXFNC 32  /* 31*32 = 992 dihedrals */
#define NMULT 7  /* including mult[0] , which is reserved for later use */
/* The following comes from experimental/deg.c and has more than sufficient precision for
 * a long double (ie. 0.0174532925199432957691L). However, we deliberately define it as
 * double (no suffix; see K&R A2.5.3) and leave the truncation to the compiler. */
#define DEG_F 0.0174532925199432957691
#define F_DEG   57.2957795130823208767
#define MINRANGE 0.001  /* Smallest allowable range for any measurement to be numberically stable */
/* We want to exlude all-around scans with step size 30 deg (the minimum acceptable step
 * size for a 6-fold with fixed phase or a 3-fold with variable phase) but not 45 deg
 * (which would be a genuine gap no matter who you ask). 45 deg = pi/4 rad ==> highest
 * possible fitting coefficients are 4 and 3 if the minimum is within the range. */
#define MAXRANGE 315.0
#define SMALLSHIFT 3.0
const char shiftwarn[] = " WARNING: phase shift > 3deg!";
const char oobwarn[] = " fitted reference value out of acceptable range";
/* The range at which the resultant phase of a 2-harmonic correction term on top of an
 * initial guess harmonic can be considered accurate, in multiples of the scan range.
 * Since physically spoken, these resultants can get arbitrary large for an arbitrary
 * close initial guess amplitude, this number is only limited by the numerical precision
 * of the LLS procedure; setting it lower would just produce spurious warnings. */
/* TODO: although we have been able to artificially trigger large imprecision in the
 * past, those results cannot be used to pin down this number. We've had a case with a
 * multiplier of 200 that did *not* cause any imprecision with respect to the output
 * rounding, but the number can and should probably be set far larger; 200 is bound to
 * yield the aforementioned spurious warnings. */
#define PRESRANGE 200.0
/* The allowed deviation of a fitted harmonic without initial guess with respect to the
 * center of the scan range will be the distance from the center of the scan range to its
 * edge multiplied by a factor 1 + 2 * EXTRANGE. For example, if the user did a 3-point
 * scan with step size 3deg on an angle, a fitted phase that deviates 3 * (1 + 3) = 12
 * degrees from the center of the scan range will be accepted. Or in other words, for a
 * minimized QM value of 109, reference values in the 97-121 range will be accepted,
 * which sounds about right. */
/* TODO: EXTRANGE does not get checked if an initial guess is present; in that case,
 * only PRESRANGE applies. This is a result of the organic growth of this
 * program and obviously doesn't make much sense (and neither does the way different
 * range checks are applied in different stages of the process in light of initial guess
 * harmonics). The range checking features throughout the program should probably be
 * overhauled to reflect the new reality. */
#define EXTRANGE 1.5

#define PHASEFLAG 128
#define BANG 129  /* bond or angle ==> bit 1 for 0-fold and 128 for phase */
#define COSINES 126  /* all the cosine functions */

enum datatyp { INT, FLOAT, DOUBL, STRING };  /* INT must come first so that it is 0 */
enum ouinon { YNVOID=-1, NO, YES };  /* NO must be 0 so that we can treat is as bool */
enum biasmode { UNIFORM=1, ADAPT, ABS };  /* Corresponds to what is asked from the user */
/* enum phaseflag { NOMULT, PH_UNDEF, PH_FIX, PH_VAR }; /* obsolete version with NOMULT=0 */
enum phaseflag { PH_UNDEF, PH_FIX, PH_VAR };

struct dihedral {
    char *name;
    struct dihedral *next;  /* next in equivalence group */
    int tmult;  /* -1: no target */
    double targ1,targ2;  /* float in "cosf" branch */
    char hash[DHASH];
};

struct eqgroup {
    struct dihedral *first;
    double delta1,delta2,dmin,dmax;  /* float in "cosf" branch */
    int size;
    unsigned char mult;
    /* bit-fields are in principle bad for performance but these are flags that
     * are only set and queried a few times and only during I/O. */
    unsigned int merged : 1;
    unsigned int ph_undef : 1;
};

struct normal {
    double avg;
    int cnt;
} *norm;  /* global because allocated in main() and used in normalize_weight() */

int *groupvec=NULL,*maxgrpv;  /* global for the same reason */
double *weighvec=NULL;        /* global for the same reason */
int ghi;                      /* global for the same reason */
struct eqgroup *eqgrp,*eqgmax;  /* allocated in main() and used in findgrp() listeqg() */
struct dihedral *dih;      /* global because allocated in main() and used in listeqg() */
fint np;        /* global because used in main() and packmult() */
void *options;  /* global because used in main() and outsideopt() */
bool interact=false;
const char helpstr[] = "Performs least-squares fit of bonded molecular mechanics parameters to\n"
"conformational energies and corresponding measurements, applying\n"
"robustness-enhancing restraints described in [1] K. Vanommeslaeghe,\n"
"M. Yang, A. D. MacKerell Jr., J. Comput. Chem. 2015, DOI:10.1002/jcc.23897 .\n"
"Please cite the above paper when publishing results obtained with this program!\n"
"\n"
"Usage: Print help or version information: %s [-h|--version]\n"
"       Run fully interactively: %s [-i] <prm file> [<ene file>]\n"
"       Run from command-line: %s [-i] <qm energies> <mm energies>\n"
"  <prm file> [<ene file>] [-w <weights>] [-g <groups>]\n"
"  [-u|-t|--abs] [-b <bias>] [--no-comp] [-e <measurements>]...\n"
"  [-p '<parameter>' [-q '<parameter>' | -m '<multiplicities>' [-f|-a]]]...\n"
"\n"
"Where: -h, -?, --help, --HELP displays this help message and exits\n"
"       --version displays version and licensing information, then exits\n"
"       <qm energies>: ASCII file containing 1 QM energy / line (m(*) lines)\n"
"       <mm energies>: ASCII file containing 1 MM energy / line (m(*) lines)\n"
"       -i, --interactive requests interactive mode, which also is automatically\n"
"  assumed when 1 or 2 command-line arguments are given.\n"
"       -w, --conf-weights <weights>: specify ASCII file containing\n"
"  1 weight factor / line (m(*) lines)\n"
"       -g, --groups <groups>: specify ASCII file containing\n"
"  1 integer group ID / line (m(*) lines)\n"
"       -u --uniform requests uniform bias (see [1] above). This is the default.\n"
"       -t --target-adapted requests target-adapted (see [1] above).\n"
"       --abs requests \"absolute value\" bias. This is not the constant bias\n"
"  in [1], but rather an experimental variation of the target-adapted bias that\n"
"  did not perform well enough to make in into the paper. Not recommended - only\n"
"  for compatibility with runs performed during alpha testing.\n"
"       -b, --glob-bias <bias>: specify global bias fraction. Default: %g .\n"
"       --no-comp disables bias compensation. Not recommended.\n"
"       -e --measure <measurements>: specify ASCII file containing all m(*)\n"
"  measurements of a specific degree of freedom (DF; bond length, angle,...),\n"
"  1/line. The first line contains the atom types that define the parameter\n"
"  associated with this DF, optionally along with its initial guess, for a total\n"
"  of m+1 lines. This flag is repeated as many times as there are relevant DF,\n"
"  carefully including *all* DF associated with parameters being fitted, even DF\n"
"  that are not actively scanned!\n"
"       -p --parameter '<parameter>': specify a parameter being fitted by its\n"
"  defining atom types (enclosed in quotes). Starts a parameter section\n"
"  consisting of either a -q flag or a parameter definition containing 0 or more\n"
"  of the other flags detailed below.\n"
"       -q --equivalent-to '<parameter>': make the parameter specified by\n"
"  the -p preceding -q equivalent to the parameter specified by -q .\n"
"       -m --multiplicities '<multiplicities>': fit one or more multiplicities\n"
"  (e.g. 246) ranging from 0 to 6, where 0 is a harmonic potential that should\n"
"   be used exclusively for bonds, angles and improper dihedrals, and is the\n"
"  default for these types of parameters. Currently, a parameter section can\n"
"  contain no more than one -m flag.\n"
" /     -f --fix-phase : fix phases for selected multiplicities.\n"
"/      -a --vary-phase : allow phases for selected multiplicities to vary.\n"
"\\ In the absence of -f or -a flags, bonds are assumed variable and\n"
" \\dihedrals and impropers fixed.\n"
"\n"
"(*) m is the number of conformations, as discussed in [1] above.\n"
"\n"
"Note: specifying options in a different order than above is not supported.\n"
"Tips: - Long command lines can be replaced with orderly input files by using\n"
"the 'xargs' UNIX tool as follows: xargs ./lsfitpar <lsfit_input_file\n"
"      - If you're on a terminal that doesn't scroll and you're seeing only this\n"
"part, consider piping the output to \"less\" or redirecting it to a file. :-P\n";
char empty='\0';  /* The pointer to this doubles as a sentinel value. */

/* TODO LATER: eliminate more of the duplicated code */
bool readvect (const char *filename, enum datatyp kind,
		void *vec, fint *dim, struct dihedral *dihe, char *disp);
struct eqgroup *findgrp(const char *buf, int swit);
bool parhash(char *hash, char att[][ATTL], int nat);
int get_list (char *line, enum datatyp kind, void *array, int minsz, int maxsz,
		int minint, int maxint, char *disp_s, char *disp_p);
void normalize_weight (double *tvec);  /* a tad too heavy for static inlining */
void finalize_mult (struct eqgroup *eqg, int *mult, int nmult, enum phaseflag phase);
char packmult (int *mult, int nmult, bool phasevar);
/* no inlining for functions that are involved in any kind of I/O (disk is not that fast) */
const char *outsideopt (const char exact, const char fuzzy, const char *disp);
void listeqg (FILE *unit, char *pre, struct eqgroup *eqg, int indent);
bool promptenter (char *buf, const char *filename, enum datatyp kind,
		void **vect, fint *dim, char *disp, char *action);
bool promptread (char *buf, const char *filename, enum datatyp kind,
		void *vec, fint *dim, struct dihedral *dihe, char *disp);
FILE *fopen_clean(const char *filename, char *mode,
			char *function, char *disposition);
size_t datatypsize (enum datatyp kind);
enum ouinon yesno (char *line, char *disp, enum ouinon def);
bool strnupper (char *string, int n);
/* can't inline this because it will be passed to a library function */
int intcmp (const void *, const void *);  /* Compare 2 integers */
int floatcmp (const void *, const void *);  /* Compare 2 double (float in cosf branch) */
int numl(const char *s);  /* no inline because only called in I/O */

int main (int argc, const char *argv[]) {
    char *fn[MAXFNC];  /* Array with pointers to filename chuncks */
    char *curfnp;      /* pointer to current position within filename chunck */
    unsigned int fnfree;  /* first free position in current filename chunck */
    int fns;  /* Actual number of filename chuncks; maximum = MAXFNC */
    char **cppcur,**cppmax;

    const char *qmef,*mmef,*ccp,*comm;
    char *cpa,*cpb,*cpc;
    int *curgrpv,*eqextra,*ipa,*ipb;  /* do not rename */
    /* TODO: some of these may have become redundant. However, we don't want to replace
     * curdih with a general-purpose one because that would encumber merging branches. */
    double *targvec,*qme,*mme,*dihmat,*phasetmp,*curdih,*cosmat,*parvec,*amat,*work,*nbvec,*nbtmp,*solvec,*biasvec,*dpa,*dpb,*dpc,*dpd,*dpe,*dpf,*dpg,*dph,*dpi;  /* TODO: if applicable, rename dpa to curtarg, dpb to curweigh and dpc to targmax when finished */
    FILE *outfile=NULL;
    const opt_t *opts;
    struct dihedral *dipa,*dipb,*dipc,*maxdih;  /* do not rename */
    struct eqgroup *eqga,*eqgb,*eqgc,*eqgd;  /* do not rename */
    struct normal *qmnorm,*curqmn,*curnor,*maxnor;

    /* ns,np,nt = number of scan points, parameters, target data points */
    fint ns,npp,np2,nt,lwork,hdim;
#ifdef DEBUG
    fint hdimm;
#endif
    fint info=0;  /* need to initialize this to avoid uninit bug */
    size_t scansz;
    /* d1,d2,delta will be float in "cosf" branch */
    double biasfract,d,d1,d2,da,db,delta,e,x,y,sum,t1,t2;
    int multemp[NMULT];
    int ndih=0,negrp,maxeq0size;
    int exitstatus,i,j,k,l,m,n,o,p;
    enum ouinon phasedef,biascomp;
    enum biasmode bmod;
    enum phaseflag pf;
    bool f,g,h,b;
    char line[MAXLINE];
    unsigned char u,v;

    char *dihnum;
    char dihdisp[] = DIHSTR;

    mmef=qmef=NULL;
    dihnum = dihdisp + DIHOFF;

    /* correct merging of stdout and stderr is essential */
    setvbuf(stderr, NULL, _IOLBF, BUFSIZ);

    options = gopt_sort(&argc, argv, gopt_start(
	gopt_option('h', 0, gopt_shorts('h','?'), gopt_longs("help","HELP")),
	gopt_option( 1 , 0, gopt_shorts(0), gopt_longs("version")),
    /* Verbosity not implemented yet */
    /*  gopt_option('q', 0, gopt_shorts('q'), gopt_longs("quiet")),
	gopt_option('v', GOPT_REPEAT, gopt_shorts('v'), gopt_longs("verbose")), */
	gopt_option('i', 0, gopt_shorts('i'), gopt_longs("interactive")),
    /* long version of ambiguous form deliberately undocumented */
	gopt_option('w', GOPT_REPEAT | GOPT_ARG, gopt_shorts('w'), gopt_longs("weight")),
	gopt_option( 2 , GOPT_ARG, gopt_shorts(0), gopt_longs("conf-weights")),
    /* Not implemented yet: */
    /* gopt_option( 3 , GOPT_REPEAT | GOPT_ARG, gopt_shorts(0), gopt_longs("mult-weight")), */
    /* long version of ambiguous form deliberately undocumented */
	gopt_option('b', GOPT_REPEAT | GOPT_ARG, gopt_shorts('b'), gopt_longs("bias")),
	gopt_option( 4 , GOPT_ARG, gopt_shorts(0), gopt_longs("glob-bias")),
    /* Not implemented yet: */
    /* gopt_option( 5 , GOPT_REPEAT | GOPT_ARG, gopt_shorts(0), gopt_longs("mult-bias")), */
	gopt_option('g', GOPT_ARG, gopt_shorts('g'), gopt_longs("groups")),
	gopt_option('u', 0, gopt_shorts('u'), gopt_longs("uniform")),
	gopt_option('t', 0, gopt_shorts('t'), gopt_longs("target-adapted")),
	gopt_option( 6 , 0, gopt_shorts(0), gopt_longs("abs")),
	gopt_option( 7 , 0, gopt_shorts(0), gopt_longs("no-comp")),
	gopt_option('e', GOPT_REPEAT | GOPT_ARG, gopt_shorts('e'), gopt_longs("measure")),
	gopt_option('p', GOPT_REPEAT | GOPT_ARG, gopt_shorts('p'), gopt_longs("parameter")),
	gopt_option('q', GOPT_REPEAT | GOPT_ARG, gopt_shorts('q'), gopt_longs("equivalent-to")),
	gopt_option('m', GOPT_REPEAT | GOPT_ARG, gopt_shorts('m'), gopt_longs("multiplicities")),
	gopt_option('f', GOPT_REPEAT, gopt_shorts('f'), gopt_longs("fix-phase")),
	gopt_option('a', GOPT_REPEAT, gopt_shorts('a'), gopt_longs("vary-phase"))));

    if (gopt(options, 'h')) {
	/* The following would be for the final helpstr as in the file builtinhelp */
	/* if ((i = strlen(ccp=argv[0]) - 3) < 0) i = 0;
	/* printf (helpstr,ccp,ccp,ccp,ccp,  /* next line: 8 */
	/* 	i,&empty,i,&empty,i,&empty,i,&empty,i,&empty,i,&empty,i,&empty,i,&empty,
	/* 	i,&empty,i,&empty,i,&empty,i,&empty,i,&empty,i,&empty,
	/* 	DEFBIAS,ANGIMPF,BONDF); */
	ccp=argv[0];
	printf (helpstr,ccp,ccp,ccp,DEFBIAS);
	return 0;
    }
    if (! (i=gopt(options,1))) {  /* abuse i to mean "--version" */
	interact = gopt(options, 'i');
	switch (argc) {
	  case 5:
	    if (! (outfile = fopen_clean(argv[4],"w","","delta"))) return 1;
	    /* ACHTUNG: fall through to case 4 */
	  case 4:
	    qmef=argv[1];
	    mmef=argv[2];
	    argv[1]=argv[3];
	    break;  /* out of switch */
	  case 3:
	    if (! (outfile = fopen_clean(argv[2],"w","","delta"))) return 1;
	    /* ACHTUNG: fall through to case 2 */
	  case 2:
	    interact = true;
	    break;  /* out of switch */
	  case 1: case 0:
	    fprintf(stderr,"ERROR: insufficient arguments; run with -h "
	                    "for help (%i lines).\n",numl(helpstr));
	    return 1;
	  default:
	    fprintf(stderr,"ERROR: too many arguments; run with -h "
	                    "for help (%i lines).\n",numl(helpstr));
	    return 1;
    }   }
    if (i || interact) {  /* TODO once we have verbosity: || verbosity >= 0 */
	printf ("%s version %s\nreleased %s\n"
"Please cite the above version information as well as the following paper\n"
"when publishing results obtained with this program: K. Vanommeslaeghe,\n"
"M. Yang, A. D. MacKerell Jr., J. Comput. Chem. 2015, DOI:10.1002/jcc.23897 .\n\n"
" Copyright (C) 2011-2015 Kenno Vanommeslaeghe\n\n"
"This program is free software: you can redistribute it and/or modify\n"
"it under the terms of the GNU Affero General Public License as\n"
"published by the Free Software Foundation, either version 3 of the\n"
"License, or (at your option) any later version. Please review the\n"
"file COPYING that you should have received along with this program;\n"
"it contain a preamble with important information.\n\n"
"This program is distributed in the hope that it will be useful,\n"
"but WITHOUT ANY WARRANTY; without even the implied warranty of\n"
"MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the\n"
"GNU Affero General Public License for more details.\n\n",
	    fullname,version,reldate);
	if (i) return 0;
	  else fflush(stdout);
    }
    if (! ((ndih=gopt(options,'e')) || interact)) {  /* order matters! */
	fprintf(stderr,"ERROR: no measurements specified in non-interactive mode; "
			"run with -h for help (%i lines).\n",numl(helpstr));
	return 1;
    }

    /***********************\
    * Read QM and MM energy *
    \***********************/
    if (! (qme = malloc(sizeof(double) * OPT_MAXT))) {
	fprintf(stderr,"ERROR: not enough memory to initialize program!\n");
	return 1;
    }
    ns=0;
    if (promptread(line,qmef,DOUBL,qme,&ns,NULL,"QM energy ")) return 1;
    /* Shrink qme to fit. It's of no use to check if ns > OPT_MAXT because
     * we'd be dead already (ie. readvect would have segfaulted) */
    /* To be safe, we have to store the return value of realloc back into qme
     * for in case realloc decides to move the block elsewhere. While this may
     * seem unlikely, valgrind 3.10's realloc() actually does this! */
    if (ns != OPT_MAXT) qme=realloc(qme,sizeof(double) * ns);
    if (! (mme = malloc(sizeof(double) * ns))) {
	fprintf(stderr,"ERROR: can't allocate memory for MM energies!\n");
	return 1;
    }
    if (promptread(line,mmef,DOUBL,dpb=mme,&ns,NULL,"MM energy ")) return 1;

    if (outfile) {  /* allocate a separate targvec so that we can preserve qme */
	if (! (dpc = targvec = malloc(sizeof(double) * ns))) {
	    fprintf(stderr,"ERROR: can't allocate memory for target vector!\n");
	    return 1;
	}
	for (dpd = (dpa=qme) + ns; dpa<dpd; *(dpc++) = *(dpa++) - *(dpb++));
    } else {  /* overwrite qme with QM-MM so that we can free mme */
	for (dpd = (dpa=targvec=qme) + ns; dpa<dpd; *(dpa++) -= *(dpb++));
	free(mme);
#ifdef DEBUG
	mme=NULL;  /* Prevent it from being freed again at the end to shut up valgrind */
#endif
    }
    if ((ccp=outsideopt(2,'w',"file with per-conformation weight factors")) == &empty)
	return 1;  /* next line: outsideopt() returned NULL if no weights on cmd line */
    if (promptenter(line,ccp,DOUBL,(void **) &weighvec,&ns,"weight ","weights")) return 1;
    /* Not safe: for (dpb = (dpa=weighvec) + ns; dpa < dpb; *dpa++ = sqrt(*dpa)); */
    /* A failed sqrt() sets errno but a successful one doesn't reset it. Therefore,
     * we could check for it outside the loop, like this: if (errno == EDOM) { ... }
     * ... except that -ffast-math suppresses sqrt()'s ability to set errno. And
     * checking for NaNs still works but needs to be done inside the loop which
     * negates most of the performance advantage, plus it's even less portable
     * because NaN is not defined in all C implementations; see
     * ../cgenff/experimental/README.root_neg . The following is really the only
     * way out; it's slower but we don't care much as we're still in I/O. */
    if (weighvec) for (dpb = (dpa=weighvec) + ns; dpa < dpb; ) {
	    if ((d=*dpa) < 0) {
		fprintf(stderr,"ERROR: negative weight factors not allowed!\n");
		return 1;
	    }
	    *dpa++ = sqrt(d);
	}

    /*************************************************************************\
    * Read group file and apply normalization and (optionally) weight factors *
    \*************************************************************************/
    ccp=NULL;
    gopt_arg(options, 'g', &ccp);  /* next line: groupvec initialized to NULL above */
    if (promptenter(line,ccp,INT,(void **) &groupvec,&ns,"RMSE group ","group fitting")) return 1;
    if (groupvec) {
	/* The right thing to do would be to translate the external group symbols
	 * into internal group integers but it's unlikely any sensible user will
	 * break this simpler implementation */
	j=MAXRGRP;  /* find lowest and highest group id */
	for (maxgrpv = (ipa=groupvec) + ns; ipa < maxgrpv; ) {
	    i=*(ipa++);
	    if (i < j) j = i;
	    if (i > ghi) ghi = i;
	}
	j = -j;
	if ((ghi+=j+1) < 1) {
	    fprintf(stderr,"ERROR: no groups read!\n");
	    return 1;
	}
	if (ghi > MAXRGRP) {
	    fprintf(stderr,"ERROR: read %i groups (maximum %i)!\n",ghi,MAXRGRP);
	    return 1;
	}  /* next line: translate group ids so that lowest is 0 */
	for (ipa=groupvec; ipa < maxgrpv; *(ipa++)+=j);
    } else {  /* groupvec must be allocated even if we don't use group fitting */
	if (! (groupvec = malloc(sizeof(int) * ns))) {
	    fprintf(stderr,"ERROR: can't allocate memory for RMSE alignment!\n");
	    return 1;
	}
	for (maxgrpv = (ipa=groupvec) + ns; ipa < maxgrpv; *(ipa++)=0);
	ghi=1;  /* # of groups */
    }
    if (! (norm = malloc(sizeof(struct normal) * ghi))) {
	fprintf(stderr,"ERROR: can't allocate memory for group normalization!\n");
	return 1;
    }
    /* This is a different ball game than aligning QM and MM - we'll do that later */
    normalize_weight (targvec);
    sum=0.0;
#ifdef DEBUG
    /* The debug code in this loop dumps the whole target vector (vertically) after
     * running normalize_weight (which applies weight, groups,...) This is potentially a
     * lot of lines (think of 1000-point scans), we're currently pretty confident
     * normalize_weight is sane, and we're not planning any changes to it. For these
     * reasons, it's still commented out even when compiling with -D DEBUG . We might
     * want to enable it if we ever implement a -vvv flag. */
    /* printf("%i:\n",ns);  /* DEEP DEBUG */
#endif
    for (dpb = (dpa=targvec) + ns; dpa < dpb; ) {
	d = *(dpa++);
#ifdef DEBUG
	/* printf("%g\n",d);  /* DEEP DEBUG */
#endif
	sum += d * d;
    }
    printf("Initial RMSE = %g\n",sqrt(sum/ns));
    /****************************************************\
    * Read measurements and determine equivalence groups *
    \****************************************************/
    /* If we're non-interactive, we'd already have issued an error if ndih==0,
     * but if we are interactive, it's still possible that -e options were
     * present so that we don't need to query... */
    while (! ndih) {  /* allow the user to retry when giving invalid input */
        printf("Number of measurements: ");
	fflush(stdout);
	/* We first used a simple scanf statement but then discovered this was a bad
	 * idea; see http://c-faq.com/stdio/scanfinterlace.html and the comment in
	 * get_list. */
	/* Yes, get_list will work on a list of length 1. Most of the code for
	 * handling longer lists will be skipped by conditional clauses, so from a
	 * performance point of view, it's not even terrible. */
	/* The last string should actually read "numbers of measurements", but we
	 * determined that most of the warning messages that might realistically be
	 * encountered will look better if we keep it singular. */
	switch (get_list(line,INT,&ndih,1,1,1,MAXDIH,
	        "number of measurements","number of measurements")) {
	  case -3: return 1;  /* EOF */
	  case -2: case -1: case 0: ndih=0;  /* force retry */
    }   }
    /* We will only know the number of *parameters* after all measurements have been
     * read, equivalences have been established, and multiplicities have been assigned,
     * so we are forced to temporarily store the measurement data. */
    if (! (curdih = dihmat = malloc(sizeof(double) * ns * ndih))) {
	fprintf(stderr,"ERROR: can't allocate memory for measurements!\n");
	return 1;
    }
    if (! (dipa = dih = malloc(sizeof(struct dihedral) * ndih))) {
	fprintf(stderr,"ERROR: can't allocate memory for measurement metadata!\n");
	return 1;
    }
    fns=0;
    fnfree=FNCHUNK;  /* Trigger allocation on next measurement */
    ccp=NULL;
    if (! interact) opts=options;
    for (i=1; i<=ndih; i++) {
	if (! interact) for (;;opts++) {
		switch (opts->key) {
		  case '\0':  /* ndih was set to gopt(options,'e') ... */
		    fprintf(stderr,"ERROR: consistency error #3; "
				    "please contact a developer!\n");
		    return 1;
		  case 'e':
		    ccp=opts++->arg;  /* ++ for next iteration of outer (ndih) for */
		    break;  /* out of switch (and therefore out of for) */
		  default: continue;  /* only way to repeat the inner loop */
		}
		break;
	    }
	sprintf(dihnum,DIHF,i);
	if (promptread(line,ccp,DOUBL,curdih,&ns,dipa,dihdisp)) return 1;
	curdih+=ns;  /* next line: abuse comm for a bit */
	j = strlen(comm = ccp ? ccp : line) + 1;  /* 1 for \0 */
	if ((fnfree+=j) > FNCHUNK) {
	    if (fns == MAXFNC) {
		fprintf(stderr,"ERROR: maximum memory allocation for "
				"file names exceeded!\n");
		return 1;
	    }
	    if (! (fn[fns++] = curfnp = malloc(sizeof(char) * FNCHUNK))) {
	        fprintf(stderr,"ERROR: can't allocate memory for "
				"measurement file names!\n");
	        return 1;
	    }
	    fnfree=j;
	}
	strncpy(dipa->name=curfnp,comm,j);
	dipa++->next=dih;  /* we use dih as a sentinel value */
	curfnp+=j;
    }
    if (! (eqgb = eqgrp = malloc(sizeof(struct eqgroup) * ndih))) {
	fprintf(stderr,"ERROR: can't allocate memory for equivalence groups!\n");
	return 1;
    }
    /* Recognition of equivalence groups */
    /* Respecting the input order is way more intuitive than sorting by hash and we
     * have to manually compare hashes no matter what in order to determine eqgrp */
    /* qsort(dih,ndih,sizeof(struct dihedral),dihcmp_h); */
    eqga = eqgb + ndih;
    /* dipa is first measurement metadata record in the next chain */
    /* dipb is the current measurement metadata record in the chain */
    /* dipc is the previous measurement metadata record in the chain */
    for (maxdih = (dipa=dih) + ndih; dipa < maxdih;) {
	if (eqgb == eqga) {
	    fprintf(stderr,"ERROR: consistency error #5; please contact a developer!\n");
	    return 1;
	}
	cpa=(eqgb->first=dipc=dipa)->hash;
	m = dipa->tmult;
	t1=dipa->targ1;
	t2=dipa->targ2;
			 /* Need to do this here because both dipa and dipc may change */
	cpb=dipa->name;  /* before the following (inner) loop's last iteration. */
	eqgb->mult = PHASEFLAG;  /* can't happen so we use it as a sentinel value */
	eqgb->ph_undef = 1;  /* in the initial state, all phases are undefined */
	eqgb->merged=0;
	i=1;
	f=true;  /* no more unassigned measurements found */
	for (dipb=dipa+1; dipb<maxdih; dipb++)
	    if (dipb->next == dih) {  /* still using dih as a sentinel value */
    /* Next line OK because sprintf writes trailing \n that cause strncmp to stop */
		if (! strncmp(cpa,dipb->hash,DHASH)) {
		    if (m != dipb->tmult || (m != -1 && (t1 != dipb->targ1 || t2 != dipb->targ2))) {
			fprintf(stderr,"ERROR: different bias target for %s and %s !\n",
					cpb,dipb->name);
			return 1;
		    }
		    dipc->next=dipb;
		    dipc=dipb;
		    i++;
		} else if (f) {
		    dipa=dipb;
		    f=false;
	    }   }
        dipc->next=NULL;
	eqgb++->size=i;
	if (f) break;
    }
    negrp = (eqgmax=eqgb) - eqgrp;

    /*************************************************\
    * Ask the user to define extra equivalence groups *
    \*************************************************/
    if ((n=gopt(options,'q')) || interact) {  /* order matters because we need n */
	if (! (eqextra = malloc(sizeof(int) * interact ? MAXGREQ : 2))) {
	    fprintf(stderr,"ERROR: can't allocate memory for equivalence groups!\n");
	    return 1;
	}
	if (n) {
	    opts=options;
	    ipb  = ipa = NULL;     /* only to satisfy (ipa==ipb) below; overwritten */
	}  /* once n hits 0; NULL is most likely to make any logical mistake fatal. */
	for (;;) {  /* repeat to allow user to define multiple extra equivalences */
	    if (n) {
		for (;;) {
		    switch (opts->key) {
		      case '\0':  /* n is necessarily != 0 */
			fprintf(stderr,"ERROR: -q does not follow -p !\n");
			return 1;
		      case 'p':
			ccp=opts++->arg;
			if (opts->key == 'q') break;  /* ONLY way to exit loop alive */
			continue;
		      default:
			opts++;
			continue;
		    }
		    break;
		}
		if (! ((eqgb=findgrp(ccp,'p')) &&
			(eqgc=findgrp(opts++->arg,'q')))) return 1;
		if (eqgb == eqgc) {  /* will rarely be triggered */
		    /* In principle, we could safely let the user get away with a non-
		     * fatal warning and immediately trigger the next iteration of the
		     * outer loop with "continue". However, this condition is likely
		     * triggered by either a typo or a fundamental flaw in the user's
		     * logic, and in both cases, I feel their interest is better served
		     * with a fatal error that forces them to correct whatever's wrong. */
		    fprintf(stderr,"ERROR: equivalencing parameter '%s' with itself\n"
				    "not permitted!\n",ccp);
		    return 1;
		}
		if (eqgb < eqgc) {  /* somewhat less likey */
		    eqga=eqgb;
		    eqgb=eqgc;
		} else eqga=eqgc;  /* somewhat more likely */
		j=2;  /* Needs to be set every time because of --j below */
		n--;  /* TODO: once we have a verbosity flag, make it output the list
		       * of equivalence groups when n hits 0 */
	    } else {  /* next line is only checked repeatedly if we're interactive... */
		if (! interact) break;  /* out of outer for */
		printf("\nThe following equivalence groups are in effect:\n");
		/* Temporary abuse of eqga; gets overwritten later, before it matters */
		for (eqga=eqgrp; eqga<eqgmax; listeqg(stdout,"  e",eqga++,3));
		printf("Input space-separated list of equivalence groups to merge\n"
		        "or press ENTER to keep the above equivalence groups.\n");
		fflush(stdout);
		switch (j=get_list(line,INT,eqextra,2,MAXGREQ,1,negrp,
		                        "equivalence group","equivalence groups")) {
		  case -2:  /* specified more than MAXGREQ integers in MAXLINE characters... */
		    fprintf(stderr,"ERROR: consistency error #2; please contact a developer!\n");
		    /* ACHTUNG: fall through to case -3 */
		  case -3: return 1;  /* EOF */
		  case -1: continue;  /* insufficient input - retry */
		  case 0: goto exit_eq; /* user pressed ENTER. Harmless per K&R SECTION 3.8. */
		}
		ipb  = (ipa=eqextra) + j;
		eqga=&eqgrp[(*(ipa++))-1];
		eqgb=&eqgrp[(*(ipa++))-1];  /* yes, ++ again! */
	    }
	    /* apply equivalences - magic begins here */
	    /* eqextra: user-supplied list of equivalence groups to be merged */
	    /* dipa: first parameter of PREVIOUS equivalence group */
	    /* eqga: first equivalence group */
	    /* eqgb: "from" equilvalence group */
	    /* eqgc: "to" equivalence group */
	    /* eqgd: next equivalence group to be merged */
	    /* ipa: integer ID of NEXT equivalence group (after eqgd!) */
	    /* ipb: integer ID of maximum equivalence group */
	    m = (dipa = eqga->first)->tmult;
	    i = eqga->size + (eqgd=eqgc=eqgb)->size;
	    t1 = dipa->targ1; /* first can never be NULL */
	    t2 = dipa->targ2;
	    cpb = dipa->name;
	    eqga->merged=1;
	    for (; eqgb<eqgmax; eqgb++)
	        if (eqgb==eqgd) {  /* merger of 2 groups; triggered on first pass */
	            for (;;) if (dipb = dipa->next) dipa=dipb;
	                  else break;  /* out of inner for on previous line */
	            dipa->next = dipb = eqgd->first;
	            dipa = dipb;  /* first can never be NULL so we can skip a cycle */
	            if (m != dipb->tmult || (m != -1 && (t1 != dipb->targ1 || t2 != dipb->targ2))) {
	                fprintf(stderr,"ERROR: different bias target for %s and %s !\n",
	                                cpb,dipb->name);
	                return 1;  /* Last line: same error message is OK because */
	            }             /* stage is different and it mentions file names. */
	            if (ipa==ipb) eqgd=NULL;
	              else i += (eqgd=&eqgrp[(*(ipa++))-1])->size;
	        } else *(eqgc++) = *eqgb;  /* shift unaffected equivalence groups to avoid gaps */
	    eqga->size=i;
	    eqgmax = eqgrp + (negrp -= --j);
	}
	exit_eq:
	free(eqextra);
    }  /* next line: shrink; recalculate eqgmax for in case realloc() moves the block */
    eqgmax = (eqgrp = realloc(eqgrp,sizeof(struct eqgroup) * negrp)) + negrp;

    /*****************************************************************************\
    * Ask the user for a null bias. Zero is allowed (but not really recommended). *
    \*****************************************************************************/
    switch ((i=gopt(options,'u')) + (j=gopt(options,'t')) + gopt(options,'6')) {
      case 0:
	if (interact) for (;;) {  /* allow the user to retry when giving invalid input */
		printf("Biasing mode (1 = uniform (default), 2 = adaptive, "
			"3 = abs (experimental)): ");
		fflush(stdout);
		switch (get_list(line,INT,&bmod,1,1,1,3,"biasing mode","biasing modes")) {
		  case -3: return 1;  /* EOF */
		  case -2: case -1: continue;  /* retry; only way to repeat the loop */
		  case 0:
		   bmod=UNIFORM;
		   break;  /* out of switch for good measure */
		}
		break;
	    }
	  else bmod=UNIFORM;
	break;  /* out of switch */
      case 1:
	if (i) bmod=UNIFORM;
	  else if (j) bmod=ADAPT;
	  else bmod=ABS;
	break;  /* out of switch */
      default:
	fprintf(stderr,"ERROR: options -u, -t and --abs are mutually exclusive!\n");
	return 1;
    }
    if ((ccp=outsideopt(4,'b',"global bias fraction")) == &empty) return 1;
    if (ccp) {
	if (! sscanf(ccp,"%lf",&biasfract)) {
	    fprintf(stderr,"ERROR: argument --glob-bias has invalid format!\n");
	    return 1;
	}
    } else if (interact) for (;;) {  /* allow the user to retry when giving invalid input */
	    /* We can't fill in a default value if the user presses ENTER because in a later
	     * version, ENTER will mean that the user wants to fill in per-multiplicity biases,
	     * and we don't want to break compatibility at that point. For now, we just
	     * re-prompt on ENTER. When we implememt per-multiplicity biases, we can insert a
	     * case 0: biaspermult=true; break;
	     * TODO: the above discussion becomes invalid if we decide to only allow
	     * per-multiplicity biases on the command-line, so ENTER may come to mean
	     * "default" after all. */
	    /* TODO: change terminology ("null bias" --> "bias fraction") in comments to avoid
	     * confusing future maintainers. */
	    printf("Bias fraction? (default: %g) ",DEFBIAS);
	    fflush(stdout);
	    /* 0,0,"" are arguments that are only used if kind == INT */
	    switch (get_list(line,DOUBL,&biasfract,1,1,0,0,"","bias fraction")) {
	      case -3: return 1;  /* EOF */
	      case -2: case -1: continue;  /* retry; only way to repeat the loop */
	      case 0:
	       biasfract = DEFBIAS;
	       break;  /* out of switch for good measure */
	    }
	    break;
	}
      else biasfract = DEFBIAS;
    /* TODO: validation of the actual numerical value (i.e. 0 < biasfract < 1)? */
    if (gopt(options,7)) biascomp=NO;
      else if (interact) {
	/* TODO: check if all interactive defaults actually work (some seem to not) */
        printf("Enable bias compensation? (recommended; default: \"yes\") ");
        fflush(stdout);
        if ((biascomp=yesno(line,"bias compensation",YES)) == YNVOID) return 1;
    } else biascomp=YES;
    if (interact) {
	/* TODO: change to the following when implementing fixed-phase angles */
	/* printf("Use default phases? (recommended; not available for valence angles) "); */
	printf("Use default phases? (recommended; default: \"yes\") ");
	fflush(stdout);
	if ((phasedef=yesno(line,"variability of phases",YES)) == YNVOID) return 1;
    }
    /***************************************\
    * Ask the user to choose multiplicities *
    \***************************************/
    /* TODO LATER: once we allow the user to assign per-multiplicity biases, we'll need
     * to allocate this here to store their input prior to transfering it to biasvec . */
    /* if (! (biastmp = malloc(sizeof(double) * negrp * 7))) {
     *     fprintf(stderr,"ERROR: can't allocate temporary storage for bias factors!\n");
     *     return 1;
     * } */
    np=0;
    ccp=NULL;
    eqga=NULL;
    for (opts=options;;opts++) {
	switch (v=opts->key) {
	  case 'q':
	    if (eqga) {  /* should have triggered "-q does not follow -p" above */
		fprintf(stderr,"ERROR: consistency error #6; "
				"please contact a developer!\n");
		return 1;
	    }
	    ccp=NULL;
	    continue;
	  case 'p':
	    if (eqga) {
		finalize_mult(eqga,multemp,i,pf);
		eqga=NULL;
	    }
	    ccp=opts->arg;
	    continue;
	  case 'm':
	    if (! ccp) {
		fprintf(stderr,"ERROR: -m flag outside parameter section!\n");
		return 1;
	    }
	    /* TODO LATER: if we allow for multiple multiplicity sections per parameter,
	     * it may be possible to simply replace the two following if conditions with
	     * these two lines, but it's far from sure; see the discussion above the
	     * actual finalize_mult() function. */
	    /* if (eqga) finalize_mult(eqga,multemp,i,pf);
	     *   else if (! (eqga=findgrp(ccp,'p'))) return 1; */
	    if (eqga) {
		fprintf(stderr,"ERROR: only one -m flag currently allowed "
				"per parameter section!\n");
		return 1;
	    }
	    if (! (eqga=findgrp(ccp,'p'))) return 1;
	    ipa=multemp;
	    for (comm=opts->arg;;) {  /* abuse comm for a couple of lines */
		if (! (u=*comm++)) break;
		if (u >= '0' && u <= '6') *ipa++ = u - '0';
		  else if (u != ' ') fprintf(stderr,"Warning: extraneous character "
				"'%c' in multiplicity specification ignored.\n",u);
	    }  /* Previous line: this is pretty much the only command-line parsing issue
	       * that I decided to make nonfatal, to allow for some freedom in the
	      * specification ("1-2-3-4-6", "1/2/3/4/6", tab characters,...) Still, I
	     * wonder if I won't regret later that I exposed the user to this risk. Also
	     * TODO LATER: verbosity. */
	    if (! (i=ipa-multemp)) {
		fprintf(stderr,"ERROR: argument -m has invalid format!\n");
		return 1;
	    }
	    qsort(multemp,i,sizeof(int),intcmp);
	    pf = PH_UNDEF;
	    continue;
	  case 'f': case 'a':  /* case 'w': case 'b': */ /* see "TODO LATER" below */
	    if (! eqga) {
		fprintf(stderr,"ERROR: -%c flag outside multiplicity section!\n",v);
		return 1;
	    }
	    switch (v) {
	      case 'f':
	        pf = PH_FIX;
	        continue;
	      case 'a':
	        pf = PH_VAR;
	        continue;
	    /* TODO LATER: the following cases can be commented in once we start
	     * allowing per-multiplicity weight factors, resp. biases. Don't forget
	     * to adjust the built-in help too; an old version of it already
	     * foresaw the introduction of these features. */
	    /*case 'w':
	     *  continue;
	     *case 'b':
	     *  continue; */
	    }  /* next line: very blatantly impossible, but the */
	       /* switch would look very unclean without it. */
	    fprintf(stderr,"ERROR: consistency error #10; "
			    "please contact a developer!\n");
	    return 1;
	  case '\0':
	    if (eqga) finalize_mult(eqga,multemp,i,pf);
	    break;  /* out of switch and for */
	  default:
	   /* Again, we could make this (as well as some of the other
	    * command-line parsing errors we've thrown up to this point)
	    * nonfatal, but as discussed for the "equivalencing parameter with
	    * itself" error, that isn't really in the user's best interest. */
	    if (ccp) {
		fprintf(stderr,"ERROR: extraneous flag -%c in %s section!\n",
				v, eqga ? "multiplicity" : "parameter");
		return 1;
	    }
	    continue;
	}
	break;
    }
    gopt_free(options);
    /* TODO LATER: this loop hasn't reached critical mess yet, but it's getting there,
     * with a control flow that's starting to look like spaghetti and some duplicated
     * branching. If we ever get to refactoring this program (which would definitely be
     * of utility), then this would presumably become one of the priorities. */
    maxeq0size=0;
    for (eqga = eqgrp; eqga < eqgmax; eqga++) {
	u=eqga->mult;
	if (! (dipb=eqga->first)->hash[LASTATT]) {  /* not a dihedral */
	    np+=2;  /* 2 harmonics */
	    if (u != BANG) {  /* TODO: make finer-grained when implementing fixed-phase angles */ 
		/* Yet again, we could very easily make this "Warning: multiplicities
		 * on non-dihedral equivalence group %i (%s etc.) reset to 0\n", but that
		 * just increases the risk a fatal typo will go undetected. */
		if (u != PHASEFLAG) {
		    listeqg(stderr,"ERROR: nonzero multiplicity or fixed phase"
				    "on non-dihedral\n  e",eqga,3);
		    return 1;
		}
		eqga->mult = BANG;
	    }
	    if ((l=eqga->size) > maxeq0size) maxeq0size=l;
	    continue;
	}
	/* TODO LATER: this may get a lot more difficult once we have per-multiplicity
	 * biases and/or weights */
	if (! eqga->ph_undef) {
	    /* If the phase were defined on the command-line, its definition must
	     * necessarily have been part of a valid multiplicity section, which
	     * would have set eqga->mult to a different value than PHASEFLAG, so: */
	    if (u == PHASEFLAG) {  /* sentinel value meaning undefined */
		fprintf(stderr,"ERROR: consistency error #11; "
			    "please contact a developer!\n");
		return 1;
	    }
	    if (u & 1 && (l=eqga->size) > maxeq0size) maxeq0size=l;
	    continue;
	}
	listeqg (stdout,"\nE",eqga,1);
	if (u == PHASEFLAG) {  /* sentinel value meaning undefined */
	    if (! interact) {
		fprintf(stderr,"ERROR: missing multiplicity specification "
				"for dihedral parameter!\n");
		return 1;
	    }
	    for (;;) {  /* allow the user to retry when giving invalid input */
		printf("Input space-separated list of multiplicities "
		    "for this equivalence group\n"
		    "or press ENTER to deactivate the equivalence group.\n");
		fflush(stdout);
		switch (j=get_list(line,INT,multemp,1,7,0,6,
					"multiplicity","multiplicities")) {
		  case -3: return 1;  /* EOF */
		  case -2: case -1: continue;  /* retry; only way to repeat the loop */
		}
		break;  /* (j==0) = ENTER or (j>0) = valid multiplicities */
	    }
	    if (! j) {  /* ENTER: deactivate eq group. Not possible from command-line
			   and not really good practice in interactive interface; why
			   specify measurements for this eq group in the first place? */
		eqga->mult = '\0';
		continue;  /* can only "continue" outer loop outside inner loop above */
	}   }
	f=false;
	if (dipb->tmult != -1 && dipb->targ2) {
	    /* Unusually verbose warning because this is a very serious pitfall */
	    fprintf(stderr,
	        "WARNING: phase on this equivalence group automatically allowed to vary based\n"
	        "on odd target phase. If this is not what you expected, then the first line\n"
	        "in %s is probably wrong!\n",dipb->name);
	    f=true;
	} else if (interact && ! phasedef) {
	    printf("Allow phase to vary? (default: \"no\") ");
	    fflush(stdout);
	    switch (yesno(line,"phase variability",NO)) {
	      case YNVOID: return 1;
	      case YES: f=true;
	}   }
	if (u != PHASEFLAG) {     /* If u were undefined, we'd have to transfer the */
	    if (f) u |= PHASEFLAG; /* new multiplicities from multemp to eqga->mult */
	    continue;
	}  /* Next line: order matters! */
	if ((eqga->mult=packmult(multemp,j,f)) & 1 &&
				    (l=eqga->size) > maxeq0size) maxeq0size=l;
    }

    /*************************\
    * Actual math starts here *
    \*************************/
    /* nt on next line: sic; parvec will be fed as I/O to dgels */
    if (! (parvec = malloc(sizeof(double) * (nt = ns + np)))) {
	fprintf(stderr,"ERROR: can't allocate memory for parameter vector!\n");
	return 1;
    }
    if (! (cosmat = malloc((scansz = sizeof(double) * ns) * np))) {
	fprintf(stderr,"ERROR: can't allocate memory for cosine functions!\n");
	return 1;
    }
    /* The following is mostly guaranteed only to overflow if the architecture itself
     * can't handle it, because we start with sizeof, which is of type size_t */
    if (! (amat = malloc(sizeof(double) * nt * np))) {
	fprintf(stderr,"ERROR: can't allocate memory for linear system!\n");
	return 1;
    }
    lwork=-1;  /* workspace query */
    dgels_(&notrans,&nt,&np,&one,amat,&nt,targvec,&nt,&d,&lwork,&info);
    if (info) {  /* Workspace query cannot fail */
	fprintf(stderr,"ERROR: consistency error #7; please contact a developer!\n");
	return 1;
    }
    lwork = (fint) d;
    if (! (work = malloc(sizeof(double) * lwork))) {
	fprintf(stderr,"ERROR: can't allocate workspace for solving linear system!\n");
        return 1;
    }
    /* TODO: This should be moved once we have per-multiplicity null biases. */
    if (! (nbvec = malloc(sizeof(double) * np))) {
	fprintf(stderr,"ERROR: can't allocate memory for biasing restraints!\n");
	return 1;
    }
    if (! (biasvec = malloc(sizeof(double) * np))) {
	fprintf(stderr,"ERROR: can't allocate memory for bias factors!\n");
	return 1;
    }
    if (maxeq0size)
        if (! (phasetmp = malloc(sizeof(double) * ns * maxeq0size))) {  /* float in "cosf" branch */
            fprintf(stderr,"ERROR: can't allocate memory for range check!\n");
            return 1;
        }
    printf("Calculating cosines...\n");
    fflush(stdout);
    /* Just a quick summary: whenever there's a 0-fold present, eqga->delta1 and delta2
     * are set by scanning the range. However, they are ignored by all other
     * multiplicities, ie. the phase of cosine functions is determined by da and db. */
    i=1;
    dpc = cosmat;
    dpd = biasvec;
    /* TODO LATER: once we allow the user to assign per-multiplicity biases... */
    /* dpf = biastmp; */
    /* This program loops over all the multiplicities in all the equivalence groups many
     * times. To curb this proliferation of boilerplate code, we could make a
     * loopovermult() function that takes a pointer to a second "inner loop" function as
     * an argument. However, doing so would significantly increase complexity in the
     * sense that many arguments need to be passed to the inner loop function if we want
     * to avoid making a whole bunch of variables global or using no-so-portable nested
     * functions / lexical scoping. In the end, the boilerplate code is not so bad
     * because it's only a few lines, so for now, we honor KISS over modularity. The
     * truly evil solution would be to use macros - muhahaha! */
    for (eqga = eqgrp; eqga < eqgmax; eqga++) {
	g = (! ((u=eqga->mult) & PHASEFLAG));
	b = ((cpa=(dipb=eqga->first)->hash)[ANGATT] == '\0');  /* it's a bond */
	k = dipb->tmult;  /* In the block of nested if()s, -1 means no target, but after */
	    /* the loop over j, it can also mean that the target matches a multiplicity. */
	/* TODO: the following nested block of if() conditions has grown organically and
	 * has now reached critical mess; therefore, a clean re-implementation is
	 * desired. See ROADMAP for details. */
	if (h = (cpa[LASTATT] != '\0')) {  /* it's a dihedral */
	    if (u & 1) {       /* if there's a 0-fold...                            */
		dpe=phasetmp;  /* ...then search for the largest gap in phase space */
		for (dipa=dipb;dipa;dipa=dipa->next)
		    /* in the cosf branch, we need a float equivalent for scansz */
		    /* memcpy(dpe,dihmat+(dipa-dih)*ns,scansz);
		     *dpe+=ns; */
		    for (dpb = (dpa=dihmat+(dipa-dih)*ns) + ns; dpa < dpb;
		                    *(dpe++) = fmod(*(dpa++)+540.0,360.0) - 180.0);
		qsort(phasetmp,dpe-phasetmp,sizeof(double),floatcmp);
		delta = (da=d1=*phasetmp) + 360.0 - (d2=*(dpe-1));
		for (dpa=phasetmp+1;dpa<dpe;) {
		    if ((d = (db=*(dpa++)) - da) > delta) {
			delta=d;
			d1=db;
			d2=da;
		    }
		    da=db;
		}
		if (delta > 360.0-MINRANGE) {  /* TODO LATER: should be finer-grained */
		    fprintf(stderr,"ERROR: measurement range too small for "
		                                    "equivalence group %i !\n",i);
		    return 1;
		}
		if (delta < 360.0-MAXRANGE) {
		    fprintf(stderr,"ERROR: measurement range too large for fitting harmonic potential\n"
		                                    "to equivalence group %i !\n",i);
		    return 1;
		}
		da = fmod(d2+360.0,360.0) - 180;  /* da,db: complimentary points (where cusps are) */
		db = fmod(d1+360.0,360.0) - 180;
		if (db < da) db += 360.0;     /* else range check will go wrong */
		if (delta < 180.0) {
		    if (!g) fprintf(stderr,"warning: range > 180deg for fitting harmonic potential "
		            "with variable phase to\nequivalence group %i. Make sure to understand "
		            "the issues!\n",i); /* The null bias works against phase fitting when d1 */
		      /* and d2 are close together though it doesn't affect the resultant amplitude. */
		    d1=da;  /* We have to use the complimentary points for the two    */
		    d2=db;  /* "basis functions" but they also still count as limits! */
		} else if (d2 < d1) d2 += 360.0;  /* else addition will go wrong */
		if (! k) {  /* there is a target and it's a 0-fold */
		    /* Shift ranges so that they're centered around the target. The first
		     * line contains a redundancy (d1 occurs 2 times; this is not the case
		     * for the alternative method we're using below for bonds and angles)
		     * but saves computations by the time we get to da, db . That said,
		     * there remain redundancies between the stuff immediately above and
		     * below this point and simplification is likely possible (TODO!) */
		    d1 = fmod(d1+(d=dipb->targ2-(d2+d1)/2)+540,360) - 180;
		    d2 = fmod(d2+d+540,360) - 180;
		    da = fmod(da+d+540,360) - 180;
		    db = fmod(db+d+540,360) - 180;
		    if (d2 < d1) d2 += 360.0;  /* else addition will go wrong */
		    if (db < da) db += 360.0;  /* else range check will go wrong */
		}
		eqga->delta1 = d1;  /* Won't work for cosines because phasor addition */
		eqga->delta2 = d2;  /* implementation assumes orthogonality.          */
		eqga->dmin = da;    /* Therefore, cosines must ignore d1, d2 .        */
		eqga->dmax = db;
	    } else {
		eqga->dmax = eqga->dmin = 0.0;
		if (u > BANG) {  /* there are cosines AND phase is variable */
		    if (k > 0) {  /* the target is also a cosine function */
			db = (da = -dipb->targ2) + 45.0;  /* Immediate reuse of da and db */
			da -= 45.0;
		    } else {
			db = 45.0;
			da = -45.0;
		    }
		} else db = 0.0;  /* da isn't used in this case */
	    }
	} else {  /* it's a bond or angle (phase automatically variable) */
	    for (dipa=dipb;dipa;dipa=dipa->next) {
		d1 = HUGE_VAL;
		d2 = -HUGE_VAL;
		for (dpb = (dpa=dihmat+(dipa-dih)*ns) + ns; dpa < dpb; ) {
		    d = *(dpa++);
		    if (d < d1) d1=d;
		    if (d > d2) d2=d;
	    }   }
	    if (k > 0) {  /* bond or angle with multiplicity 1 or greater */
		fprintf(stderr,"ERROR: consistency error #9; please contact a developer!\n");
		return 1;  /* TODO: actually, this blatantly can't happen, and we have */
	    }              /* not been very consistent about issuing these errors... */
	    if (!k) {  /* multiplicity 0 ==> there is a target */
		/* Shift range so that it's centered around the target. The way we're */
		/* doing it here is mathematically equivalent to what we're doing with
		 * 0-fold dihedrals above, but with less redundant calculations. */
		d1 = (e=dipb->targ2) - (d=(d2-d1)/2);  /* target - half the distance */
		d2 = e + d;                            /* target + half the distance */
		da = e - (d *= PRESRANGE);
		db = e + d;
#ifdef DEBUG
	        fprintf(stderr,"DEBUG: d1 = %g ; d2 = %g\n",d1,d2);
#endif
	    } else {
		/* Even if there is no target, the correctly fitted MM phase is often
		 * outside the scan range centered around the QM minimum, owing to the
		 * other forces in the molecule. Hence EXTRANGE , which actually makes
		 * more sense than the PRESRANGE feature above. */
		da = d1 - (d = EXTRANGE * (d2-d1));
		db = d2 + d;
	    }
	    /* TODO: the following can probably be put outside of the "if" condition once
	     * we've completed the planned overhaul of the range checking. */
	    eqga->delta1 = d1;
	    eqga->delta2 = d2;
	    eqga->dmin = da;
	    eqga->dmax = db;
	}
	for (j=0;j<NMULT;j++) {  /* Actual cosine calculation */
	    if (u & 1) {
		if (j == k) k = -1;  /* there's a multiplicity matching the target */
		for (f=false;;f=true) {
		    for (dpb = (dpa=dpc) + ns; dpa < dpb; *dpa++ = 0.0);
		    for (dipa=dipb;dipa;dipa=dipa->next) {
			curdih=dihmat+(dipa-dih)*ns;  /* apply cosine function */
			              /* More performant and easier to vectorize than */
			switch (j) {  /* if the switch were inside the inner loop.    */
			  case 0:  /* harmonic function */
			    if (f) delta=-d2;
			      else delta=-d1;
			    if (b) for (dpa = dpc; dpa < dpb; ) {
				    d = *(curdih++) + delta;
				    /* Faster fix for bond problem if there are few eq groups */
			            /* *(dpa++) += d * d * BONDF; */
			            *(dpa++) += d * d;
				}
			      else for (dpa = dpc; dpa < dpb; ) {  /* 0-fold always has delta */
				    /* Confirmed in source that CHARMM maps to [-180:180] AFTER
				     * applying delta. Possibly CHARMM-specific. Performing the mapping
				     * BEFORE will give trouble with incorrectly defined impropers (we
				     * now know that this is not the cause of the "tie fighter" and
				     * that said incorrectly defined impropers are stable under typical
				     * MD conditions, though not under SGLD) and could hypothetically
				     * give trouble with linear (sp) angles if negative measurements
				     * are allowed (unlikely). The only advantage would be the fitting
				     * of 0-folds to 360deg dihedral scans, but the combination of the
				     * 1-fold periodicity and the discontinuity should make 0-folds
				     * unsuitable both for single and double bonds. */
			            /* d = (*(curdih++) + delta) * DEG_F; */
			            d = (fmod(*(curdih++) + delta + 540.0, 360.0) - 180.0) * DEG_F;
				    /* Faster fix for bond problem if there are few eq groups */
			            /* *(dpa++) += d * d * ANGIMPF; */
			            *(dpa++) += d * d;
			        }
			    break;  /* out of switch */
			  case 1:  /* save a multiplication by 1 */
			    if (f)
			        for (dpa = dpc; dpa < dpb; *(dpa++) += cos((*(curdih++)+da)*DEG_F));
			      else for (dpa = dpc; dpa < dpb; *(dpa++) += cos((*(curdih++)+db)*DEG_F));
			    break;  /* out of switch */
			  default:
			    if (f)
			        for (dpa = dpc; dpa < dpb; *(dpa++) += cos((*(curdih++)*j+da)*DEG_F));
			      else for (dpa = dpc; dpa < dpb; *(dpa++) += cos((*(curdih++)*j+db)*DEG_F));
			    break;  /* out of switch for good measure */
		    }   }
		    /* Faster fix for bond problem if there are many eq groups. Need to crawl through
		     * the array one additional time, but not doing redundant multiplications. */
		    if (!j) {
		        if (b) x = BONDF;
		          else x = ANGIMPF;
		        for (dpa = dpc; dpa < dpb; *(dpa++) *= x);
		    }
		    /* In the case of varying phases, doing this here will lead to a different
		     * "normalization" (more correctly alignment) for the components, ie. a
		     * different constant being added to each energy function. The offset of the
		     * resultant energy function will similarly be a linear combination of the
		     * offsets of the components, without affecting the non-constant part of the
		     * energy function. Also, when applying this correctly, the resultant is still
		     * "normalized" (ie. RMSD aligned, ie. its integral is 0). If we were to
		     * attempt to align the components *together*, the non-constant part would
		     * still be sound, but the resultant would generally *not* be aligned, so that
		     * would be worse. The same thing (correct non-constant part but alignment
		     * lost) will also inevitably happen when implementing zero-alignment (the
		     * resultant integral is still the weighted sum of the component integrals, but
		     * unless all the minima are at the same location, no such simple statement can
		     * be made about the value at the absolute minimum). That's why it's generally
		     * nontrivial to implement the zero aligned least-squares fitting (it would
		     * require multiple fits, each including extra equations that ensure the
		     * alignment and that are exactly satisfied using DGGLSE (or conceivably
		     * DGGGLM)). Anyhow TODO: change terminology (normalize --> align) in code to
		     * avoid confusing future maintainers. */
		    normalize_weight (dpc);
		    dpc=dpb;
		    /* TODO LATER: once we allow the user to assign per-multiplicity biases, we
		     * need to do this selection earlier, and just have the following line here: */
		    /* *dpd++ = dpf[j]; */
		    /* Dihedrals always get biased - bonds and angles only if they have a target */
		    *dpd++ = (h || dipb->tmult == j) ? biasfract : 0.0;
		    if (f || (j && g)) break;  /* out of for */
	    }   }
	    u >>= 1;
	}
	if (k != -1) {
	    fprintf(stderr,"ERROR: %i-fold in %s but not among multiplicities for "
			    "equivalence group %i !\n",k,dipb->name,i);
	}
	/* TODO LATER: once we allow the user to assign per-multiplicity biases... */
	/* dpf += 7; */
	i++;
    }
    if (maxeq0size) free (phasetmp);
    /* TODO LATER: once we allow the user to assign per-multiplicity biases... */
    /* free(biastmp);  /* So far, this appears to be the first instance where
                        * we would be fragmenting memory. */
    printf("Applying biasing restraint...\n");
    fflush(stdout);
    /* Re-use phasetmp . Won't get away with this if we move it someplace else. Also,
     * phasetmp is float in "cosf" branch which may or may not be appropriate for this
     * purpose (haven't thought this out). */
    if (! (phasetmp = malloc(sizeof(double) * ((np * (npp=np+1)) >> 1)))) {
        fprintf(stderr,"ERROR: can't allocate memory for storing <f|g>!\n");
        return 1;
    }
    if (bmod != UNIFORM)
	if (! (nbtmp = malloc(sizeof(double) * np))) {
	    fprintf(stderr,"ERROR: can't allocate memory for storing <f|t>!\n");
	    return 1;
	}
    if (np & 1) {  /* np odd */
	hdim = npp >> 1;
#ifdef DEBUG
	hdimm = hdim - 1;
#endif
	o = 0;
	np2 = np;
	m = -((hdim-1) * np + hdim);
    } else {      /* np even */
	hdim = np >> 1;
#ifdef DEBUG
	hdimm = 0;  /* we can only use hdimm==0 as a parity test if np is guaranteed */
#endif             /* > 1, which is not the case in this program, so we use o instead. */
	o = 1;
	np2 = npp;
	m = -(hdim * npp + hdim);
    }
    /* No need to initialize phasetmp as long as beta=0.0 (LAPACK has a branch for it) */
    dsfrk_(&notrans,&lo,&trans,&np,&ns,&alpha,cosmat,&ns,&beta,phasetmp);
#ifdef DEBUG
    fprintf(stderr,"DEBUG: R.R =\n");
    for (i=0;i<np;i++) {
	for (k=0;k<i;k++) fprintf(stderr,"         ");
	for (;k<np;k++) if (o)  /* np even */
		fprintf(stderr,"%8.3g ",phasetmp[i<hdim ? i*npp+k+1 : (k-hdim)*npp+i-hdim]);
	      else               /* np odd */
		fprintf(stderr,"%8.3g ",phasetmp[i<hdim ? i*np+k : (k-hdimm)*np+i-hdim]);
	putc('\n',stderr);
    }
#endif
    if (bmod != UNIFORM) {
        dgemv_(&trans,&ns,&np,&alpha,cosmat,&ns,targvec,&one,&beta,nbtmp,&one);
	dpg=nbtmp;
    }
    dpc=cosmat;
    dpd=amat;
    dpe=nbvec;
    dph=biasvec;
    i=0;  /* in the paper, we called this k for the purpose of calulating b_k */
    p=1;  /* This counts the equivalence groups for the purpose of issuing warnings */
    for (eqga = eqgrp; eqga < eqgmax; eqga++) {
	g = (((u=eqga->mult) & PHASEFLAG) == 0);  /* == 0 is poor man's conversion to bool */
	for (j=0;j<NMULT;j++) {
	    if (u & 1) for (f=false;;f=true) {
		    memcpy(dpd,dpc,scansz);
		    for (dpb = (dpa=(dpd+=ns)) + np; dpa < dpb; *dpa++ = 0.0);
		    sum=0.0;
		    switch (bmod) {
		      case ADAPT:
			dpf=nbtmp;
			x = *dpg++;
			break;  /* out of switch */
		      case ABS:
			dpf=nbtmp;
			x = fabs(*dpg++);
			break;  /* out of switch for good measure */
		    }
		    dpi=biasvec;
		    if (h = j && g) n = -1;  /* assignment, not == */
                      else if (f) n = i - 1;
		      else n = i + 1;
		    l = i + o;
		    switch (bmod) {  /* Keep more branching out of inner loop. The code duplication is hideous, But the */
		      case ADAPT:   /* alternative would be even worse, with multiple conditions in the inner loop. This
				   /* could be worked around by setting a pointer to a one-line function, but the latter
				  /* would rely on a bunch of local variables, and copying them to the stack in the inner
				 /* loop would be terrible for performance. The only winning move is not to play. */
			if (i < hdim) {  /* Pray that branch prediction is smart enough... */
			    /* we can get away with k++ because K&R 2.11 (p.51): "only one of expr2 and expr3 is evaluated" */
			    for (k=0;k<i;dpf++)    sum += (k==n ? x : *dpf) * *dpi++ * phasetmp[k++ * np2 + l];
			    l = i * np2 + o;
			    for (   ;k<np;dpf++)   sum += (k==n ? x : *dpf) * *dpi++ * phasetmp[l + k++];
			} else {
			    for (k=0;k<hdim;dpf++) sum += (k==n ? x : *dpf) * *dpi++ * phasetmp[k++ * np2 + l];
			    l = i * np2 + m;
			    for (   ;k<i;dpf++)    sum += (k==n ? x : *dpf) * *dpi++ * phasetmp[l + k++];
			    l = m + i;
			    for (   ;k<np;dpf++)   sum += (k==n ? x : *dpf) * *dpi++ * phasetmp[k++ * np2 + l];
			}  /* Phew! Now I wish I had gone with dsprktest.c instead. */
#ifdef DEBUG
			fprintf(stderr,"DEBUG: sum = %17.12g\n",sum/x);
#endif
			/* Multiplication of denominator by (1.0 - *dph++) is not expensive
			 * enough or repeated often enough to justify precomputing the values. */
			dpd[i++] = *dpe++ = sqrt(fabs(sum / (x * (1.0 - *dph++))));
			break;  /* out of switch */
		      case ABS:
			if (i < hdim) {
			    /* TODO: The way the fabs()s were introduced is an ugly hack; if they
			     * really improve things, we may want to refactor the code a bit more. */
			    for (k=0;k<i;dpf++)    sum += *dpi++ * (k==n ? x * phasetmp[k++ * np2 + l] : fabs(*dpf * phasetmp[k++ * np2 + l]));
			    l = i * np2 + o;
			    for (   ;k<np;dpf++)   sum += *dpi++ * (k==n ? x * phasetmp[l + k++]       : fabs(*dpf * phasetmp[l + k++]));
			} else {
			    for (k=0;k<hdim;dpf++) sum += *dpi++ * (k==n ? x * phasetmp[k++ * np2 + l] : fabs(*dpf * phasetmp[k++ * np2 + l]));
			    l = i * np2 + m;
			    for (   ;k<i;dpf++)    sum += *dpi++ * (k==n ? x * phasetmp[l + k++]       : fabs(*dpf * phasetmp[l + k++]));
			    l = m + i;
			    for (   ;k<np;dpf++)   sum += *dpi++ * (k==n ? x * phasetmp[k++ * np2 + l] : fabs(*dpf * phasetmp[k++ * np2 + l]));
			}
			/* Yes, we still need fabs() here, because phasetmp[] (and hence x * phasetmp[]) can be negative,
			 * and we cannot eliminate this possibility without completely messing up the harmonics. */
			dpd[i++] = *dpe++ = sqrt(fabs(sum / (x * (1.0 - *dph++))));
			break;  /* out of switch */
		      case UNIFORM:       /* UNIFORM without ABS makes no sense as explained in the commit log for CVS rev. */
			if (i < hdim) {   /* 1.56 ; indeed, it was tried in rev. 1.55 and the harmonics were wrong. */
			    for (k=0;k<i;)    sum += *dpi++ * (k==n ? phasetmp[k++ * np2 + l] : fabs(phasetmp[k++ * np2 + l]));
			    l = i * np2 + o;
			    for (   ;k<np;)   sum += *dpi++ * (k==n ? phasetmp[l + k++]       : fabs(phasetmp[l + k++]));
			} else {
			    for (k=0;k<hdim;) sum += *dpi++ * (k==n ? phasetmp[k++ * np2 + l] : fabs(phasetmp[k++ * np2 + l]));
			    l = i * np2 + m;
			    for (   ;k<i;)    sum += *dpi++ * (k==n ? phasetmp[l + k++]       : fabs(phasetmp[l + k++]));
			    l = m + i;
			    for (   ;k<np;)   sum += *dpi++ * (k==n ? phasetmp[k++ * np2 + l] : fabs(phasetmp[k++ * np2 + l]));
			}
			/* TODO: it's a bit inconsistent to have this warning and debug
			 * message for the uniform bias but not for the other ones... */
			/* TODO: this warning is OK while we're still beta testing, but it won't
			 * mean much to non-expert users, so we want to hide it behind a -v flag. */
	                if (sum < 0.0) fprintf(stderr,"Warning: b_k for %s multiplicity %i on equivalence group %i\n"
				"is square root of opposite of negative number.\n",f ? "2nd" : "1st",j,p);
#ifdef DEBUG
			fprintf(stderr,"DEBUG: sum = %17.12g\n",sum);
#endif
			dpd[i++] = *dpe++ = sqrt(fabs(sum / (1.0 - *dph++)));
			break;  /* out of switch for good measure */
		    }
		    dpd=dpb;
		    dpc+=ns;
		    if (f || h) break;
		}
	    u >>= 1;
	}
	p++;
    }
    if (bmod != UNIFORM) free(nbtmp);
    free(phasetmp);
    printf("Solving linear system...\n");
    fflush(stdout);
    /* See also http://www.nag.co.uk/lapack-ex/node44.html */
    memcpy(parvec,targvec,scansz);  /* Next line: target for null bias should be zero. */
    for (dpb = (dpa=parvec+ns) + np; dpa < dpb; *dpa++ = 0.0);
    dgels_(&notrans,&nt,&np,&one,amat,&nt,parvec,&nt,work,&lwork,&info);
    if (info) {
	fprintf(stderr,"ERROR: failed to solve linear system!\n");
	return 1;
    }  /* solvec is a now-unused row of amat and will be used to temporarily
        * store the potential due to the parmeters that are frozen out. */
    solvec = amat + nt * (np-1);  
    printf("Verifying harmonic functions...\n");
    fflush(stdout);
    np2 = 0;
    /* TODO: consolidation of generic variables */
    dpf = parvec;
    for (dph = (dpa=solvec) + ns; dpa < dph; *dpa++ = 0.0);
    i=1;
    for (eqga = eqgrp; eqga < eqgmax; eqga++) {
	k = (dipb=eqga->first)->tmult;
	t1 = dipb->targ1;
	g=(((u=eqga->mult) & PHASEFLAG) == 0);
	for (j=0;j<NMULT;j++) {  /* Actual cosine calculation */
	    if (u & 1) {
		if (j) {  /* cosine function */
		    if (g) dpf++;
		      else dpf+=2;
		} else if (dipb->hash[ANGATT]) {  /* harmonic and not a bond */
		    d=*(dpf++);  /* d, e : force constants / amplitudes */
		    e=*(dpf++);  /* delta1, delta2 not in radians! */
		    delta = (d*eqga->delta1+e*eqga->delta2)/(x=d+e);
#ifdef DEBUG
		    fprintf(stderr,"DEBUG: k1 = %g ; k2 = %g ; k = %g ; delta = %g\n",d,e,x,delta);
#endif
		    f = true;  /* enables realignment */
		    if (delta < eqga->dmin || delta > eqga->dmax) {  /* out of range */
			/* If we're here, we know that j == 0. If k == 0 (not higher or
			 * lower) && targ1 != 0, then we have a target for the harmonic, so
			 * any reference value is possible and should only be a warning. */
			if (! k && t1) fprintf(stderr,"Warning: fitted reference value for "
			            "equivalence group %i out of acceptable range.\n",i);
			  else if (dipb->hash[LASTATT]) {  /* a dihedral */
			    fprintf(stderr,"WARNING: fitted reference value for "
			            "equivalence group %i out of acceptable\nrange; "
			            "refitting without harmonic potential.\n",i);
			    /* trigger warning in the prm file & refit */
			    eqga->dmax = eqga->dmin = HUGE_VAL;
			    np2 += 2;  /* 2 more degrees of freedom that fall out */
			    eqga->delta2 = eqga->delta1 = 0.0;
			    f = false;  /* do not realign */
			} else {
			    fprintf(stderr,"ERROR: harmonic potential for equivalence "
			            "group %i out of acceptable range.\n",i);
			    return 1;
		    }   }
		    if (f) {  /* harmonic, not a bond, and fitted ref value is in range */
			if (g) {  /* Fixed phase ==> not an angle ==> realign dihedral */
			    /* Faster fix for bond problem if there are no harmonic functions */
			    /* sum = x * ANGIMPF;  /* abuse sum */
			    delta = fabs(delta);
			    for (dipa=dipb;dipa;dipa=dipa->next) {
				curdih=dihmat+(dipa-dih)*ns;
				if (delta > 90.0)  /* closer to 180 than to 0 ==> shift to 180 */
				    for (dpa = solvec; dpa < dph; ) {
				        d = (fmod(*(curdih++) + 360.0, 360.0) - 180.0) * DEG_F;
				        /* Faster fix for bond problem if there are no harmonic functions */
				        /* *(dpa++) += d * d * sum; */
				        *(dpa++) += d * d * x;
				    }     /* next line: closer to 0 than to 180 ==> shift to 0 */
				  else for (dpa = solvec; dpa < dph; ) {
				        d = (fmod(*(curdih++) + 540.0, 360.0) - 180.0) * DEG_F;
				        /* Faster fix for bond problem if there are no harmonic functions */
				        /* *(dpa++) += d * d * sum; */
				        *(dpa++) += d * d * x;
			    }       }
			    eqga->dmin = HUGE_VAL;  /* trigger refit */
			    np2 += 2;  /* 2 more degrees of freedom that fall out */
			    eqga->delta1 = x;
			    eqga->delta2 = delta;
		    }   }
		} else dpf+=2;
	    }
	    u >>= 1;
	}
	i++;
    }
    normalize_weight (solvec);  /* needs to happen no matter what for calculating RMSE */
    if (np2) np2 = np -= np2;  /* 0 if everything is frozen. If so, no refit required. */
      else np2 = -1;  /* nothing frozen */
    if (np2 > 0) {
	printf("Freezing harmonic functions...\n");
	fflush(stdout);
	h=false;
	dpg = dpc = cosmat;
	dpd = amat;
	dpe = nbvec;
	nt = np2 + ns;  /* np2 depends on range check and is required for apllying null biases */
	k = i = 0;
        for (eqga = eqgrp; eqga < eqgmax; eqga++) {
	    g = (((u=eqga->mult) & PHASEFLAG) == 0);  /* Fixed phase requested */
	    b = (eqga->dmin == HUGE_VAL);  /* tagged for refit */
	    for (j=0;j<NMULT;j++) {
		if (u & 1) {
		    if (b && !j) {  /* tagged for refit and a harmonic */
		        dpe+=2;  /* was dpe++ , but a harmonic has 2 null biases */
		        h=true;
		        dpc += 2*ns;
		    } else for (f=false;;f=true) {
		            if (h) memcpy(dpg,dpc,scansz);
		            memcpy(dpd,dpc,scansz);
		            for (dpb = (dpa=(dpd+=ns)) + np; dpa < dpb; *dpa++ = 0.0);
		            dpd[k++]=*(dpe++);
		            dpd=dpb;
		            dpc+=ns;
		            dpg+=ns;
		            if (f || (j && g)) break;
		}       }
		u >>= 1;
	}   }
	printf("Solving linear system with fixed harmonic functions...\n");
	fflush(stdout);
	dpb=targvec;
	/* Faster fix for bond problem if there are harmonic functions. Was: */
	/* for (dpa = solvec; dpa < dph; *(dpb++) -= *(dpa++));  /* subtract from target! */
	/* Subtract solvec (potential due to frozen parmeter) from targvec.
	 * *= because we do use solvec later (eg. in dgemv_(...,&alpha,solvec,&one)) .
	 * As all frozen parameters in solvec are 0-fold dihedrals, it seems safe to assume
	 * they're impropers, hense ANGIMPF . */
	for (dpa = solvec; dpa < dph; *(dpb++) -= (*(dpa++) *= ANGIMPF));
	memcpy(parvec,targvec,scansz);  /* Next line: target for null bias should be zero. */
	for (dpb = (dpa=parvec+ns) + np; dpa < dpb; *dpa++ = 0.0);
	dgels_(&notrans,&nt,&np,&one,amat,&nt,parvec,&nt,work,&lwork,&info);
	if (info) {
	    fprintf(stderr,"ERROR: failed to solve linear system!\n");
	    return 1;
	}  /* Faster fix for bond problem if there are harmonic functions */
    } else for (dpa = solvec; dpa < dph; *(dpa++) *= ANGIMPF);
    free (nbvec);  /* This should be moved once we have per-multiplicity null biases. */
    /* TODO LATER: iterative adaptation of null bias to fit maximum k */
    /* When we implement iterative adaptation of null bias to fit maximum k, we'll have
     * to copy cosmat to amat and targvec to parvec repeatedly, which seems
     * inefficient. This might be a reason to revert to dspsv_ , but probably not,
     * because constructing MtM and MtV repeatedly doesn't seem that much better. */
    free(work);
    /* Apply bias compensation */
    if (biascomp) {
	dpc = biasvec;
	for (dpb = (dpa=parvec) + np; dpa < dpb; *dpa++ /= (1.0 - *dpc++));
    }
    free(biasvec);

    /****************\
    * Output results *
    \****************/
    /* Note that we cannot obtain an RMSD without null bias penalty other than by
     * calculating solvec */
    switch (np2) {
      case 0:   /* everything frozen */
	if (!biascomp) break;  /* if (biascomp) fall through and recalculate solvec */
      case -1:  /* nothing frozen ==> no refit ==> solution is simple matrix-vector product */
        dgemv_(&notrans,&ns,&np,&alpha,cosmat,&ns,parvec,&one,&beta,solvec,&one);
	break;  /* out of switch */
      default:  /* some degrees of freedom frozen ==> refit happened ==> add solvec to solution */
	dgemv_(&notrans,&ns,&np,&alpha,cosmat,&ns,parvec,&one,&alpha,solvec,&one);
    }
    dpc=solvec;  /* solvec is now the actual solution! */
    sum=0.0;
    if (outfile) {
	/* Align QM energy - small but important differences with normalize_weight */
	/* memo: ghi: # of groups */
	if (! (qmnorm = malloc(sizeof(struct normal) * ghi))) {
	    fprintf(stderr,"ERROR: can't allocate memory for normalizing QM energy!\n");
	    return 1;
	}
	if (! (dpa = dpd = malloc(sizeof(double) * ghi))) {
            fprintf(stderr,"ERROR: can't allocate memory "
			    "for aligning QM and MM energy!\n");
            return 1;
        }  /* grind... grind... grind... */
        for (maxnor = (curnor=qmnorm) + ghi; curnor < maxnor; ) {
	    curnor->avg=0.0;
	    curnor++->cnt=0;
	    *(dpa++)=HUGE_VAL;  /* greater than any other value in following comparisons */
	}
		  /* Find minimum QM energy for each goup and store in dpd (which we */
	dpa=qme;  /* just initialized to HUGE_VAL). Also make sums for group averages. */
	for (curgrpv=groupvec; curgrpv < maxgrpv; ) {
	    if ((d=*(dpa++)) < *(dpb=&dpd[i=*(curgrpv++)])) *dpb=d;
	    (curnor = &qmnorm[i])->avg += d;
	    curnor->cnt++;
	}
	dpa=dpd;  /* next line: calculate normalization offsets (re-using maxnor) */
	for (curnor=qmnorm; curnor < maxnor; curnor++) {
	    curnor->avg = curnor->avg / curnor->cnt + (*dpa=-*dpa);
	    dpa++;
	}
	dpa=qme;
	for (curgrpv=groupvec; curgrpv < maxgrpv; *(dpa++) += dpd[*(curgrpv++)]);
	free(dpd);
	/* Align MM energy - more similar to normalize_weight but still different */
	for (maxnor = (curnor=norm) + ghi; curnor < maxnor; curnor->avg=0.0, curnor++->cnt=0);
	dpa=mme;
	for (curgrpv=groupvec; curgrpv < maxgrpv; ) {
	    (curnor = &norm[*(curgrpv++)])->avg += *(dpa++);
	    curnor->cnt++;
	}  /* next line: calculate normalization offsets for groups (re-using maxnor) */
	curqmn=qmnorm;
	for (curnor=norm; curnor < maxnor; curnor++)
	    curnor->avg = curqmn++->avg - curnor->avg / curnor->cnt;
	free(qmnorm);
	dpa=mme;
        for (curgrpv=groupvec; curgrpv < maxgrpv; *(dpa++) += norm[*(curgrpv++)].avg);
	dpd = qme;
	dpe = mme;
	if (weighvec) dpf = weighvec;
	for (dpb = (dpa=targvec) + ns; dpa < dpb; ) {
	    d = *(dpa++) - (e=*(dpc++));  /* targvec - solvec */
	    sum += d * d;
	    if (weighvec) e /= *(dpf++);
	    x=*(dpe++);
	    fprintf (outfile,"%10.4f%10.4f%10.4f%10.4f\n",*(dpd++),x,x+e,d);
	}
	fclose(outfile);
    } else
	/* the following gives the RMSE *including* null bias penalty, which generally is
	 * not useful to the user. */
	/* for (dpb = (dpa=parvec+np) + ns; dpa < dpb; ) {
	 *     d = *(dpa++);
	 *     sum += d * d;
	 * } */
	for (dpb = (dpa=targvec) + ns; dpa < dpb; ) {
	    d = *(dpa++) - *(dpc++);
	    sum += d * d;
	}
    printf("Final RMSE = %g\n",sum = sqrt(sum/ns));
    fflush(stdout);  /* for good measure */
    if (! (outfile = fopen_clean(argv[1],"w","","parameter"))) return 1;
    fprintf(outfile,"* Parameters generated by analogy by\n* %s version %s\n"
		    "* RMSE = %g\n*\n\nBONDS\n",fullname,version,sum);
    h=false;
    k=ANGATT;
    for (l=1; l<5; l++) {  /* Woe, such lazy programming! */
	switch (l) {
	  case 2:
	    fprintf(outfile,"\nANGLES\n");
	    k=LASTATT;
	    break;  /* out of switch */
	  case 3:
	    fprintf(outfile,"\nDIHEDRALS\n");
	    k=DHASH-1;
	    break;  /* out of switch */
	  case 4:
	    fprintf(outfile,"\nIMPROPERS\n");
	    break;  /* out of switch for good measure */
	}
	i=1;
        dpa=parvec;
        for (eqga = eqgrp; eqga < eqgmax; eqga++) {
	    /* Need to compute this outside of the condition because K&R p.22 & 41:
	     * "Expressions connected by && or || are evaluated left to right, and
	     * evaluation stops as soon as the truth or falsehood is known." */
	    t2 = (dipa=eqga->first)->targ2;
	    /* can't have g (fixed phase; dihedrals only) and targ2 at the same time */
	    if ((g = (((v=eqga->mult) & PHASEFLAG) == 0)) && t2) {
		fprintf(stderr,"ERROR: consistency error #8; please contact a developer!\n");
		return 1;
	    }
	    f = false;  /* true if we want to skip */
	    if ((cpc=dipa->hash)[k]) f=true;  /* skip for current l */
	      else if (v & 1) {  /* there's a 0-fold */
		if (l == 3) f=true;  /* wait until l=4 (impropers) */
		  else {  /* don't skip for current l */
	            d1 = eqga->delta1;
	            d2 = eqga->delta2;
	    }   }
	    b = ((da=eqga->dmin) == HUGE_VAL);  /* tagged for refit */
	    /* TODO: it's lame (and asking for trouble) to go through all this just to ensure
	     * that dpa points to the right number. We should probably revert this part to
	     * CVS version 1.38, then equip struct eqgroup with a pointer or index. */
	    if (f) {  /* skip, but do update dpa */
		u=v;
		for (j=0;j<NMULT;j++) {
		    if (u & 1) {
			if (j) {  /* cosine (never gets eliminated from fit) */
			    if (g) dpa++;
			      else dpa+=2;
			} else if (!b) dpa+=2;  /* harmonic that was not eliminated */
		    }
		    u >>= 1;
		}
		continue;
	    }  /* Now f gets repurposed... */
	    f = ((db=eqga->dmax) == HUGE_VAL);  /* out of bounds warning */
	    m = dipa->tmult;   /* t2 has been assigned above but we want to keep m */
	    t1 = dipa->targ1;  /*and t1 here, after the "continue" a few lines up. */
	    cpa=&empty;  /* trigger non-match on first iteration */
	    for (;dipa;dipa=dipa->next) /* find unique dihedrals in group */
		if (strncmp(cpa,cpb=dipa->hash,DHASH)) {  /* not the same */
		    /* Within the current loop, we may output the exact same set of multiplicities
		     * multiple times if there are different parameters in the same eq group. */
		    dpb=dpa;
		    u=v;
		    for (j=0;j<NMULT;j++) {
			if (u & 1) {
			    comm=&empty;  /* extra warning message in prm file */
			    /* We're doing the calculation over and over. The alternative would be to
			     * first do the calcualtion for all the multiplicities, saving the results
			     * to a data structure, and then loop over all the dihedrals in the eq
			     * group, reading the amplitudes and phases back from the data structure.
			     * This would only be a bit more work, but it may just as well decrease
			     * efficiency given the rarity of manually defined equivalencies. At any
			     * rate, trying to shave a few cycles from an output routine is madness. */
			    if (j) {  /* cosine */
				if (g) {  /* fixed phase */
				    if ((x = j == m ? *(dpb++) + t1 : *(dpb++)) < 0) {
				        x=-x;
				        delta=180.0;
				    } else delta=0.0;
				} else {  /* variable phase */
				    d=*(dpb++);  /* amplitudes */
				    e=*(dpb++);  /* Phasor addition for theta1=0 and theta2=90: easy! */
				    x = sqrt(d*d + e*e);
				    /* we do need to use atan2 (rather than atan) to recover the sign */
				    delta = fmod(atan2(e,d)*F_DEG+t2+495.0,360.0) - 180.0;
				        /* Naive implementation of general phasor addition. TODO: */
				             /* contains mathematical redundancies; simplify (and */
				    if (j == m) {        /* eliminate y from main() if possible)! */
					d = x*cos(delta*=DEG_F) + t1*cos(y=t2*DEG_F);
					e = x*sin(delta) + t1*sin(y);
					x = sqrt(d*d + e*e);
					delta = fmod(atan2(e,d)*F_DEG+540.0,360.0) - 180.0;
				}   }
			    } else {  /* harmonic */
				if (b) {  /* Was realigned & frozen, so retrieve the values from before the freezing... */
				    x = d1;
				    delta = d2;  /* ...except that delta = fabs(delta) */
				    if (f) comm = oobwarn;
				} else {  /* Must necessarily have been variable phase */
				    d=*(dpb++);  /* d, e : force constants */
				    e=*(dpb++);  /* We cannot add targ[12] here because the check below is about the quality of the fit */
				    delta = (d*d1+e*d2)/(x=d+e);  /* d1, d2 not in radians! */
#ifdef DEBUG
				    fprintf(stderr,"DEBUG: k1 = %g ; k2 = %g ; k = %g ; delta = %g\n",d,e,x,delta);
#endif
				    /* don't need to check range again if we didn't refit */
				    if (np2 > 0 && (delta < da || delta > db)) {
					if (! m && t1) fprintf(stderr,"Warning: harmonic potential for equivalence "
				                "group %i out of acceptable range after refitting.\n",i);
				          else {
					    fprintf(stderr,"ERROR: harmonic potential for equivalence "
				                "group %i out of acceptable range after refitting.\n",i);
				            return 1;   /* Can be avoided with truly iterative refitting */
				}   }   }  /* but expected to be extremely rare in realistic situations. */
				/* This is a bit unelegant but will disappear once we have per-multiplicity weight factors;
				 * at that time, we might even want to perform the multiplication inside the if(b). */
				x *= cpb[ANGATT] ? ANGIMPF : BONDF;
				if (g) {  /* fixed phase harmonic dihedral (bonds and angles cannot be fixed) */
				    delta = fabs(delta);  /* same as in "verifying harmonic functions" section */
				    if (delta > 90.0) {   /* but for refit and/or with extra shiftwarn. */
				        if (delta < 180.0-SMALLSHIFT) {
					    h = true;
					    comm=shiftwarn;
					}
				        delta = 180.0;
				    } else {
				        if (delta > SMALLSHIFT) {
					    h = true;
					    comm=shiftwarn;
					}
				        delta = 0.0;
				    }
				    if (!m) {
					/* For harmonic potentials, t1 (or amp in readvect()) is not wrapped
					 * positive, otherwise the cusps would be switched with harmonic extrema. */
					if (delta != t2) {     /* This is also the reason for the non-additivity. */
					    /* fprintf(stderr,"ERROR: fitted harmonic with phase 180 for equivalence group "
					     *             "%i\nnot additive with phase %7.2f in %s !\n",i,t2,dipa->name); */
					    fprintf(stderr,"ERROR: harmonic with phase %7.2f in %s\nnot additive "
					        "with fitted phase 180 for equivalence group %i !\n",t2,dipa->name,i);
					    return 1;
					}
					x += t2;
				    }
				} else if (!m) {  /* variable phase harmonic (bond, angle or dihedral) */
				    /* The fitted phase is guaranteed to be wrapped between -180.0 and 180.0, but
				     * the same can't be said for the target phase from the measurement file.
				     * However, we can't be sure wrapping the latter (here or in readvect()) won't
				     * have a different effect than the creator of the force field intended; since
				     * harmonics with extreme phases are mostly unexplored territory, there might
				     * even exist differences between programs. Therefore, we leave this responsi-
				     * bility to the user; the below warning and error should be clear enough. */
				    /* Legitimate use cases may occur for bonds, angles and well-behaved impropers
				     * in strained structures when fitting conformations close to the equilibrium
				     * but far from the reference value, so we can't make this condition fatal. */
				    if (t2 < d1 || t2 > d2) fprintf(stderr,"Warning: harmonic with reference "
					    "value %7.2f in %s\noutside of measurement data range "
					    "[%7.2f:%7.2f] for equivalence group %i\n",t2,dipa->name,d1,d2,i);
				    /* May not be safe with regard to x; we don't take the chance. */
				    /* delta = (x*delta+t1*t2)/(x+=t1) */
				    delta = (x*delta+t1*t2)/(e=x+t1);  /* recycle e instead */
				    x = e;
				    if (delta < da || delta > db) {
					/* This is getting so far-fetched we won't even bother
					 * giving more details in the error message. */
				        fprintf(stderr,"ERROR: final harmonic potential for "
				                "equivalence group %i out of acceptable range.\n",i);
				        return 1;
			    }   }   }
			    if (! cpb[ANGATT])
			        fprintf (outfile,"%s %7.2f %10.4f ! RMSE = %g\n",cpb,x,delta,sum);
			      else if (! cpb[LASTATT])
			        fprintf (outfile,"%s  %6.2f    %6.2f ! RMSE = %g\n",cpb,x,delta,sum);
			      else fprintf (outfile,"%s   %8.4f  %1i  %7.2f ! RMSE = %g%s\n",
						     cpb,x,j,delta,sum,comm);
			}
			u >>= 1;
		    }
		    cpa=cpb;
		}
	    dpa=dpb;  /* finally we can move on to the next eq group */
	    cpc[DHASH-1] = cpc[LASTATT] = '\1';  /* dangerous vandalism */
	    i++;
    }   }
    fprintf(outfile,"\nEND\n");
    if (h) fprintf(stderr,"warning: some 0-folds have a phase shifts > %g degrees; "
						"fit may be inaccurate\n",SMALLSHIFT);
    printf("Done.\n");

#ifdef DEBUG
    /* Freeing all the memory before exiting is redundant, but we do it anyway to shut up
     * valgrind, which got confused by all the pointer manipulations we're doing. */
    /* Iteratively free memory for dihedral file names. */
    for (cppmax = (cppcur=fn) + fns; cppcur < cppmax; free(*(cppcur++)));
    fclose(outfile);  /* Opening files allocates memory and closing them frees it */
    /* TODO: this probably needs to be reorganized */
    free(amat);
    free(cosmat);
    free(parvec);
    free(eqgrp);
    free(dih);
    free(dihmat);
    free(norm);
    free(groupvec);
    free(weighvec); /* K&R APPENDIX B: "[free] does nothing if [its argument] is NULL" */
    free(mme);      /* " " */
    if (targvec != qme) free(targvec);
    free(qme);
#endif  /* DEBUG */
    return 0;
}

void normalize_weight (double *tvec) {
    int *curgrpv;
    double *dpa,*dpb;
    struct normal *curnor,*maxnor;

    for (maxnor = (curnor=norm) + ghi; curnor < maxnor; curnor->avg=0.0, curnor++->cnt=0);
    /* sum energy differences */
    dpa=tvec;
    for (curgrpv=groupvec; curgrpv < maxgrpv; ) {
	(curnor = &norm[*(curgrpv++)])->avg += *(dpa++);
	curnor->cnt++;
    }  /* next line: calculate normalization offsets for groups (re-using maxnor) */
    for (curnor=norm; curnor < maxnor; curnor++) curnor->avg = - curnor->avg / curnor->cnt;
    dpa=tvec;
    if (weighvec) {
	dpb=weighvec;
        for (curgrpv=groupvec; curgrpv < maxgrpv; dpa++)  /* sic */
	    *dpa = (*dpa + norm[*(curgrpv++)].avg) * *(dpb++);
    } else for (curgrpv=groupvec; curgrpv < maxgrpv; *(dpa++) += norm[*(curgrpv++)].avg);
}

/* TODO LATER: right now, this is a stub, but if we allow for multiple multiplicity
 * sections per parameter, we may have to merge the bitfields somewhere (we can't simply
 * "or" them because we have to check for duplicated bits). Then again, struct eqgroup
 * may look very different at that point, and it's far from sure that we'll actually have
 * to merge them or that we'll do it inside this function. */
void finalize_mult (struct eqgroup *eqg, int *mult, int nmult, enum phaseflag phase) {
    eqg->mult = packmult(mult,nmult,phase==PH_VAR);
    eqg->ph_undef = phase==PH_UNDEF;
}

char packmult (int *mult, int nmult, bool phasevar) {
    int *himult;
    int k,l;
    char u=1;

    k=*(himult=&mult[nmult-1]);  /* we know that get_list sorted multemp */
    for (;himult>=mult;) {
        u <<= (k - (l=*(himult--)));
        u |= 1;
        k = l;
    }
    if (k) u <<= k;
    if (phasevar) {
        u |= PHASEFLAG;
        np += nmult*2;
    } else if (k) np += nmult;  /* k := ((k := mult[0]) != 0) */
      else np += nmult+1;
    return u;
}

const char *outsideopt (const char exact, const char fuzzy, const char *disp) {
    const opt_t *opts;
    const char *exitstatus=NULL;
    char key;

    gopt_arg(options, exact, &exitstatus);
    if (! gopt(options,fuzzy)) return exitstatus; /* NULL is a valid return value */
    for (opts=options;;opts++) {
	switch (key=opts->key) {
	  case 'p': case '\0': return exitstatus; /* loop's only non-error exit point */
	}
	if (key == fuzzy) {
	    if (exitstatus) {
	        fprintf(stderr,"ERROR: %s specified two or more times!\n",disp);
	        return &empty;  /* sentinel */
	    }
	    exitstatus=opts->arg;
}   }   }  /* repeat loop because we can have multiple fuzzies */

void listeqg (FILE *unit, char *pre, struct eqgroup *eqg, int indent) {
    struct dihedral *dip;
    dip=eqg->first;  /* caller knows eqg->first and (eqg-eqgrp+1) but this is I/O */
    if (eqg->merged) fprintf(unit,"%squivalence group %i :\n",pre,(int) (eqg-eqgrp+1));
      else fprintf(unit,"%squivalence group %i (%s):\n",pre,(int) (eqg-eqgrp+1),dip->hash);
    for (;dip;dip=dip->next)
	fprintf(unit,"%*s| measurement file #%i = %s\n",indent,&empty,(int) (dip-dih+1),dip->name);
}

/* warning: promptenter does not require disp to have a trailing space, while
 * all the other functions do! */
bool promptenter (char *buf, const char *filename, enum datatyp kind,
		void **vect, fint *dim, char *disp, char *action) {
    void *vec;
    size_t size;

    if (! filename) {
	if (interact) {
	    printf("%sfile (ENTER = no %s): ",disp,action);
	    fflush(stdout);
	    if (! fgets(buf,MAXLINE,stdin)) {
	        fprintf(stderr,"ERROR: end of input stream "
			    "while reading %sfile name!\n",disp);
	        return true;
	    }
	filename=buf;
	} else return false;  /* *vect left at input value */
    }
    if (*filename == '\n') return false;
    if (! (vec=*vect)) {
        if (! (size=datatypsize(kind))) return true;
        if (! (vec = malloc(size * *dim))) {
	    fprintf(stderr,"ERROR: can't allocate memory for %s!\n",action);
	    return true;
    }   }
    if (readvect(filename,kind,vec,dim,NULL,disp)) return true;
    *vect=vec;
    return false;
}

/* caller may not need filename but then one may think of it as "passing a buffer" */
bool promptread (char *buf, const char *filename, enum datatyp kind,
		void *vec, fint *dim, struct dihedral *dihe, char *disp) {
    if (filename) return readvect(filename,kind,vec,dim,dihe,disp);
    if (! interact) {
	fprintf(stderr,"ERROR: %sfile not specified "
			"in non-interactive mode!\n",disp);
	return true;
    }
    for (;;) {
        printf("%sfile: ",disp);
	fflush(stdout);
        if (! fgets(buf,MAXLINE,stdin)) {
	    fprintf(stderr,"ERROR: end of input stream "
			    "while reading %sfile name!\n",disp);
	    return true;
	}
	if (*buf != '\n') return readvect(buf,kind,vec,dim,dihe,disp);
	fprintf(stderr,"warning: no name given for %sfile\n",disp);
}   }

/* The part with "kind" is not the most efficient thing to do,
 * but we're just reading files... */
bool readvect (const char *filename, enum datatyp kind,
		void *vec, fint *dim, struct dihedral *dihe, char *disp) {
    void *cur,*max;
    char *cp;
    FILE *file;
    double d;
    fint n,m;
    size_t size;
    float f;
    int i,j,k;

    if (cp=strchr(filename,'\n')) *cp='\0';
    if (! (file = fopen_clean(filename, "r", "readvect", disp))) return true;
    /* See you in ~100 lines, folks! */
    if (dihe) {  /* next 3 lines: allocate on stack only if (dihe) */
	char *hash = dihe->hash;
	char line[MAXLINE],att[4][ATTL],c;
	double fc,amp,phase;  /* float in "cosf" branch */
	unsigned int mult;
	if (! fgets(line,MAXLINE,file)) {
	    fprintf(stderr,"readvect ERROR: %sfile %s is empty!\n",disp,filename);
	    return true;
	}  /* Next line: %lf becomes %f in "cosf" branch */
	if ((i=sscanf(line,"%" ATTF " %" ATTF " %" ATTF " %" ATTF " %lf %u %lf"
		,att[0],att[1],att[2],att[3],&amp,&mult,&phase)) > 4) k=4;
	  else k=i;  /* Next line: use sscanf as a lazy man's tool to determine whether
	       * it's a number (and if it is, keep it for later). %f in "cosf" branch. */
	for (j=0;j<k;j++) if (sscanf(att[j],"%lf",&fc) == 1) break;
	/* Does not *need* to be initialized, but it saves a few lines of code,
	 * is less error-prone, and this is an I/O routine... */
	dihe->targ2 = dihe->targ1 = 0.0;
	switch (j) {  /* Position (field) of first number - first one is 0 */
	  case 2:
	    switch (i) {
	      case 2:
		dihe->tmult = -1;
		break;  /* out of inner switch (i) */
	      case 3:
		fprintf(stderr,"readvect ERROR: incomplete bond parameter "
			"in %sfile %s !\n",disp,filename);
		return true;
	      default:  /* i necessarily >3 */
		/* Since we don't know yet which reference angles the program will pick,
		 * we we simply store the force constant and reference value. This turns
		 * out to be easier for final processing anyway. Same for other harmonics
		 * incl. impropers. */
		if (sscanf(att[3],"%lf",&dihe->targ2) != 1) {
		    fprintf(stderr,"readvect ERROR: 4th field of bond parameter "
			"is not a number in %sfile %s !\n",disp,filename);
		    return true;
		}		  /* It is tempting to divide by BONDF here, but that would */
		dihe->targ1 = fc; /* frustrate implementing per-parameter BONDF later. Same */
		dihe->tmult = 0;   /* thing for other harmonics further down this function. */
		break;  /* Out of switch for good measure */
	    }
	    break;  /* out of switch */
	  case 3:
	    switch (i) {
	      case 3:
		dihe->tmult = -1;
		break;  /* out of inner switch (i) */
	      case 4:
		fprintf(stderr,"readvect ERROR: incomplete bond parameter "
			"in %sfile %s !\n",disp,filename);
		return true;
	      case 7:  /* In all other cases, we silently ignore extraneous fields, but
			this one is a particularly nasty trap for the uninformed user. */
		fprintf(stderr,"readvect warning: suspected Urey-Bradley in %sfile %s "
			"ignored; please define it as a separate bond!\n",disp,filename);
	      default:  /* ACHTUNG! Fell through! */
		dihe->tmult = 0;
		dihe->targ1 = fc;
		dihe->targ2 = amp;
		break;  /* Out of switch for good measure */
	    }
	    break;  /* out of switch */
	  case 4:
	    switch (i) {
	      case 4:
		dihe->tmult = -1;
		break;  /* out of inner switch (i) */
	      case 5: case 6:
		fprintf(stderr,"readvect ERROR: incomplete dihedral parameter "
			"in %sfile %s !\n",disp,filename);
		return true;
	      default:  /* can only be 7 but default is lighter and more extensible */
		if (mult < 0 || mult > 6) {
		    fprintf(stderr,"readvect ERROR: invalid multiplicity %u "
		            "in %sfile %s !\n",mult,disp,filename);
		    return true;
		}
		dihe->tmult = mult;  /* next line: targ2 already initialized to 0; */
				   /* 180 would be identified as "odd target phase" */
		if (phase == 180.0) dihe->targ1 = -amp;
		  else {
		    dihe->targ1 = amp;
		    dihe->targ2 = phase;
		}
		break;  /* Out of switch for good measure */
	    }
	    break;  /* out of switch */
	  default:
	    fprintf(stderr,"readvect ERROR: invalid first line "
			    "in %sfile %s !\n",disp,filename);
	    return true;
	}
	if (parhash(hash,att,j)) {
	    fprintf(stderr,"ERROR: consistency error #4; please contact a developer!\n");
	    return true;  /* can't happen because just caught in above "default:" */
	}
	if (j == 2) hash[LASTATT]='\0';  /* For output routine; the sprintf() call in */
    }                                    /* parhash() already put a \0 at ANGATT-1 . */
    cur=vec;
    if (! (size=datatypsize(kind))) return true;
    /* this usage of size is permitted per http://c-faq.com/ptrs/castincr.html */
    if (n=*dim) max = (char *) vec + size*n;
      else max = (char *) vec + size*OPT_MAXT;
    /* Reminder: scanf works in terms of whitespace-separated records, where
     * newlines are treated as whitespace. */
    for (;;) {
	switch (kind) {
	  case INT:
	    j=fscanf(file,"%i",&i);
	    break;  /* switch */
	  case FLOAT:
	    j=fscanf(file,"%f",&f);
	    break;  /* switch */
	  case DOUBL:
	    j=fscanf(file,"%lf",&d);
	    break;  /* switch */
	}
	if (j==EOF) break;
        if (j) {
	    if (cur==max) {
                fprintf(stderr,"readvect ERROR: more than %li records "
			    "in %sfile %s !\n",(long int) n,disp,filename);
		return true;
	    }
	    switch (kind) {
	      case INT:
	        *((int *) cur) = i;
		break;
	      case FLOAT:
	        *((float *) cur) = f;
		break;
	      case DOUBL:
	        *((double *) cur) = d;
		break;
	    }
	    cur = (char *) cur + size;
	} else fprintf(stderr,"readvect warning: invalid record in %sfile %s "
			"ignored\n",disp,filename);
    }
    fclose(file);
    m = ((char *) cur - (char *) vec) / size;
    if (n) {
        if (n != m) {
            fprintf(stderr,"readvect ERROR: %li records read (%li expected) "
			    "in %sfile %s !\n",(long int) m,(long int) n,disp,filename);
	    return true;
	}
    } else *dim=m;
    return false;
}

struct eqgroup *findgrp(const char *buf, int swit) {
    struct dihedral *dipa;
    struct eqgroup *eqga;
    char atty[4][ATTL];
    char query[DHASH];

    /* Yes, we checked this is safe */
    if (parhash(query,atty,sscanf(buf,"%" ATTF " %" ATTF " %" ATTF " %" ATTF " %" ATTF,
					atty[0],atty[1],atty[2],atty[3],query))) {
	fprintf(stderr,"ERROR: wrong # atom types in argument -%c '%s'!\n",swit,buf);
	return NULL;
    }
    for (eqga = eqgrp; eqga < eqgmax; eqga++)
	for (dipa=eqga->first;dipa;dipa=dipa->next)     /* user might have equivalenced */
	     if (! strncmp(query,dipa->hash,DHASH)) return eqga; /* different dihedrals */
    fprintf(stderr,"ERROR: parameter not found in argument -%c '%s'!\n",swit,buf);
    return NULL;
}

/* Not much of a hash function - just performs alphabetic horizontal sorting
 * and stores result as single human-readable string (like in dihmunge).
 * Which is all that's needed for the present program. */
bool parhash(char *hash, char att[][ATTL], int nat) {  /* K&R p112 */
    int i;

    switch (nat) {
      case 2:
	if (strncmp(att[0],att[1],HASHSZ) > 0)
	    sprintf(hash,"%-" ATTF " %-" ATTF,att[1],att[0]);
	  else sprintf(hash,"%-" ATTF " %-" ATTF,att[0],att[1]);
	return false;  /* break out of switch and function at the same time */
      case 3:
	if (strncmp(att[0],att[2],HASHSZ) > 0)
	    sprintf(hash,"%-" ATTF " %-" ATTF " %-" ATTF,att[2],att[1],att[0]);
	  else sprintf(hash,"%-" ATTF " %-" ATTF " %-" ATTF,att[0],att[1],att[2]);
	return false;  /* break out of switch and function at the same time */
      case 4:
	if (! (i = strncmp(att[1],att[2],HASHSZ))) i = strncmp(att[0],att[3],HASHSZ);
	if (i > 0) sprintf(hash,"%-" ATTF " %-" ATTF " %-" ATTF " %-" ATTF,
			    att[3],att[2],att[1],att[0]);
	  else sprintf(hash,"%-" ATTF " %-" ATTF " %-" ATTF " %-" ATTF,
			    att[0],att[1],att[2],att[3]);
	return false;  /* break out of switch and function at the same time */
      default:
	return true;
}   }

/* caller may not need line but then one may think of it as "passing a buffer" */
int get_list (char *line, enum datatyp kind, void *array, int minsz, int maxsz,
		int minint, int maxint, char *disp_s, char *disp_p) {
    char *cp;
    void *vpa,*vpb,*vpc;
    int *ipa,*ipb,*ipc;
    double d;
    size_t size;
    float f;
    int i,j,k;

    if (! fgets(cp=line,MAXLINE,stdin)) {
        fprintf(stderr,"ERROR: end of input stream while reading %s!\n",disp_p);
        return -3;
    }
    if (line[0] == '\n') return 0;
    if (! (size=datatypsize(kind))) return -3;
    /* At fist glance, it seems sensible to make the following routine shorter by
     * exploiting (1) sscanf's feature to ignore leading whitespace and (2) sscanf's
     * %n directive. In practice, both shortcuts may cause trouble. Leading whitespace
     * may contain the terminal newline, which we'd need to distinguish from an invalid
     * input (sscanf will return 0 in both cases). Invalid inputs will not be matched by
     * %i and will not be inclused in %n's count, causing an endless loop. See also:
     * http://c-faq.com/stdio/scanfjam.html
     * http://c-faq.com/stdio/scanfprobs.html
     * The code below, while long and ugly, is the only 100% foolproof solution I could
     * come up with. */
    for (vpc = (char *) (vpa=array) + size*maxsz;;) {  /* repeat to read in every element of line */
	for (j=0;!j;) switch (*cp) {  /* skip whitespace */
	      case  ' ': case '\t': cp++; break;  /* out of switch */
	      case '\n': case '\0': j=2;  break;
	      default:              j=1;  break;
	    }
	if (j==2) break;
	switch (kind) {
	  case INT:
	    k=sscanf(cp,"%i",&i);
	    break;  /* switch */
	  case FLOAT:
	    k=sscanf(cp,"%f",&f);
	    break;  /* switch */
	  case DOUBL:
	    k=sscanf(cp,"%lf",&d);
	    break;  /* switch */
	}
	if (k != 1) {
	    fprintf(stderr,"warning: incorrectly formatted input field ignored\n");
	    j=3;
	} else if ((! kind) && (i < minint || i > maxint)) {
	    fprintf(stderr,"warning: %s %i out of range %i-%i ignored\n",
			    disp_s,i,minint,maxint);
	    j=3;
	}
	if (j == 1) {
	    if (vpa == vpc) {  /* it should not be possible to go out of bounds */
	        fprintf(stderr,"warning: more than %i %s specified\n",maxsz,disp_p);
	        return -2;
	    }
	    switch (kind) {
	      case INT:
	        *((int *) vpa) = i;
		break;
	      case FLOAT:
	        *((float *) vpa) = f;
		break;
	      case DOUBL:
	        *((double *) vpa) = d;
		break;
	    }
	    vpa = (char *) vpa + size;
	}
	for (j=0;!j;) switch (*cp) {  /* skip non-whitespace string */
	      case  ' ': case '\t': j=1;  break;
	      case '\n': case '\0': j=2;  break;
	      default:              cp++; break;  /* out of switch */
	    }
	if (j==2) break;
    }
    if ((j = ((char *) vpa - (char *) array) / size) < minsz) {
        fprintf(stderr,"warning: less than %i valid %s specified\n",minsz,disp_p);
        return -1;
    }
    if (! kind) {  /* ie. (kind == INT) */
	qsort(ipb=array,j,sizeof(int),intcmp);
	/* check for repeats */
	k=*(ipb++);
	ipc=NULL;
	for(ipa = (int *) array + j; ipb < ipa; ipb++) {
	    if ((i=*ipb) == k) {
	        ipc=ipb++;
	        break;
	    }
	    k=i;
	}
	if (ipc) {
	    fprintf(stderr,"warning: repeated %s ignored\n",disp_p);
	    for(;ipb<ipa;) if ((i=*(ipb++)) != k) *(ipc++) = k = i;
	    j = ipc - (int *) array;
	}
    }
    if (j < minsz) {
        fprintf(stderr,"warning: less than %i unique %s specified\n",minsz,disp_p);
        return -1;
    }
    return j;
}

size_t datatypsize (enum datatyp kind) {
    switch (kind) {
      case INT:
	return sizeof(int);
      case FLOAT:
	return sizeof(float);
      case DOUBL:
	return sizeof(double);
      default:
	fprintf(stderr,"ERROR: consistency error #1; please contact a developer!\n");
	return 0;
}   }

enum ouinon yesno (char *line, char *disp, enum ouinon def) {
    char answer[4];  /* English: "yes\0" */
    char *cp;

    for (;;) {
	if (! fgets(cp=line,MAXLINE,stdin)) {  /* "checking", "verifying" */
	    fprintf(stderr,"ERROR: end of input stream while confirming %s!\n",disp);
	    return YNVOID;
	}
	switch(sscanf(line,"%3s",answer)) {
	  case 0:
	    if (def != YNVOID) return def;
	    break;  /* out of switch */
	  case 1:
	    strnupper (answer,3);
	    /* K&R p.21 & 41: Expressions connected by && or || are evaluated left to
	     * right, and evaluation stops as soon as the truth or falsehood is known. */
	    if (*answer == 'N' && (*(cp=answer+1) == '\0' ||
		    (*cp == 'O' && *++cp == '\0'))) return NO;
	    if (*answer == 'Y' && (*(cp=answer+1) == '\0' ||
		    (*cp == 'E' && *++cp == 'S' && *++cp == '\0'))) return YES;
	    break;  /* out of switch */
	}
	printf("Unrecognised answer; please try again! ");
	fflush(stdout);
}   }  /* Can't happen */

bool strnupper (char *string, int n) {
    char *cp;
    char c;
    for (cp = string + n; string < cp; string++) switch (c = *string) {
	  case '\n': case '\0': return;
	  default: *string = toupper(c);  /* Watch out, no breaks */
}   }  /* no meaningful return value */

int intcmp (const void *a, const void *b) {
    return *((const int *)a) - *((const int *)b);
}  /* prev line: int subtraction faster than comparison on most modern architectures */

int floatcmp (const void *a, const void *b) {
    const double ia = *((const double *)a);  /* float in "cosf" branch */
    const double ib = *((const double *)b);  /* float in "cosf" branch */
    /* If the compiler is any smart at all, directly comparing floats should go
     * a lot faster than doing a floating point subtraction; see
     * http://www.cygnus-software.com/papers/comparingfloats/comparingfloats.htm
     * Also keep in mind that for unsorted data, in half of the cases, we get
     * away with only 1 comparison. */
    if (ia < ib) return -1;
    if (ia > ib) return 1;
    return 0;
}

FILE *fopen_clean(const char *filename, char *mode,
			char *function, char *disposition) {
    FILE *exitstatus;

    if (*filename) {
	if (exitstatus = fopen(filename,mode)) return exitstatus;
	fprintf(stderr,"%s ERROR: can't open %sfile \"%s\" !\n",
		    function,disposition,filename);
    } else fprintf(stderr,"%s ERROR: no name given for %sfile!\n",function,disposition);
    return NULL;
}

int numl(const char *s) {
    int i=0;
    for (;;) {
	switch (*s++) {
	  case '\0': break;  /* out of switch (and therefore out of for) */
	  case '\n': i++;  /* ACHTUNG: fall through to default: */
	  default: continue;  /* only way to repeat the loop */
	}
	break;
    }
    return i;
}
