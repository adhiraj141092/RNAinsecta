/* miranda.h									  */
/*--------------------------------------------------------------------------------*/
/* miRanda- An miRNA target scanner, aims to predict mRNA targets for microRNAs,  */
/* using dynamic-programming alignment and thermodynamics                         */
/*                                                                                */
/* Copyright (C) (2003) Memorial Sloan-Kettering Cancer Center, New York          */
/*                                                                                */
/* Distributed under the GNU Public License (GPL)                                 */
/* See the files 'COPYING' and 'LICENSE' for details                              */
/*                                                                                */
/* Authors: Anton Enright, Bino John, Chris Sander and Debora Marks               */
/* Email: mirnatargets@cbio.mskcc.org - reaches all authors                       */
/*                                                                                */
/* Written By: Anton Enright (enrighta@mskcc.org)                                 */
/*                                                                                */
/* Please send bug reports to: miranda@cbio.mskcc.org                             */
/*                                                                                */
/* If you use miRanda in your research please cite:                               */
/* Enright AJ, John B, Gaul U, Tuschl T, Sander C and Marks DS;                   */
/* (2003) Genome Biology; 5(1):R1.                                                */
/*                                                                                */
/* This software will be further developed under the open source model,           */
/* coordinated by Anton Enright and Chris Sander:                                 */
/* miranda@cbio.mskcc.org (reaches both).                                         */
/*--------------------------------------------------------------------------------*/
/*
 * Copyright (C) (2003) Memorial Sloan-Kettering Cancer Center
 *
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation; either
 * version 2 of the License, or (at your option) any later
 * version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
 */


/* Setup Traceback Flags */
#define DIAG 0
#define UP 1
#define LEFT -1
int CURR;

/* Setup Strands */
#define FORWARD 0
#define REVERSE 1

/* Adjustable Algorithm Parameters */
double scale;
int nomodel;
double gap_open;
double gap_extend;
double average;
double stdev;
double z_threshold;
double score_threshold;
double energy_threshold;
int do_shuffle;
int no_energy;
int shuffle_window;
int total_shuffles;
int verbosity;
unsigned int uniform;
int outfile;
int truncated;
int total_hits;

/* Output Filehandle */
FILE *fpout;

/* Generic structure to store individual hit information */
typedef struct hit_struct
{

  double score;
  int query_start;
  int query_end;
  int ref_start;
  int ref_end;
  char *alignment[3];
  char *rest[6];

}
hit_struct;

/* Generic structure to store consecutive UTR hit information */
typedef struct final_score
{
  double scan_score;
  int no_hits;
  double max_hit;
  double max_score;
  double total_score;
  char positional[1000];
}
final_score;

/* Score structure for non-optimal path detection */
typedef struct score_struct
{

  double score;
  int path;
  int i;
  int j;

}
score_struct;

/* Declare Functions */
double max (double, double, double);
double score (char, char);
void traceback (double **, int **, char *, char *, int, int, int,
		hit_struct *, int);
long readinseq (long, FILE *, char *, char *, char *);
void shuffle (char *, int, int);
void build_matrix (double **, int **, int **, char *, char *, int, int,
		   score_struct *);
double build_matrix_quick (double **, int **, char *, char *, int, int);

void revstring (char s[]);
void irand (int);
int nrand (int);
void printhit (char *, char *, hit_struct *, char *, char *, int, double, double, FILE *);
void print_license();
void print_small_license();
double get_energy (hit_struct *);
double vfold (char *);
void clear_hit (hit_struct *, int, int);
void clear_matrix (double **, int, int, int, int);
void initialize_bases ();
int cmpscores (const void *, const void *);
int find_targets (FILE *, FILE *, FILE *, char *);
double do_alignment (double **, int **, int **, char *, char *,
		     score_struct *, hit_struct *, int, int, int,
		     final_score *, int, char *, char *, FILE *);
void print_banner();
void print_usage();
void print_options();
