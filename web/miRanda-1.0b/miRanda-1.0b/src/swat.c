/* swat.c */
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


#include <config.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <ctype.h>
#include <string.h>
#include "utils.h"
#include "miranda.h"
#include "scmatrix.h"

void
initialize_bases ()
{
  bases['C'] = 0;
  bases['G'] = 1;
  bases['A'] = 2;
  bases['T'] = 3;
  bases['U'] = 4;
  bases['X'] = 5;
  bases['N'] = 5;
  bases['c'] = 0;
  bases['g'] = 1;
  bases['a'] = 2;
  bases['t'] = 3;
  bases['u'] = 4;
  bases['x'] = 5;
  bases['n'] = 5;

}

double
build_matrix_quick (double **m1, int **m2, char *sequence1,
		    char *sequence2, int seqlen1, int seqlen2)
{

  double penalty1 = 0;
  double penalty2 = 0;
  int i, j = 0;
  double maxi = 0;

  for (i = 1; i <= seqlen1; i++)
    {
      for (j = 1; j <= seqlen2; j++)
	{

	  if (m2[i][j - 1] == LEFT)
	    {
	      penalty1 = gap_extend;
	  } else
	    {
	      penalty1 = gap_open;
	    }

	  if (m2[i - 1][j] == UP)
	    {
	      penalty2 = gap_extend;
	  } else
	    {
	      penalty2 = gap_open;
	    }

	  m1[i][j] =
	    max (m1[i - 1][j - 1] +
		 (score (sequence1[i - 1], sequence2[j - 1])),
		 m1[i][j - 1] + (penalty1), m1[i - 1][j] + (penalty2));
	  m2[i][j] = CURR;

	  if (m1[i][j] < 0)
	    {
	      m1[i][j] = 0;
	    }

	  if (m1[i][j] > maxi)
	    {
	      maxi = m1[i][j];
	    }

	}
    }
  return (maxi);
}


void
build_matrix (double **m1, int **m2, int **m3, char *sequence1,
	      char *sequence2, int seqlen1, int seqlen2,
	      score_struct * scores)
{

  double penalty1 = 0;
  double penalty2 = 0;
  int i, j = 0;
  int path = 0;

  for (i = 1; i <= seqlen1; i++)
    {
      for (j = 1; j <= seqlen2; j++)
	{

	  if (m2[i][j - 1] == LEFT)
	    {
	      penalty1 = gap_extend;
	  } else
	    {
	      penalty1 = gap_open;
	    }

	  if (m2[i - 1][j] == UP)
	    {
	      penalty2 = gap_extend;
	  } else
	    {
	      penalty2 = gap_open;
	    }

	  if (i < (strlen (sequence1) - 10))
	    {
	      m1[i][j] =
		max (m1[i - 1][j - 1] +
		     (score (sequence1[i - 1], sequence2[j - 1])),
		     m1[i][j - 1] + (penalty1), m1[i - 1][j] + (penalty2));
	      m2[i][j] = CURR;
	      m3[i][j] = 0;
	  } else
	    {
	      m1[i][j] =
		max (m1[i - 1][j - 1] +
		     (scale * score (sequence1[i - 1], sequence2[j - 1])),
		     m1[i][j - 1] + (scale * penalty1),
		     m1[i - 1][j] + (scale * penalty2));
	      m2[i][j] = CURR;
	      m3[i][j] = 0;
	    }

	  if (m1[i][j] < 0)
	    {
	      m1[i][j] = 0;
	  } else
	    {
	      if (m2[i][j] == DIAG)
		{
		  if (m3[i - 1][j - 1] <= 0)
		    {
		      total_hits++;
		      m3[i][j] = path = total_hits;
		  } else
		    {
		      m3[i][j] = path = m3[i - 1][j - 1];
		    }
	      } else if (m2[i][j] == LEFT)
		{
		  if (m3[i][j - 1] <= 0)
		    {
		      total_hits++;
		      m3[i][j] = path = total_hits;
		  } else
		    {
		      m3[i][j] = path = m3[i][j - 1];

		    }
	      } else
		{
		  if (m3[i - 1][j] <= 0)
		    {
		      total_hits++;
		      m3[i][j] = path = total_hits;
		  } else
		    {
		      m3[i][j] = path = m3[i - 1][j];
		    }
		}

	      if (m1[i][j] >= scores[path].score)
		{
		  scores[path].score = m1[i][j];
		  scores[path].path = path;
		  scores[path].i = i;
		  scores[path].j = j;
		}
	    }
	}
    }
}

void
traceback (double **m1, int **m2, char *sequence1, char *sequence2, int i,
	   int j, int remove, hit_struct * hit_ptr, int length)
{

  double score = 0;

  if (length == 0)
    {
      hit_ptr->query_end = i;
      hit_ptr->ref_end = j;
    }
  if (m1[i][j] > 0)
    /*if (i>0 && j >0) */
    {

      if (m2[i][j] == DIAG)
	{
	  length++;
	  traceback (m1, m2, sequence1, sequence2, i - 1, j - 1, 1, hit_ptr,
		     length);

	  score =
	    (double)
	    match[bases[(int) sequence1[i - 1]]][bases
						 [(int) (sequence2[j - 1])]];
	  if (i >= (strlen (sequence1) - 10))
	    {
	      score = score * scale;
	    }

	  hit_ptr->score += score;
	  hit_ptr->alignment[0][length - 1] = sequence1[i - 1];
	  hit_ptr->alignment[2][length - 1] = sequence2[j - 1];

	  if (score > 2)
	    {
	      hit_ptr->alignment[1][length - 1] = '|';
	  } else if ((score >= 0) && (score <= 4))
	    {
	      hit_ptr->alignment[1][length - 1] = ':';
	  } else
	    {
	      hit_ptr->alignment[1][length - 1] = ' ';
	    }
	  if (remove)
	    {
	      m1[i][j] = 0;
	    }
      } else if (m2[i][j] == UP)
	{

	  if (m2[i - 1][j] == UP)
	    {
	      if (i < (strlen (sequence1) - 10))
		{
		  hit_ptr->score += gap_extend;
	      } else
		{
		  hit_ptr->score += (gap_extend * scale);
		}
	  } else
	    {
	      if (i < (strlen (sequence1) - 10))
		{
		  hit_ptr->score += gap_open;
	      } else
		{
		  hit_ptr->score += (gap_open * scale);
		}
	    }


	  length++;
	  traceback (m1, m2, sequence1, sequence2, i - 1, j, 1, hit_ptr,
		     length);
	  hit_ptr->alignment[0][length - 1] = sequence1[i - 1];
	  hit_ptr->alignment[1][length - 1] = ' ';
	  hit_ptr->alignment[2][length - 1] = '-';

	  if (remove)
	    {
	      m1[i][j] = 0;
	    }
      } else
	{

	  if (m2[i][j - 1] == LEFT)
	    {
	      if (i < (strlen (sequence1) - 10))
		{
		  hit_ptr->score += gap_extend;
	      } else
		{
		  hit_ptr->score += (gap_extend * scale);
		}
	  } else
	    {
	      if (i < (strlen (sequence1) - 10))
		{
		  hit_ptr->score += gap_open;
	      } else
		{
		  hit_ptr->score += (gap_open * scale);
		}
	    }

	  length++;
	  traceback (m1, m2, sequence1, sequence2, i, j - 1, 1, hit_ptr,
		     length);
	  hit_ptr->alignment[0][length - 1] = '-';
	  hit_ptr->alignment[1][length - 1] = ' ';
	  hit_ptr->alignment[2][length - 1] = sequence2[j - 1];

	  if (remove)
	    {
	      m1[i][j] = 0;
	    }
	}
  } else
    {
      hit_ptr->query_start = i;
      hit_ptr->ref_start = j;
      hit_ptr->alignment[0][length] = '\0';
      hit_ptr->alignment[1][length] = '\0';
      hit_ptr->alignment[2][length] = '\0';
    }
}

int
build_sub_matrix (int **matrix)
{

  int i = 0;
  int j = 0;

  for (i = 0; i < 6; i++)
    {
      for (j = 0; j < 6; j++)
	{
	  matrix[toupper (baselist[i])][toupper (baselist[j])] =
	    (int) match[i][j];
	}
    }

  return (1);
}


double
score (char nt1, char nt2)
{
  return (double) match[bases[(int) nt1]][bases[(int) nt2]];
}
