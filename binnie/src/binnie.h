/*
 * binnie.h - constants shared by all binnie code
 *
 * Copyright (c) 2013 Genome Research Ltd. 
 * Author: Joshua C. Randall <jcrandall@alum.mit.edu>
 *
 * This file is part of BridgeBuilder. 
 *
 * BridgeBuilder is free software: you can redistribute it and/or modify it under 
 * the terms of the GNU General Public License as published by the Free Software 
 * Foundation; either version 3 of the License, or (at your option) any later 
 * version. 
 * 
 * This program is distributed in the hope that it will be useful, but WITHOUT 
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS 
 * FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.  
 * 
 * You should have received a copy of the GNU General Public License along with 
 * this program. If not, see <http://www.gnu.org/licenses/>.
 *
 */
#ifndef BINNIE_H
#define BINNIE_H


/* gnulib headers */
#include <stdbool.h>
#include "size_max.h" /* to ensure SIZE_MAX is available */


/* 
 * if true, don't compare RG tags between original and bridge when matching reads
 * (i.e. in case they are missing from the bridge-mapped file) 
 */
bool ignore_rg;


/* option defaults */
#define BINNIE_DEFAULT_BUFFER_SIZE  1000000
#define BINNIE_DEFAULT_BUFFER_BASES 10000


/* hash parameters */
#define BINNIE_TABLESIZE SIZE_MAX


/* output bins */
#define BINNIE_UNCHANGED    0
#define BINNIE_BRIDGED      1
#define BINNIE_REMAP        2


/* exit codes */
#define BINNIE_EXIT_SUCCESS           	 0
#define BINNIE_EXIT_ERR_ARGS          	 1
#define BINNIE_EXIT_ERR_IN_FILES      	 2
#define BINNIE_EXIT_ERR_OUT_FILES     	 3 
#define BINNIE_EXIT_ERR_UID          	 4
#define BINNIE_EXIT_ERR_READ_ORIG     	 5
#define BINNIE_EXIT_ERR_READ_BRIDGE   	 6
#define BINNIE_EXIT_ERR_SEGMENT_INDEX 	 7
#define BINNIE_EXIT_ERR_ORIG_TRUNCATED	 8
#define BINNIE_EXIT_ERR_UNEXPECTED_MATES 9
#define BINNIE_EXIT_ERR_NULL         	 10
#define BINNIE_EXIT_ERR_NOT_NULL         11
#define BINNIE_EXIT_ERR_BUFFER_NOT_EMPTY 12
#define BINNIE_EXIT_ERR_BAM_UNSORTED     13
#define BINNIE_EXIT_ERR_INVALID_BIN      14
#define BINNIE_EXIT_ERR_WRITE            15
#define BINNIE_EXIT_ERR_BUFFER_REMOVE    16
#define BINNIE_EXIT_ERR_BRIDGE_SORT	 17

#endif
