/*
 * brindley.h - constants shared by all brindley code
 *
 * Copyright (c) 2014 Genome Research Ltd. 
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
#ifndef BRINDLEY_H
#define BRINDLEY_H


/* gnulib headers */
#include <stdbool.h>
#include "size_max.h" /* to ensure SIZE_MAX is available */


/* hash parameters */
#define BRINDLEY_TABLESIZE 4294967296 // 2^32


/* exit codes */
#define BRINDLEY_EXIT_SUCCESS           	 0
#define BRINDLEY_EXIT_ERR_ARGS          	 1
#define BRINDLEY_EXIT_ERR_IN_FILES      	 2
#define BRINDLEY_EXIT_ERR_OUT_FILES     	 3 
#define BRINDLEY_EXIT_ERR_WRITE            15

#endif
