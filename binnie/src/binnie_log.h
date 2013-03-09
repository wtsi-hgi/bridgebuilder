/*
 * binnie_log.h - logging functions
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
#ifndef BINNIE_LOG_H
#define BINNIE_LOG_H

#include <stdbool.h>


/* verbosity level (0-3; increasing with each -v option): 0 is silent, 3 is maximum verbosity */
unsigned int verbosity;


void blog(unsigned int level, const char *msgfmt, ...);


#ifdef DEBUG
/* debug flag: if true, print debugging messages to stderr */
bool debug_flag;
void DLOG(const char *msgfmt, ...);
#else
/* (void)sizeof will slurp up variadic functions and don't get evaluated at run-time */
#define DLOG (void)sizeof
#endif


#endif
