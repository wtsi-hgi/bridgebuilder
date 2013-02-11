/*
 * binnie_files.h - handles BAM/SAM input and output (via htslib)
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

#ifndef _BINNIE_FILES_H
#define _BINNIE_FILES_H 1

#include <htslib/sam.h>

htsFile *binnie_open_out(const char *filename);
htsFile *binnie_open_in(const char *filename);
void binnie_close(htsFile *fp);

#endif /* _BINNIE_FILES_H */
