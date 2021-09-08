/******************************************************************************
 * I/O code.
 *
 * Copyright 2017, 2019, Alexander Jones.
 *
 * This program is free software; you can redistribute it and/or modify it under
 * the terms of the GNU General Public License as published by the Free Software
 * Foundation; either version 3 of the License, or (at your option) any later
 * version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
 * or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public
 * License for more details.
 *
 * You should have received a copy of the GNU General Public License along with
 * this program; see the file LICENSE.  If not, see http://www.gnu.org/licenses/
 * or write to the Free Software Foundation, Inc., 51 Franklin St, Fifth Floor,
 * Boston, MA 02110-1301, USA.
 ******************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "rho.h"

int readComposites(Composite **composites, char *filename) {
	Composite *newComposites;
	Composite *comp;
	int length;
	void *tmp;
	char num[MAX_NUM_SIZE];
	int aliquotSeq, aliquotTerm, factor;
	int floyd, brent, yafu, montgomery;
	FILE *compositeFile = fopen(filename, "r");
	fscanf(compositeFile, "%d", &length);
	tmp = malloc(length*sizeof(Composite));
	if (!tmp) {
		exit(1);
	}
	newComposites = (Composite *)tmp;
	for (int i = 0; i < length; i++) {
		comp = &(newComposites[i]);
		fscanf(compositeFile, "%s", num);
		fscanf(compositeFile, "%d:%d %d", &aliquotSeq, &aliquotTerm, &factor);
		strncpy(comp->number, num, MAX_NUM_SIZE);
		comp->aliquotSeq = aliquotSeq;
		comp->aliquotTerm = aliquotTerm;
		for (int j = 0; j < NUM_POLYS; j++) {
			fscanf(compositeFile, "%d %d %d %d", &floyd, &brent, &yafu, &montgomery);
			comp->maxLengths[j].floyd = floyd;
			comp->maxLengths[j].brent = brent;
			comp->maxLengths[j].yafu = yafu;
			comp->maxLengths[j].montgomery = montgomery;
		}
	}
	fclose(compositeFile);
	*composites = newComposites;
    return length;
}
