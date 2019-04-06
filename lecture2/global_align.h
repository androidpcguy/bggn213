/*
*/

#ifndef __GLOBAL_ALIGN_H_
#define __GLOBAL_ALIGN_H_

int const MATCH_SCORE		=  2;
int const MISMATCH_SCORE	= -1;
int const GAP_SCORE		= -2;

typedef struct matrix_element {
	int score;
	unsigned char direction;
} matrix_element;

void print_alignment();
void calculate_cell_score();

#endif
