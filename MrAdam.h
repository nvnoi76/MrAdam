#ifndef _MRADAM_H
#define _MRADAM_H

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h> 

#define DBL_RAND_MAX (double)(RAND_MAX)
#define M_PI   3.14159265358979323846

//So tham so hoi quy
#define NUMPARAM  20

#define	A11	p[0]
#define	A12	p[1]
#define	A2	p[2]
#define	A3	p[3]
#define	A4	p[4]
#define	A5	p[5]
#define	A6	p[6]
#define	A7	p[7]
#define	A8	p[8]
#define	A9	p[9]
#define	A10	p[10]
#define	B11	p[11]
#define	B12	p[12]
#define	B2	p[13]
#define	B3	p[14]
#define	B4	p[15]
#define	B5	p[16]
#define	B6	p[17]
#define	B7	p[18]
#define	B8	p[19]

//So tham so bang du lieu
#define NUMCOL	  14
//So dong bang du lieu
#define NUMROW	  315

extern double gdata[NUMROW][NUMCOL];

#define	MR		gdata[i][0]
#define	PA		gdata[i][1]
#define	WOPT	gdata[i][2]
#define	WW		gdata[i][3]
#define	XM3		gdata[i][4]
#define	SS		gdata[i][5]
#define	QU		gdata[i][6]
#define	PI		gdata[i][7]
#define	LL		gdata[i][8]
#define	N200	gdata[i][9]
#define	XMD		gdata[i][10]
#define	MAU		gdata[i][11]
#define	SUBMAU	gdata[i][12]
#define	ID		gdata[i][13]
#endif 
