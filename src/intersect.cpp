#include <cstdlib>
#include <cmath>
#include <cassert>
#include <cstring>
#include <iostream>
#include <fstream>
#include "node.h"
#include "elt.h"
#include "grid.h"
#include "inside.h"
#include "polyg.h"

#include "intersect.h"

using namespace std;

/** returns index of edge of a that is shared with b,
    or NOT_FOUND if a and b do not share an edge */
int neighbour_idx(const Elt& a, const Elt& b)
{
	for (int i = 0; i < a.n; i++)
	{
		for (int j = 0; j < b.n; j++)
		{
			assert(squaredist(a.vertex[ i       ], b.vertex[ j       ]) > EPS*EPS ||
			       squaredist(a.vertex[(i+1)%a.n], b.vertex[(j+1)%b.n]) > EPS*EPS);
			if (   squaredist(a.vertex[ i       ], b.vertex[ j           ]) < 1e-13*1e-13 &&
			       squaredist(a.vertex[(i+1)%a.n], b.vertex[(j+b.n-1)%b.n]) < 1e-13*1e-13)
			{
				return i;
			} 
		}
	}
	return NOT_FOUND;
}


/** 
If `a` and `b` are neighbours (in the sense that they share an edge) 
then this information will be stored in `a` (but not in `b`)
*/
void set_neighbour(Elt& a, const Elt& b)
{
	if (b.id.ind == a.id.ind) return;
	int idx = neighbour_idx(a, b);
	if (idx != NOT_FOUND)
		a.neighbour[idx] = b.id.ind;
}

/** return true if `a` and `b` share an edge */
bool isNeighbour(const Elt& a, const Elt& b)
{
	return neighbour_idx(a, b) != NOT_FOUND; 
}

/* computes intersection between elements a and b */
void intersect(Elt *a, Elt *b)
{
	int na = a->n; /* vertices of a */
	int nb = b->n; /* vertices of b */
	Coord *c   = new Coord[na+nb];
	Coord *c2  = new Coord[na+nb];
	Coord *xc  = new Coord[na+nb];
	Coord *xc2 = new Coord[na+nb];
	Coord gc, gc2;
	double *d = new double[na+nb];
	double *d2 = new double[na+nb];
	double are, are2;
	Ipt ipt[NMAX*NMAX];
	Ipt ipt2[NMAX*NMAX];
	ptsec(a, b, ipt);
	/* make ipt2 transpose of ipt */
	for (int ii = 0; ii < na; ii++)
		for (int jj = 0; jj < nb; jj++)
			ipt2[jj*na+ii] = ipt[ii*nb+jj];
	list<Sgm> iscot;
	recense(a, b, ipt, iscot, 0);
	recense(b, a, ipt2, iscot, 1);

	int nseg = iscot.size();
	int nc = 0;
	int nc2 = 0;
	while (iscot.size() && nc < 2)
		nc = assemble(iscot, c, d, xc);
	while (iscot.size() && nc2 < 2)
		nc2 = assemble(iscot, c2, d2, xc2);
	/*assert(nseg == nc + nc2 || nseg == 1); // unused segment*/
        if (not(nseg == nc + nc2 || nseg == 1))
        {
cerr << "WARNING. NOT assert(nseg == nc + nc2 || nseg == 1). ";
        }
	if (nc == 1) nc = 0;
	if (nc2 == 1) nc2 = 0;
	gc = barycentre(xc, nc);
	gc2 = barycentre(xc2, nc2);
	orient(nc, xc, c, d, gc);

	Coord pole = srcGrid.pole;
	if (pole == ORIGIN) pole = tgtGrid.pole;
	const double MINBASE = 1e-11;
	if (nc == 2) /* nc is the number of vertices of super mesh element */
	{
		double base = arcdist(xc[0], xc[1]);
cerr << "DID ARRIVE " << base << xc[0] << xc[1] << endl;
		gc = midpoint(gc, midpointSC(xc[0], xc[1]));
		/* intersection area `are` must be zero here unless we have one great and one small circle */
		are = alun(base, fabs(scalarprod(xc[0], pole)));
	}
	else
	{
		are = airbar(nc, xc, c, d, pole, gc);
	}
	if (nc2 == 2)
	{
		double base = arcdist(xc2[0], xc2[1]);
cerr << "DID ARRIVE " << base << xc2[0] << xc2[1] << endl;
		assert(base > MINBASE);
		gc2 = midpoint(gc2, midpointSC(xc2[0], xc2[1]));
		are2 = alun(base, fabs(scalarprod(xc2[0], pole))); // 0
	}
	else
	{
		are2 = airbar(nc2, xc2, c2, d2, pole, gc2);
	}

	if (nc > 1)
	{
		/* create one super mesh polygon that src and dest point to */
		Polyg *is = new Polyg;
		is->x = gc;
		is->area = are;
		is->id = b->id;
		is->src_id = b->src_id;
		is->n = nc;
		(a->is).push_back(is);
		(b->is).push_back(is);
	}
	if (nc2 > 1)
	{
		Polyg *is = new Polyg;
		is->x = gc2;
		is->area = are2;
		is->id = b->id; /* intersection holds id of corresponding source element (see Elt class definition for details about id) */
		is->src_id = b->src_id;
		is->n = nc2;
		(a->is).push_back(is);
		(b->is).push_back(is);
	}
	delete [] c;
	delete [] c2;
	delete [] xc;
	delete [] xc2;
	delete [] d;
	delete [] d2;
}
