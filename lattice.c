#include "lattice.h"


void _set_Q(Lattice *pLat) {

   char info[80];

   assert((MIN_LATTICE_EDGE_LEN - EPS < pLat->edgeLen) && (pLat->edgeLen < MAX_LATTICE_EDGE_LEN + EPS));

   double edgeLen = pLat->edgeLen;

   int ptsOnEdge = (int)(edgeLen * LATTICE_DENSITY);

   pLat->Q = ptsOnEdge * ptsOnEdge * ptsOnEdge;

sprintf(info, "lattice edge length: %2.1f, lattice density: %d, points on edge: %d, lat.Q=%d\n", edgeLen, LATTICE_DENSITY, ptsOnEdge, pLat->Q);
write_to_log(info);
}

void _set_p(Lattice *pLat) {

   assert((MIN_LATTICE_EDGE_LEN - EPS < pLat->edgeLen) && (pLat->edgeLen < MAX_LATTICE_EDGE_LEN + EPS));
   assert(pLat->Q > 0);

   double edgeLen = pLat->edgeLen;

   int q = pLat->Q;

   pLat->p = (double**)malloc(sizeof(double*)*q);
   assert(pLat->p != NULL);

   for (int i = 0; i < q; i++) {
       pLat->p[i] = (double*)malloc(sizeof(double)*3);
       assert(pLat->p[i] != NULL);
   }

   double step = 1./(LATTICE_DENSITY - 1);
   double halfEdgeLen = edgeLen/2.;

   int ptsOnEdge = (int)(edgeLen * LATTICE_DENSITY);

   for (int i = 0; i < q; i++) {
	
         int indX = i % ptsOnEdge;
         int indY = (i % (ptsOnEdge * ptsOnEdge)) / ptsOnEdge;
         int indZ = i / (ptsOnEdge * ptsOnEdge); 

         pLat->p[i][0] = indX * step - halfEdgeLen; 
         pLat->p[i][1] = indY * step - halfEdgeLen; 
         pLat->p[i][2] = indZ * step - halfEdgeLen; 
   }
}

void _set_vscores(Lattice *pLat) {

   assert(pLat->Q > 0);
   assert(pLat->p != NULL);

   int q = pLat->Q;

   for (int i = 0; i < q; i++)
      assert(pLat->p[i] != NULL);

   pLat->vscores = (double*)malloc(sizeof(double) * q);
   assert(pLat->vscores != NULL);

   for (int i = 0; i < q; i++) 
      pLat->vscores[i] = pLat->p[i][0] * pLat->p[i][0] + pLat->p[i][1] * pLat->p[i][1] + pLat->p[i][2] * pLat->p[i][2];
}


Lattice init_lattice(double edgeLen) {

   char info[80];

   assert((MIN_LATTICE_EDGE_LEN - EPS < edgeLen) && (edgeLen < MAX_LATTICE_EDGE_LEN + EPS));

   Lattice lat;

   lat.edgeLen = edgeLen;

   lat.Q = -1;

   _set_Q(&lat);
   _set_p(&lat);
   _set_vscores(&lat);

   return lat;
}

void free_lattice(Lattice lat) {
  
   if (lat.vscores != NULL) free(lat.vscores);

   if (lat.p == NULL) {
       lat.Q = -1;
        lat.edgeLen = -1.;
       return;
   }

   for (int i = 0; i < lat.Q; i++) 
       if (lat.p[i] != NULL) free(lat.p[i]);

   free(lat.p);

   lat.Q = -1;

   lat.edgeLen = -1.;
}
