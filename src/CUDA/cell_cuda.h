#ifndef SWIFT_CUDA_CELL_H
#define SWIFT_CUDA_CELL_H

#define CUDA_MAX_LINKS 27

struct cell_cuda {
  /* The cell location on the grid. */
  double loc[3];

  /* The cell dimensions. */
  double width[3];

  /* Max smoothing length in this cell. */
  double h_max;

  /* Index of the particle data. */
  int first_part;
  
  /* Number of particles in the cell. */
  int part_count;

  /* Indices of the next level of cells. */
  int progeny[8];

  /* Index of the parent cell. */
  int parent;

  /* Index of the super cell.*/
  int super;

  /* Minimum end of time step in this cell. */
  integertime_t ti_end_min;

  /* Maximum end of time step in this cell. */
  integertime_t ti_end_max;
  
  /* Minimum dimension of this cell */
  float dmin;

  /* Need the density links only, maximum of 27 */
  int links[CUDA_MAX_LINKS];

  /* Number of links */  
  int nr_links;

  /* IS split? */
  int split;
  
  float h_max_old;
  float dx_max_old;

};

#endif /* SWIFT_CUDA_CELL_H */
