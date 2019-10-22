/*---------------------------------------------------------------------------*\
License

    This is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This code is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with this code.  If not, see <http://www.gnu.org/licenses/>.

Copyright 2018-     Paul Kieckhefen, TUHH

\*---------------------------------------------------------------------------*/

#ifndef ONE2ONE_H
#define ONE2ONE_H

#include <mpi.h>

class One2One
{
 public:
  One2One(MPI_Comm);
  ~One2One();

  void setup
  (
    int nsrc_procs,
    int *src_procs,
    int ndst_procs,
    int* dst_procs,
    int nlocal
  );

  template <typename T>
  void exchange(T *&, T *&, int data_length=1);
  template <typename T>
  void exchange(T **&, T **&, int data_length=1);
  template <typename T>
  void exchange(T **&, T *&, int data_length=1);

  int ncollected_;            // # of ids in from group

 protected:

  MPI_Comm comm_;
  int me_, nprocs_;                // rank and size

  // communication partners
  int  nsrc_procs_;                // # of off-processor IDs
  int* src_procs_;                 // procs I receive data from
  int  ndst_procs_;                // # of off-processor IDs
  int* dst_procs_;                 // procs I receive data from

  int nlocal_;                     // # particle ids I own
  int* natoms_;

  MPI_Request* request_;
  MPI_Status* status_;

  void deallocate();
};

#endif
