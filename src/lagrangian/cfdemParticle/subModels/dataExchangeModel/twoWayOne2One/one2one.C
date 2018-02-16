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

//#define O2ODEBUG

#ifdef O2ODEBUG
#include <iostream>
#endif

#include <mpi.h>
#include "one2one.H"

template <typename T>
struct mpi_type_wrapper {
  MPI_Datatype mpi_type;
  mpi_type_wrapper();
};
template <> mpi_type_wrapper<float>::mpi_type_wrapper()
: mpi_type(MPI_FLOAT) {}
template <> mpi_type_wrapper<double>::mpi_type_wrapper()
: mpi_type(MPI_DOUBLE) {}
template <> mpi_type_wrapper<int>::mpi_type_wrapper()
: mpi_type(MPI_INT) {}

/* ---------------------------------------------------------------------- */

One2One::One2One(MPI_Comm caller)
:
ncollected_(-1),
comm_(caller),
nsrc_procs_(-1),
src_procs_(nullptr),
ndst_procs_(-1),
dst_procs_(nullptr),
nlocal_(-1),
natoms_(nullptr),
request_(nullptr),
status_(nullptr)
{
  MPI_Comm_rank(comm_,&me_);
  MPI_Comm_size(comm_,&nprocs_);
}

/* ---------------------------------------------------------------------- */

One2One::~One2One()
{
  deallocate();
}

/* ----------------------------------------------------------------------
  communicate particle ids based on processor communication pattern
------------------------------------------------------------------------- */

void One2One::setup
(
  int nsrc_procs,
  int *src_procs,
  int ndst_procs,
  int *dst_procs,
  int nlocal
)
{
  // free any previous one2one info
  deallocate();

  src_procs_ = src_procs;
  nsrc_procs_ = nsrc_procs;
  dst_procs_ = dst_procs;
  ndst_procs_ = ndst_procs;
  nlocal_ = nlocal;

  // gather number of ids for reserving memory
  natoms_ = new int[nprocs_];
  MPI_Allgather // may be replaced by send/irecv
  (
    &nlocal_,
    1,
    MPI_INT,
    natoms_,
    1,
    MPI_INT,
    comm_
  );

  ncollected_ = 0;
  for (int i = 0; i < nsrc_procs_; i++)
    ncollected_ += natoms_[src_procs_[i]];

  request_ = new MPI_Request[nsrc_procs_];
  status_ = new MPI_Status[nsrc_procs_];
}

/* ----------------------------------------------------------------------
   src: what is present on this proc
   dst: what is received from other procs
   all comm according to map set up in setup(...)
------------------------------------------------------------------------- */
template  <typename T>
void One2One::exchange(T *&src, T *&dst, int data_length)
{
  mpi_type_wrapper<T> wrap;

  // post receives
  int offset = 0;
  for (int i = 0; i < nsrc_procs_; i++)
  {
  #ifdef O2ODEBUG
  std::cout<< "[" << me_ << "]"
           << " RCV " << i
           << " of "  << nsrc_procs_
           << " from: " << src_procs_[i]
           << " natoms_[src_procs_[i]] " << natoms_[src_procs_[i]]
           << " datalength " << data_length
           << " offset " << offset
           << std::endl;
  #endif
    MPI_Irecv
    (
      dst+offset,
      natoms_[src_procs_[i]]*data_length,
      wrap.mpi_type,
      src_procs_[i],
      MPI_ANY_TAG,
      comm_,
      &request_[i]
    );
    offset += natoms_[src_procs_[i]]*data_length;
  }

  // make sure all receives are posted
  MPI_Barrier(comm_);

  // blocking sends - do nonblocking instead
  //                  since doing many-2-many here?
  for (int i = 0; i < ndst_procs_; i++)
  {
  #ifdef O2ODEBUG
  std::cout<< "[" << me_ << "]"
           << " SEND to: " << dst_procs_[i]
           << " nlocal_ " << nlocal_
           << " data_length " << data_length
           << std::endl;
  #endif
    MPI_Send
    (
      src,
      nlocal_*data_length,
      wrap.mpi_type,
      dst_procs_[i],
      0,
      comm_
    );
  }

  MPI_Waitall(nsrc_procs_, request_, status_);

  // copy on-proc info instead of communicating everything
  // may yield a first step towards true one2one communication
  //for (i = 0; i < nown; i++)
  //    dest[dest_own[i]] = src[src_own[i]];

}

template void One2One::exchange<int>(int*&, int*&, int);
template void One2One::exchange<double>(double*&, double*&, int);

// there should be a way to do this without copying data
template  <typename T>
void One2One::exchange(T **&src, T **&dst, int data_length)
{
  mpi_type_wrapper<T> wrap;

  T* tmp_dst = new T[ncollected_*data_length];
  T* tmp_src = new T[nlocal_*data_length];

  for (int i = 0; i < nlocal_; i++)
    for (int j = 0; j < data_length; j++)
      tmp_src[data_length*i+j] = src[i][j];

  exchange<T>(tmp_src, tmp_dst, data_length);

  for (int i = 0; i < ncollected_; i++)
    for (int j = 0; j < data_length; j++)
      dst[i][j] = tmp_dst[data_length*i+j];

  delete [] tmp_src;
  delete [] tmp_dst;
}
template void One2One::exchange<int>(int**&, int**&, int);
template void One2One::exchange<double>(double**&, double**&, int);


template  <typename T>
void One2One::exchange(T **&src, T *&dst, int data_length)
{
  mpi_type_wrapper<T> wrap;

  T* tmp_src = new T[nlocal_*data_length];

  for (int i = 0; i < nlocal_; i++)
    for (int j = 0; j < data_length; j++)
      tmp_src[data_length*i+j] = src[i][j];

  exchange<T>(tmp_src, dst, data_length);


  delete [] tmp_src;
}
template void One2One::exchange<int>(int**&, int*&, int);
template void One2One::exchange<double>(double**&, double*&, int);

/* ---------------------------------------------------------------------- */

void One2One::deallocate()
{
  delete [] src_procs_;
  delete [] dst_procs_;
  delete [] natoms_;

  delete [] request_;
  delete [] status_;
}