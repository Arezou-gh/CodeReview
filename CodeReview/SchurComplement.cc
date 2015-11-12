#include "SchurComplement.hh"

SchurComplement::
SchurComplement (const BlockSparseMatrix<double> &system_matrix,
    const SparseDirectUMFPACK &A_inverse)
  :
    system_matrix (&system_matrix),
    A_inverse (&A_inverse),
    tmp1 (system_matrix.block(0,0).m()),
    tmp2 (system_matrix.block(0,0).m())
{}


void SchurComplement::vmult (Vector<double> &dst,
    const Vector<double> &src) const
{
  system_matrix->block(0,1).vmult (tmp1, src);
  A_inverse->vmult (tmp2, tmp1);
  system_matrix->block(1,0).vmult (dst, tmp2);
}
