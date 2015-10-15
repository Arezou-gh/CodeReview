#ifndef _SCHURCOMPLEMENT_HH
#define _SCHURCOMPLEMENT_HH

#include <deal.II/base/smartpointer.h>
#include <deal.II/lac/block_sparse_matrix.h>
#include <deal.II/lac/sparse_direct.h>
#include <deal.II/lac/vector.h> 

using namespace dealii;

class SchurComplement : public Subscriptor
{
  public:
    SchurComplement (const BlockSparseMatrix<double> &system_matrix,
        const SparseDirectUMFPACK &A_inverse);
    void vmult (Vector<double>       &dst,
        const Vector<double> &src) const;
  private:
    const SmartPointer<const BlockSparseMatrix<double> > system_matrix;
    const SmartPointer<const SparseDirectUMFPACK> A_inverse;
    mutable Vector<double> tmp1, tmp2;
};

#endif
