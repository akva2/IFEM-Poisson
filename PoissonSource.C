// $Id$
//==============================================================================
//!
//! \file PoissonSource.C
//!
//! \date Jan 29 2025
//!
//! \author Arne Morten Kvarving / SINTEF
//!
//! \brief Class for Poisson source function.
//!
//==============================================================================

#include "Poisson.h"
#include "PoissonSource.h"

#include "AnaSol.h"
#include "TensorFunction.h"


PoissonAnaSolSource::PoissonAnaSolSource (const AnaSol& aSol,
                                          const Poisson& prob)
  : anaSol(aSol), poisson(prob)
{
}


double PoissonAnaSolSource::evaluate (const Vec3& X) const
{
  return -poisson.getMaterial(X) * anaSol.getScalarSol()->hessian(X).trace();
}


PoissonAnaSolSourceVec::PoissonAnaSolSourceVec (const AnaSol& aSol,
                                                const Poisson& prob)
    : anaSol(aSol), poisson(prob)
{
  ncmp = (*aSol.getVectorSecSol())(Vec3()).dim();
}


Vec3 PoissonAnaSolSourceVec::evaluate (const Vec3& X) const
{
  Vec3 res;
  const utl::matrix3d<Real> hess = anaSol.getVectorSol()->hessian(X);
  for (size_t i = 0; i < ncmp; ++i)
    res[i] = -hess.trace(i+1);

  return res;
}
