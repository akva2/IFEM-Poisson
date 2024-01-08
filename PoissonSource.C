// $Id$
//==============================================================================
//!
//! \file PoissonSource.C
//!
//! \date Jan 4 2024
//!
//! \author Arne Morten Kvarving / SINTEF
//!
//! \brief Class for Poisson source function.
//!
//==============================================================================

#include "PoissonSource.h"

#include "AnaSol.h"


PoissonAnaSolSource::PoissonAnaSolSource (const AnaSol& aSol,
                                          const std::vector<Poisson::Kappa>& material) :
  anaSol(aSol), kappa(0.0), materials(&material)
{
  ncmp = 1;
}


PoissonAnaSolSource::PoissonAnaSolSource (const AnaSol& aSol,
                                          const double kapp) :
  anaSol(aSol), kappa(kapp), materials(nullptr)
{
  ncmp = 1;
}


Real PoissonAnaSolSource::evaluate (const Vec3& X) const
{
  const SymmTensor hess = anaSol.getScalarSol()->hessian(X);
  return -this->getKappa(X) * hess.trace();
}


Real PoissonAnaSolSource::getKappa (const Vec3& X) const
{
  if (!materials)
    return kappa;

  if (materials->at(pIdx).func)
    return (*materials->at(pIdx).func)(X);

  return materials->at(pIdx).constant;
}
