// $Id$
//==============================================================================
//!
//! \file PoissonSource.h
//!
//! \date Jan 4 2024
//!
//! \author Arne Morten Kvarving / SINTEF
//!
//! \brief Class for Poisson source function.
//!
//==============================================================================
#ifndef POISSON_SOURCE_H_
#define POISSON_SOURCE_H_

#include "Function.h"

#include "Poisson.h"

#include <vector>

class AnaSol;


/*!
  \brief Class that derives the Poisson source function from the analytic solution.
 */

class PoissonAnaSolSource : public RealFunc
{
public:
  //! \brief Constructor.
  //! \param aSol Analytic solution to use
  //! \param material Material coefficients
  PoissonAnaSolSource(const AnaSol& aSol,
                      const std::vector<Poisson::Kappa>& material);

  //! \brief Constructor.
  //! \param aSol Analytic solution to use
  //! \param kappa Constant material coefficient
  PoissonAnaSolSource(const AnaSol& aSol,
                      const double kappa);

  bool initPatch(size_t pidx) override { pIdx = pidx; return true; }

protected:
  //! \brief Evaluates the function.
  Real evaluate(const Vec3& X) const override;

  //! \brief Returns material coeffient for point.
  Real getKappa(const Vec3& X) const;

  const AnaSol& anaSol; //!< Reference to analytic solution
  const double kappa; //!< Constant material coefficient
  const std::vector<Poisson::Kappa>* materials; //!< Material coefficients
  size_t pIdx = 0; //!< Current patch index
};

#endif
