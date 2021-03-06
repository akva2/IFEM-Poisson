// $Id$
//==============================================================================
//!
//! \file main_Poisson.C
//!
//! \date 20 May 2010
//!
//! \author Einar Christensen / SINTEF
//!
//! \brief Main program for the isogeometric solver for the Poisson equation.
//!
//==============================================================================

#include "IFEM.h"
#include "SIMPoisson.h"
#include "AdaptiveSIM.h"
#include "HDF5Writer.h"
#include "XMLWriter.h"
#include "Utilities.h"
#include "Profiler.h"
#include "AppCommon.h"
#include "SIMSolverAdap.h"
#include "PoissonArgs.h"
#include <stdlib.h>
#include <string.h>
#include <ctype.h>


/*!
  \brief Sets up and launches the simulation.
  \param[in] infile The input file to process
  \param[in] checkRHS If \e true, check patches for a right-hand-side model
  \param[in] ignoredPatches List of patches to ignore in the simulation
  \param[in] fixDup If \e true, collapse co-located nodes into a single node
  \param[in] vizRHS If \e true, save the load vector to VTF for visualization
*/

template<class Dim, template<class T> class Solver=SIMSolver>
int runSimulator(char* infile, bool checkRHS,
                 const std::vector<int>& ignoredPatches,
                 bool fixDup, bool vizRHS, bool dumpASCII)
{
  SIMPoisson<Dim> model(checkRHS);
  Solver<SIMPoisson<Dim>> solver(model);

  if (!model.read(infile) || !solver.read(infile))
    return 1;

  // Boundary conditions can be ignored only in generalized eigenvalue analysis
  if (model.opt.eig != 4 && model.opt.eig != 6)
    SIMbase::ignoreDirichlet = false;

  // Load vector visualization is not available when using additional viz-points
  for (size_t i = 0; i < 3; i++)
    if (i >= Dim::dimension)
      model.opt.nViz[i] = 1;
    else if (model.opt.nViz[i] > 2)
      vizRHS = false;

  if (!model.preprocess(ignoredPatches, fixDup))
    return 2;

  model.setVizRHS(vizRHS);
  model.setDumpASCII(dumpASCII);

  model.setQuadratureRule(model.opt.nGauss[0],true);
  model.initSystem(model.opt.solver,1,model.opt.eig>0?0:1,0,true);

  std::unique_ptr<DataExporter> exporter;
  if (model.opt.dumpHDF5(infile)) {
    if (model.opt.discretization < ASM::Spline && !model.opt.hdf5.empty())
      IFEM::cout <<"\n ** HDF5 output is available for spline discretization only"
        <<". Deactivating...\n"<< std::endl;
    else
      exporter.reset(SIM::handleDataOutput(model, solver, model.opt.hdf5, false, 1, 1));
  }

  return solver.solveProblem(infile, exporter.get(), "Solving Poisson problem", false);
}


/*!
  \brief Main program for the NURBS-based isogeometric Poisson equation solver.

  The input to the program is specified through the following
  command-line arguments. The arguments may be given in arbitrary order.

  \arg \a input-file : Input file with model definition
  \arg -dense :   Use the dense LAPACK matrix equation solver
  \arg -spr :     Use the SPR direct equation solver
  \arg -superlu : Use the sparse SuperLU equation solver
  \arg -samg :    Use the sparse algebraic multi-grid equation solver
  \arg -petsc :   Use equation solver from PETSc library
  \arg -lag : Use Lagrangian basis functions instead of splines/NURBS
  \arg -spec : Use Spectral basis functions instead of splines/NURBS
  \arg -LR : Use LR-spline basis functions instead of tensorial splines/NURBS
  \arg -nGauss \a n : Number of Gauss points over a knot-span in each direction
  \arg -vtf \a format : VTF-file format (-1=NONE, 0=ASCII, 1=BINARY)
  \arg -nviz \a nviz : Number of visualization points over each knot-span
  \arg -nu \a nu : Number of visualization points per knot-span in u-direction
  \arg -nv \a nv : Number of visualization points per knot-span in v-direction
  \arg -nw \a nw : Number of visualization points per knot-span in w-direction
  \arg -hdf5 : Write primary and projected secondary solution to HDF5 file
  \arg -dumpASC : Dump model and solution to ASCII files for external processing
  \arg -ignore \a p1, \a p2, ... : Ignore these patches in the analysis
  \arg -eig \a iop : Eigenproblem solver to use (1...6)
  \arg -nev \a nev : Number of eigenvalues to compute
  \arg -ncv \a ncv : Number of Arnoldi vectors to use in the eigenvalue analysis
  \arg -shift \a shf : Shift value to use in the eigenproblem solver
  \arg -free : Ignore all boundary conditions (use in free vibration analysis)
  \arg -check : Data check only, read model and output to VTF (no solution)
  \arg -checkRHS : Check that the patches are modelled in a right-hand system
  \arg -vizRHS : Save the right-hand-side load vector on the VTF-file
  \arg -fixDup : Resolve co-located nodes by merging them into a single node
  \arg -1D : Use one-parametric simulation driver
  \arg -2D : Use two-parametric simulation driver
  \arg -adap : Use adaptive simulation driver with LR-splines discretization
*/

int main (int argc, char** argv)
{
  Profiler prof(argv[0]);
  utl::profiler->start("Initialization");

  std::vector<int> ignoredPatches;
  bool checkRHS = false;
  bool vizRHS = false;
  bool fixDup = false;
  bool dumpASCII = false;
  char* infile = nullptr;
  PoissonArgs args;

  int myPid = IFEM::Init(argc,argv,"Poisson solver");
  int ignoreArg = -1;
  for (int i = 1; i < argc; i++)
    if (i == ignoreArg || SIMoptions::ignoreOldOptions(argc,argv,i))
      ; // ignore the obsolete option
    else if (!strcmp(argv[i],"-dumpASC"))
      dumpASCII = myPid == 0; // not for parallel runs
    else if (!strcmp(argv[i],"-ignore"))
      while (i < argc-1 && isdigit(argv[i+1][0]))
        utl::parseIntegers(ignoredPatches,argv[++i]);
    else if (!strcmp(argv[i],"-free"))
      SIMbase::ignoreDirichlet = true;
    else if (!strcmp(argv[i],"-checkRHS"))
      checkRHS = true;
    else if (!strcmp(argv[i],"-vizRHS"))
      vizRHS = true;
    else if (!strcmp(argv[i],"-fixDup"))
      fixDup = true;
    else if (!strcmp(argv[i],"-1D"))
      args.dim = 1;
    else if (!strcmp(argv[i],"-2D"))
      args.dim = 2;
    else if (!strncmp(argv[i],"-adap",5))
      args.adap = true;
    else if (!infile) {
      infile = argv[i];
      ignoreArg = i;
      if (strcasestr(infile,".xinp")) {
        if (!args.readXML(infile,false))
          return 1;
        i = 0;
      }
    } else
      std::cerr <<"  ** Unknown option ignored: "<< argv[i] << std::endl;

  if (!infile)
  {
    std::cout <<"usage: "<< argv[0]
              <<" <inputfile> [-dense|-spr|-superlu[<nt>]|-samg|-petsc]\n"
              <<"       [-lag|-spec|-LR] [-1D|-2D] [-nGauss <n>]\n"
              <<"       [-hdf5] [-vtf <format> [-nviz <nviz>]"
              <<" [-nu <nu>] [-nv <nv>] [-nw <nw>]]\n       [-adap[<i>]]"
              <<" [-eig <iop> [-nev <nev>] [-ncv <ncv] [-shift <shf>] [-free]]"
              <<"\n       [-ignore <p1> <p2> ...] [-fixDup]"
              <<" [-checkRHS] [-check] [-dumpASC]\n";
    return 0;
  }

  if (args.adap)
    IFEM::getOptions().discretization = ASM::LRSpline;

  IFEM::cout <<"\nInput file: "<< infile;
  IFEM::getOptions().print(IFEM::cout);
  if (SIMbase::ignoreDirichlet)
    IFEM::cout <<"\nSpecified boundary conditions are ignored";
  if (fixDup)
    IFEM::cout <<"\nCo-located nodes will be merged";
  if (checkRHS && args.dim > 1)
    IFEM::cout <<"\nCheck that each patch has a right-hand coordinate system";
  if (!ignoredPatches.empty())
  {
    IFEM::cout <<"\nIgnored patches:";
    for (size_t i = 0; i < ignoredPatches.size(); i++)
      IFEM::cout <<" "<< ignoredPatches[i];
  }
  IFEM::cout << std::endl;

  utl::profiler->stop("Initialization");

  if (args.adap) {
    if (args.dim == 1)
      return runSimulator<SIM1D,SIMSolverAdap>(infile, checkRHS, ignoredPatches,
                                               fixDup, vizRHS, dumpASCII);
    else if (args.dim == 2)
      return runSimulator<SIM2D,SIMSolverAdap>(infile, checkRHS, ignoredPatches,
                                               fixDup, vizRHS, dumpASCII);
    else
      return runSimulator<SIM3D,SIMSolverAdap>(infile, checkRHS, ignoredPatches,
                                               fixDup, vizRHS, dumpASCII);
  } else {
    if (args.dim == 1)
      return runSimulator<SIM1D>(infile, checkRHS, ignoredPatches,
                                 fixDup, vizRHS, dumpASCII);
    else if (args.dim == 2)
      return runSimulator<SIM2D>(infile, checkRHS, ignoredPatches,
                                 fixDup, vizRHS, dumpASCII);
    else
      return runSimulator<SIM3D>(infile, checkRHS, ignoredPatches,
                                 fixDup, vizRHS, dumpASCII);
  }
}
