#include <iostream>
#include <vector>

#include <AMReX_MultiFab.H>
#include <AMReX_VisMF.H>
#include <AMReX_ParmParse.H>
#include <AMReX_PlotFileUtil.H>

#include "mechanism.H"
#include <GPU_misc.H>

#include <PelePhysics.H>

int
main(int argc, char* argv[])
{
  amrex::Initialize(argc, argv);

  static pele::physics::eos::EosParams<
         pele::physics::PhysicsType::eos_type>
         eos_parms;
  static pele::physics::TabFuncParams tabfunc_data;
  amrex::Print() << " Initialization of EOS (CPP)... \n";
#ifdef USE_MANIFOLD_EOS
  amrex::Print() << " Initialization of Table (CPP)... \n";
  tabfunc_data.initialize();
  eos_parms.allocate(tabfunc_data.device_tabfunc_data());
#else
  eos_parms.allocate();
#endif

  {
    amrex::ParmParse pp;
    int size = 128;
    pp.query("size",size);
    int niter = 100;
    pp.query("niter",niter);
    int do_plt = 0;
    pp.query("do_plt",do_plt);

    // Define geometry
    amrex::Array<int, AMREX_SPACEDIM> npts{AMREX_D_DECL(1, 1, 1)};

    for (int i = 0; i < AMREX_SPACEDIM; ++i) {
      npts[i] = size;
    }

    amrex::Box domain(
      amrex::IntVect(AMREX_D_DECL(0, 0, 0)),
      amrex::IntVect(AMREX_D_DECL(npts[0] - 1, npts[1] - 1, npts[2] - 1)));

    amrex::RealBox real_box(
      {AMREX_D_DECL(-1.0, -1.0, -1.0)}, {AMREX_D_DECL(1.0, 1.0, 1.0)});

    int coord = 0;

    amrex::Array<int, AMREX_SPACEDIM> is_periodic{AMREX_D_DECL(1, 1, 1)};

    amrex::Geometry geom(domain, real_box, coord, is_periodic);

    // Define BoxArray
    int max_size = 32;
    pp.query("max_size", max_size);
    amrex::BoxArray ba(domain);
    ba.maxSize(max_size);

    amrex::ParmParse ppa("amr");
    std::string pltfile("plt");
    ppa.query("plot_file", pltfile);

    amrex::DistributionMapping dm{ba};
    int num_grow = 0;

    // Data MFs
    amrex::MultiFab mass_frac(ba, dm, NUM_SPECIES, num_grow);
    amrex::MultiFab temperature(ba, dm, 1, num_grow);
    amrex::MultiFab density(ba, dm, 1, num_grow);
    amrex::MultiFab energy(ba, dm, 1, num_grow);

    const auto geomdata = geom.data();

    std::cout << "Initialize Data" << std::endl;
    {
      BL_PROFILE("PelePhysics::InitializeData()");
#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for (amrex::MFIter mfi(mass_frac, amrex::TilingIfNotGPU()); mfi.isValid();
         ++mfi) {

      const amrex::Box& bx = mfi.tilebox();
      auto const& Y_a = mass_frac.array(mfi);
      auto const& T_a = temperature.array(mfi);
      auto const& rho_a = density.array(mfi);
      auto const& e_a = energy.array(mfi);
      auto const* leosparm = eos_parms.device_eos_parm();
      amrex::ParallelFor(
        bx, [Y_a, T_a, rho_a, e_a,
             geomdata, leosparm] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
	      initialize_data(i, j, k, Y_a, T_a, rho_a, e_a, geomdata, leosparm);
        });
    }
    }

    amrex::MultiFab VarPlt(ba, dm, 4, num_grow);
    amrex::MultiFab cp(ba, dm, 1, num_grow);
    amrex::MultiFab wdot(ba, dm, NUM_SPECIES, num_grow);
    amrex::MultiFab Tout(ba, dm, 1, num_grow);
    Tout.setVal(250.0);

    std::cout << "Get Cp" << std::endl;
    {
      BL_PROFILE("PelePhysics::GetCp()");
      for (int iter = 0; iter < niter; iter++) {
#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
	for (amrex::MFIter mfi(mass_frac, amrex::TilingIfNotGPU()); mfi.isValid();
	     ++mfi) {
	  const amrex::Box& box = mfi.tilebox();
	  auto const& Y_a = mass_frac.const_array(mfi);
	  auto const& T_a = temperature.const_array(mfi);
	  auto const& cp_a = cp.array(mfi);
	  auto const* leosparm = eos_parms.device_eos_parm();
	  amrex::ParallelFor(
			     box, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
				    get_cp(i, j, k, Y_a, T_a, cp_a, leosparm);
				  });
	}
      }
    }

    std::cout << "Get T1" << std::endl;
    {
      BL_PROFILE("PelePhysics::GetT1()");
      for (int iter = 0; iter < niter; iter++) {
#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
	for (amrex::MFIter mfi(mass_frac, amrex::TilingIfNotGPU()); mfi.isValid();
	     ++mfi) {
	  const amrex::Box& box = mfi.tilebox();
	  auto const& Y_a = mass_frac.const_array(mfi);
	  auto const& e_a = energy.const_array(mfi);
	  auto const& T_a = Tout.array(mfi);
	  auto const* leosparm = eos_parms.device_eos_parm();
	  amrex::ParallelFor(
			     box, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
				    get_T_from_EY(i, j, k, Y_a, T_a, e_a, leosparm);
				  });
	}
      }
    }

    std::cout << "Get Wdot" << std::endl;
    {
      BL_PROFILE("PelePhysics::GetWdot()");
      for (int iter = 0; iter < niter; iter++) {
#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
	for (amrex::MFIter mfi(mass_frac, amrex::TilingIfNotGPU()); mfi.isValid();
	     ++mfi) {
	  const amrex::Box& box = mfi.tilebox();
	  auto const& Y_a = mass_frac.const_array(mfi);
	  auto const& T_a = temperature.const_array(mfi);
	  auto const& rho_a = density.const_array(mfi);
	  auto const& wdot_a = wdot.array(mfi);
	  auto const* leosparm = eos_parms.device_eos_parm();
	  amrex::ParallelFor(
			     box, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
				    get_wdot(i, j, k, Y_a, rho_a, T_a, wdot_a, leosparm);
				  });
	}
      }
    }

    std::cout << "Get T" << std::endl;
    {
      BL_PROFILE("PelePhysics::GetT()");
      for (int iter = 0; iter < niter; iter++) {
#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
	for (amrex::MFIter mfi(mass_frac, amrex::TilingIfNotGPU()); mfi.isValid();
	     ++mfi) {
	  const amrex::Box& box = mfi.tilebox();
	  auto const& Y_a = mass_frac.const_array(mfi);
	  auto const& e_a = energy.const_array(mfi);
	  auto const& T_a = Tout.array(mfi);
	  auto const* leosparm = eos_parms.device_eos_parm();
	  amrex::ParallelFor(
			     box, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
				    get_T_from_EY(i, j, k, Y_a, T_a, e_a, leosparm);
				  });
	}
      }
    }

    if (do_plt) {
      amrex::MultiFab::Copy(VarPlt, Tout, 0, 0, 1, num_grow);
      amrex::MultiFab::Copy(VarPlt, temperature, 0, 1, 1, num_grow);
      amrex::MultiFab::Copy(VarPlt, cp, 0, 2, 1, num_grow);
      amrex::MultiFab::Copy(VarPlt, wdot, 0, 3, 1, num_grow);

      std::string outfile = amrex::Concatenate(pltfile, 1);
      // TODO: add fct count to this output
      amrex::Vector<std::string> plt_VarsName;
      plt_VarsName.push_back("Tout");
      plt_VarsName.push_back("temperature");
      plt_VarsName.push_back("cp");
      plt_VarsName.push_back("wdot");
      amrex::WriteSingleLevelPlotfile(outfile, VarPlt, plt_VarsName, geom, 0.0, 0);
    }
  }

  eos_parms.deallocate();
#ifdef USE_MANIFOLD_EOS
  tabfunc_data.deallocate();
#endif
  amrex::Finalize();

  return 0;
}
