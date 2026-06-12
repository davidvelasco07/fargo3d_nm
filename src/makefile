# If you want your environment variables to have precedence
# you must execute the makefile with the flag "-e" (e.g.: make -e).
#------------------------OPTIONS FILE---------------------------------

SETUPSDIR=../setups
BINDIR=../bin
STDDIR=../std
SRCDIR=../src
ARCHDIR=../arch
ARCHNAME = source.tar.bz2
SCRIPTSDIR=../scripts
GIT_VERSION := $(shell git describe --abbrev=8 --dirty --always --tags 2>/dev/null)
ARCH_CREATE := $(shell mkdir -p $(ARCHDIR))
SHORTARCH := $(shell echo $(ARCHDIR) | sed -e 's|../||')

#------------------------OPTIONS FILE---------------------------------

BLACK=\033[0;30m
BLUE=\033[1;34m #  \e[m'
GREEN=\033[1;32m
CYAN=\033[1;36m
RED=\033[1;31m
PURPLE=\033[1;35m
BROWN=\033[0;33m
LIGHT_GRAY=\033[0;37m
END_COLOR=\033[0m

ifdef FARGO_COLOR
	FARGO_OPT += -DFARGO_COLOR
endif

ifeq (${BIGMEM},1)
	FARGO_OPT += -DBIGMEM
endif
ifeq (${RESCALE}, 1)
	FARGO_OPT += -DRESCALE
endif
ifeq (${UNITS}, MKS)
	FARGO_OPT += -DMKS
endif
ifeq (${UNITS}, CGS)
	FARGO_OPT += -DCGS
endif
ifeq (${NOGHOSTSX}, 1)
	FARGO_OPT += -DNOGHOSTSX
endif
ifeq (${MPICUDA}, 1)
	FARGO_OPT += -DMPICUDA -DPARALLEL
endif
ifeq (${PROFILING}, 1)
	FARGO_OPT += -DPROFILING
endif
ifeq (${PARALLEL}, 1)
	FARGO_OPT += -DPARALLEL
endif
ifeq (${DEBUG}, 1)
	FARGO_OPT += -DDEBUG
	FARGO_DEBUG = 1
endif
ifeq (${FULLDEBUG}, 1)
	FARGO_OPT += -DDEBUG
	FARGO_OPT += -DWRITEGHOSTS
	FARGO_DEBUG = 1
endif
ifeq (${BLOCKS}, 1)
	FARGO_OPT += -DNOOUTPUTS
endif
ifeq (${GPU}, 1)
	FARGO_OPT += -DGPU
endif
ifdef SETUP
	FARGO_OPT_FILE = ${SETUPSDIR}/${SETUP}/${SETUP}.opt
	VARFILE = ${SETUPSDIR}/${SETUP}/${SETUP}.var
	PARFILE = ${SETUPSDIR}/${SETUP}/${SETUP}.par
endif
ifeq (${NOFLUSH}, 1)
	FARGO_OPT += -DNOFLUSH
endif
ifndef FARGO_OPT_FILE
FARGO_OPT_FILE = ${SETUPSDIR}/fargo/fargo.opt
VARFILE = ${SETUPSDIR}/fargo/fargo.var
endif
ifndef FARGO_DISPLAY
FARGO_DISPLAY = NONE
endif
#---------------------------------------------------------------------
include ${FARGO_OPT_FILE} #The name must not contain digits
#--------------------FARGO_VARIABLES----------------------------------
ifndef FARGO_ARCH
FARGO_ARCH  = LINUX
endif
ifeq (-DGPU,${findstring -DGPU,${FARGO_OPT}})
	FARGO_GPU   = 1
else
	FARGO_GPU   = 0
endif
ifeq (-DDEBUG,${findstring -DDEBUG,${FARGO_OPT}})
	FARGO_DEBUG = 1
endif
#----------------------SYSTEM VARIABLES-------------------------------
ifndef CUDA
CUDA           = /usr/local/cuda
endif
ifndef PYTHON_INCLUDE
PYTHON_INCLUDE = /usr/include/python2.7
endif
#---------------------------------------------------------------------
#LINUX PLATFORM (GENERIC, default)
#FARGO_ARCH must be set to LINUX
CC_LINUX      = gcc 
SEQOPT_LINUX  = -O3 -ffast-math
PARAOPT_LINUX = ${SEQOPT_LINUX}
PARACC_LINUX  = mpicc
LIBS_LINUX    = -lm
INC_LINUX     = 
NVCC_LINUX    = nvcc
CUDAOPT_LINUX = -O3 -arch=sm_20 -w
PARAINC_LINUX = 
PARALIB_LINUX =
CUDAINC_LINUX = -I${CUDA}/include/
CUDALIB_LINUX = -L${CUDA}/lib64/  -lcudart
CXX_LINUX     = g++ -stdlib=libstdc++
CXXPARA_LINUX = mpic++ -stdlib=libstdc++
#---------------------------------------------------------------------
#PRINCETON DELLA CLUSTER (CUDA-aware OpenMPI, A100 GPUs)
#FARGO_ARCH must be set to DELLA
#Load modules first, e.g.:
#  module load cudatoolkit/12.6 openmpi/cuda-12.6/gcc/4.1.6
#and export CUDA=$CUDA_HOME before building.
CC_DELLA      = gcc
# -fcommon restores pre-GCC-10 tentative-definition merging, which this
# legacy code relies on (globals declared without 'extern' in global.h).
SEQOPT_DELLA  = -O3 -ffast-math -fcommon
PARAOPT_DELLA = ${SEQOPT_DELLA}
PARACC_DELLA  = mpicc
LIBS_DELLA    = -lm
ENVRANK_DELLA = OMPI_COMM_WORLD_LOCAL_RANK
INC_DELLA     =
NVCC_DELLA    = nvcc
CUDAOPT_DELLA = -O3 -arch=sm_80 -w
PARAINC_DELLA =
PARALIB_DELLA =
CUDAINC_DELLA = -I${CUDA}/include/
CUDALIB_DELLA = -L${CUDA}/lib64/ -lcudart
CXX_DELLA     = g++
CXXPARA_DELLA = mpic++
#---------------------------------------------------------------------                                                                                                                                                                                                    
#PIZ DAINT PLATFORM AT CSCS                                                                                                                                                                                                                                                      
#FARGO_ARCH must be set to DAINT                                                                                                                                                                                                                                  
CC_DAINT      = cc
SEQOPT_DAINT  = -O3
PARAOPT_DAINT = ${SEQOPT_DAINT}
PARACC_DAINT  = cc
LIBS_DAINT    = -lm
ENVRANK_DAINT = OMPI_COMM_WORLD_LOCAL_RANK 
INC_DAINT     = -I../include
NVCC_DAINT    = nvcc
CUDAOPT_DAINT = -O3 -arch=sm_60 -I$(CRAY_MPICH2_DIR)/include
PARAINC_DAINT =
PARALIB_DAINT =
CUDAINC_DAINT = 
CUDALIB_DAINT = -lcuda -lcudart
CXX_DAINT     = CC
CXXPARA_DAINT = CC
#---------------------------------------------------------------------
#S3IT PLATFORM AT UZH
#FARGO_ARCH must be set to VOLTA
CC_VOLTA      = gcc
SEQOPT_VOLTA  = -O3
PARAOPT_VOLTA = ${SEQOPT_VOLTA}
PARACC_VOLTA  = mpicc
LIBS_VOLTA    = -lm
ENVRANK_VOLTA = OMPI_COMM_WORLD_LOCAL_RANK
INC_VOLTA     = -I../include
NVCC_VOLTA    = nvcc
CUDAOPT_VOLTA = -O3 -arch=sm_37 --cudart=shared
PARAINC_VOLTA =
PARALIB_VOLTA =
CUDAINC_VOLTA = -I${CUDA_HOME}/include/
CUDALIB_VOLTA = -L${CUDA_HOME}/lib64/  -lcudart -lrt -ldl
CXX_VOLTA     = g++
CXXPARA_VOLTA = mpic++
#---------------------------------------------------------------------
#KEPLER PLATFORM WITH OPENMPI-1.7.4 AT UNAM
#FARGO_ARCH must be set to KEPLEROPENMPI
CC_KEPLEROPENMPI      = gcc 
SEQOPT_KEPLEROPENMPI  = -O3 -ffast-math
PARAOPT_KEPLEROPENMPI = ${SEQOPT_KEPLEROPENMPI}
PARACC_KEPLEROPENMPI  = mpicc
LIBS_KEPLEROPENMPI    = -lm
ENVRANK_KEPLEROPENMPI = OMPI_COMM_WORLD_LOCAL_RANK
INC_KEPLEROPENMPI     = 
NVCC_KEPLEROPENMPI    = nvcc
CUDAOPT_KEPLEROPENMPI = -O3 -arch=sm_35
PARAINC_KEPLEROPENMPI = 
PARALIB_KEPLEROPENMPI =
CUDAINC_KEPLEROPENMPI = -I${CUDA}/include/
CUDALIB_KEPLEROPENMPI = -L${CUDA}/lib64/  -lcudart
CXX_KEPLEROPENMPI     = g++ -stdlib=libstdc++
CXXPARA_KEPLEROPENMPI = mpic++ -stdlib=libstdc++
#---------------------------------------------------------------------
#KEPLER PLATFORM WITH MVAPICH2-2.0b AT UNAM
#FARGO_ARCH must be set to KEPLERMVAPICH
CC_KEPLERMVAPICH      = gcc 
SEQOPT_KEPLERMVAPICH  = -O3 -ffast-math
PARAOPT_KEPLERMVAPICH = ${SEQOPT_KEPLERMVAPICH}
PARACC_KEPLERMVAPICH  = mpicc
LIBS_KEPLERMVAPICH    = -lm
ENVRANK_KEPLERMVAPICH = MV2_COMM_WORLD_LOCAL_RANK
INC_KEPLERMVAPICH     = 
NVCC_KEPLERMVAPICH    = nvcc
CUDAOPT_KEPLERMVAPICH = -O3 -arch=sm_35
PARAINC_KEPLERMVAPICH = 
PARALIB_KEPLERMVAPICH = 
CUDAINC_KEPLERMVAPICH = -I${CUDA}/include/
CUDALIB_KEPLERMVAPICH = -L${CUDA}/lib64/ -lcudart -L/share/apps/mvapich2/lib
CXX_KEPLERMVAPICH     = g++ -stdlib=libstdc++
CXXPARA_KEPLERMVAPICH = mpic++ -stdlib=libstdc++
#---------------------------------------------------------------------
#TESLA PLATFORM WITH OPENMPI-1.7.2 AT UNAM
#FARGO_ARCH must be set to TESLAOPENMPI
CC_TESLAOPENMPI      = gcc 
SEQOPT_TESLAOPENMPI  = -O3 -ffast-math
PARAOPT_TESLAOPENMPI = ${SEQOPT_TESLAOPENMPI}
PARACC_TESLAOPENMPI  = mpicc
LIBS_TESLAOPENMPI    = -lm
ENVRANK_TESLAOPENMPI = OMPI_COMM_WORLD_LOCAL_RANK
INC_TESLAOPENMPI     = 
NVCC_TESLAOPENMPI    = nvcc
CUDAOPT_TESLAOPENMPI = -O3 -arch=sm_20
PARAINC_TESLAOPENMPI = 
PARALIB_TESLAOPENMPI =
CUDAINC_TESLAOPENMPI = -I${CUDA}/include/
CUDALIB_TESLAOPENMPI = -L${CUDA}/lib64/  -lcudart
CXX_TESLAOPENMPI     = g++
CXXPARA_TESLAOPENMPI = mpic++
#---------------------------------------------------------------------
#TESLA PLATFORM WITH MVAPICH2-1.9 AT UNAM
#FARGO_ARCH must be set to TESLAMVAPICH
CC_TESLAMVAPICH      = gcc 
SEQOPT_TESLAMVAPICH  = -O3 -ffast-math
PARAOPT_TESLAMVAPICH = ${SEQOPT_TESLAMVAPICH}
PARACC_TESLAMVAPICH  = mpicc
LIBS_TESLAMVAPICH    = -lm
ENVRANK_TESLAMVAPICH = MV2_COMM_WORLD_LOCAL_RANK
INC_TESLAMVAPICH     = 
NVCC_TESLAMVAPICH    = nvcc
CUDAOPT_TESLAMVAPICH = -O3 -arch=sm_20
PARAINC_TESLAMVAPICH = 
PARALIB_TESLAMVAPICH = -lmpich -lOpenCL
CUDAINC_TESLAMVAPICH = -I${CUDA}/include/
CUDALIB_TESLAMVAPICH = -lcudart
CXX_TESLAMVAPICH     = g++ -stdlib=libstdc++
CXXPARA_TESLAMVAPICH = mpic++ -stdlib=libstdc++
#---------------------------------------------------------------------
#TESLA C2070 PLATFORM WITH MVAPICH2-1.9 AT CORDOBA
#FARGO_ARCH must be set to MENDIETA
CC_MENDIETA      = gcc 
SEQOPT_MENDIETA  = -O3 -ffast-math
PARAOPT_MENDIETA = ${SEQOPT_MENDIETA}
PARACC_MENDIETA  = mpicc
LIBS_MENDIETA    = -lm
ENVRANK_MENDIETA = MV2_COMM_WORLD_RANK
INC_MENDIETA     = 
NVCC_MENDIETA    = nvcc
CUDAOPT_MENDIETA = -O3 -arch=sm_20
PARAINC_MENDIETA = 
PARALIB_MENDIETA = -lmpich -lOpenCL
CUDAINC_MENDIETA = -I${CUDA}/include/
CUDALIB_MENDIETA = -lcudart
CXX_MENDIETA     = g++ -stdlib=libstdc++
CXXPARA_MENDIETA = mpic++ -stdlib=libstdc++
#---------------------------------------------------------------------
#QUADRO PLATFORM (IATE COMPUTER)
#FARGO_ARCH must be set to QUADRO
CC_QUADRO      = gcc
SEQOPT_QUADRO  = -O3 #-ffast-math 
PARAOPT_QUADRO = ${SEQOPT_QUADRO}
PARACC_QUADRO  = mpicc
LIBS_QUADRO    = -lm
INC_QUADRO     = 
NVCC_QUADRO    = ${CUDA}/bin/nvcc -arch=sm_20
CUDAOPT_QUADRO = 
PARAINC_QUADRO = -I/usr/include/openmpi-x86_64/
PARALIB_QUADRO = -L/usr/lib64/openmpi/lib/
CUDALIB_QUADRO = -L${CUDA}/lib64/ -lcudart
CUDAINC_QUADRO = -I${CUDA}/include/
CXX_QUADRO     = g++ -stdlib=libstdc++
CXXPARA_QUADRO = mpic++ -stdlib=libstdc++
#---------------------------------------------------------------------
#MAC PLATFORM (FREDERIC'S COMPUTER) - Works with OSX 10.10 (Yosemite)
#and CUDA release 7.5, V7.5.19 + GPU Driver Version: 10.4.2 310.41.35f01
#FARGO_ARCH must be set to MacIntel
CC_MacIntel      = gcc
SEQOPT_MacIntel  = -O3 -arch x86_64 -Wmissing-prototypes #-Wall -Wextra
PARAOPT_MacIntel = ${SEQOPT_MacIntel}
PARACC_MacIntel  = mpicc
ENVRANK_MacIntel = OMPI_COMM_WORLD_RANK
LIBS_MacIntel    = 
INC_MacIntel     = 
NVCC_MacIntel    = ${CUDA}/bin/nvcc
CUDAOPT_MacIntel = -arch=sm_30 -m64 -O3 --compiler-options '-Wno-return-type-c-linkage' -Xcudafe --diag_suppress=set_but_not_used -Xcudafe --diag_suppress=declared_but_not_referenced
PARAINC_MacIntel = 
PARALIB_MacIntel = 
CUDALIB_MacIntel = -L${CUDA}/lib/ -lcudart
CUDAINC_MacIntel = -I${CUDA}/include/
CXX_MacIntel     = g++ -stdlib=libstdc++
CXXPARA_MacIntel = mpic++ -stdlib=libstdc++
#--------------Architecture independent names-------------------------
CC	 = ${CC_${FARGO_ARCH}}
SEQOPT   = ${SEQOPT_${FARGO_ARCH}}
PARAOPT  = ${PARAOPT_${FARGO_ARCH}}
PARACC   = ${PARACC_${FARGO_ARCH}}
LIBS     = ${LIBS_${FARGO_ARCH}}
NVCC     = ${NVCC_${FARGO_ARCH}}
ENVRANK  = ${ENVRANK_${FARGO_ARCH}}
CUDAOPT  = ${CUDAOPT_${FARGO_ARCH}}
PARAINC  = ${PARAINC_${FARGO_ARCH}}
PARALIB  = ${PARALIB_${FARGO_ARCH}}
CUDAINC  = ${CUDAINC_${FARGO_ARCH}}
CUDALIB  = ${CUDALIB_${FARGO_ARCH}}
CXX      = ${CXX_${FARGO_ARCH}}
CXXPARA  = ${CXXPARA_${FARGO_ARCH}}
INCLUDE  = ${INC_${FARGO_ARCH}} ${CUDAINC_${FARGO_ARCH}} \
	  -I${SRCDIR} -I${BINDIR}
#---------------------------------------------------------------------

MAKEFILE= makefile

MAKEFILE += ${SCRIPTSDIR}/make.py

GLOBAL   = ${MAKEFILE} ${SCRIPTSDIR}/c2cuda.py ${FARGO_OPT_FILE}

GLOBAL += fargo3d.h define.h types_def.h prototypes.h fondam.h structs.h

BOUNDARIES =  boundparser.py boundaries.txt 			    \
	      centering.txt  				    	    \
	      $(wildcard ${SETUPSDIR}/${SETUP}/$(SETUP).bound) 	    \
	      boundary_template.c

MPIDUMMY = mpi_dummy.o #dummy file, for sequential runs.
EXENAME  = ../fargo3d  #Name and location of the executable file

VPATH=${SETUPSDIR}/${SETUP}:${STDDIR}:${SCRIPTSDIR}:${SRCDIR}

#---------------#Sequential and parallel objects----------------------
MAINOBJ  = LowTasks.o psys.o main.o param.o var.o usage.o   	\
	   split.o comm.o CompPresIso.o         		\
	   boundary.o checknans.o 				\
	   vanleer_x_a.o vanleer_x_b.o 				\
           vanleer_y_a.o vanleer_y_b.o				\
	   vanleer_z_a.o vanleer_z_b.o 				\
	   divrho.o transport.o 				\
	   substep3.o CompPresAd.o CompPresPoly.o		\
	   substep2_a.o substep2_b.o condinit.o     		\
	   output.o output2d.o output_vtk.o 			\
	   output_vtk_merging.o			                \
	   algogas.o planets.o runge_kutta.o	        	\
	   stockholm.o mhd.o potential.o	           	\
	   timeinfo.o allocate.o resistivity.o		      	\
	   defout.o fresh.o  light_global_dev.o          	\
	   fargo_mhd.o substep1_x.o substep1_y.o substep1_z.o	\
	   compfields.o advect_shift.o rescale.o		\
	   momenta_x.o momenta_y.o momenta_z.o			\
	   reduction_sum.o reduction_min.o			\
	   update_a_x.o update_a_y.o update_a_z.o 		\
	   update_dens_a_x.o update_dens_a_y.o			\
	   update_dens_a_z.o multrho.o					\
	   update_b_x.o update_b_y.o update_b_z.o 		\
	   newvel_x.o newvel_y.o newvel_z.o 			\
	   compute_residual.o change_frame.o compute_vmed.o	\
           correct_vtheta.o cfl.o cfl_b.o			\
	   compute_force.o compute_slopes.o compute_star.o	\
	   compute_emf.o update_magnetic.o lorentz_force.o	\
	   fargo_ppa_a.o fargo_ppa_b.o fargo_ppa_c.o            \
	   fargo_ppa_d.o fargo_ppa_d_2d.o fargo_ppa_steep.o	\
	   vanleerx_ppa.o dh_boundary.o cpvel.o resist.o	\
	   integrate_emf.o monitor.o Compute_Mean.o 		\
	   viscosity.o init_stockholm.o visctensor_cart.o	\
	   addviscosity_cart.o addviscosity_cyl.o		\
	   addviscosity_sph.o visctensor_sph.o			\
	   visctensor_cyl.o stretch.o 				\
	   xmin_bound.o xmax_bound.o				\
	   ymin_bound.o ymax_bound.o				\
           zmin_bound.o zmax_bound.o				\
	   mon_dens.o mon_momx.o mon_momy.o mon_momz.o          \
	   mon_torque.o comm_device.o checkmute_y.o		\
	   checkmute_z.o mon_reynolds.o mon_maxwell.o           \
	   mon_bxflux.o mon_torque_excl.o                       \
	   postrestarthook.o change_arch.o select_device.o	\
	   redefoptions.o setuphook1.o				\
	   RT_main.o RT_TemperatureToEnergy.o RT_CompTemp.o     \
	   RT_MatrixIterate.o RT_RadEnergyToTemperature.o       \
	   RT_MatrixElements.o RT_ComputeEtas.o PlanetHeating.o \
	   RT_ComputeDiffCoef.o RT_CompOpacity.o                \
	   RT_boundaries.o RT_SolveMatrix.o RT_NormRHS.o	\
	   RT_heating.o RT_GuessErad.o RT_DiffErad.o		\
           fill_ghostsx.o                                       \
	   fill_ghostsx_erad.o                                  \
	   remap.o get_new_borders.o remap_copy.o               \
	   remap_geometry_stock.o J_refine_desc.o J_refine_field.o\
	   J_gridfile.o J_gridabs.o J_log.o J_gridbuild.o       \
	   J_low.o J2F3D.o J_split.o J_create.o J_create_fluid.o\
	   J_desc.o J_comm_init.o J_comm_strip.o J_overlap.o    \
	   J_comm_alloc.o J_comm_check.o J_gridread.o           \
	   J_Output.o J_comm_exec.o J_comm_fill.o               \
	   J_faceflux.o J_correct.o J_resetfluxes.o             \
	   J_subcycling.o J_iter.o J_comm_device.o              \
	   J_merge_field.o J_merge_desc.o J_readfield.o		\
	   summary.o substep4_a.o substep4_b.o	mon_force.o	\
	   LumExponentialSrc.o LumPoissonSrc.o			\
	   LumPointSrc.o LumSphericalSrc.o Luminous_sphere.o	

ifeq (-DIRRADIATION,${findstring -DIRRADIATION,${FARGO_OPT}})
MAINOBJ += RT_IRR_Intercpuscan.o RT_IRR_LocalScan.o		\
	   RT_IRR_OffsetPatches.o RT_IRR_RadialScan.o		\
	   RT_IRR_TauStellar.o RT_IRR_getbgslice.o              \
	   RT_IRR_setbgslice.o RT_IRR_StellarHeating.o
endif

#-----------------------#Cuda objects---------------------------------
GPU_OBJ  =  reduction_sum_device.o reduction_min_device.o	\
	   	xmin_bound_gpu.o xmax_bound_gpu.o		\
		ymin_bound_gpu.o ymax_bound_gpu.o		\
	        zmin_bound_gpu.o zmax_bound_gpu.o               \
	   	lowtasks_gpu.o 

ifeq (-DIRRADIATION,${findstring -DIRRADIATION,${FARGO_OPT}})
GPU_OBJ  += RT_IRR_localscan_gpu.o
endif

#<CUDA OBJECTS> #Do not remove!!!
GPU_OBJBLOCKS =  CompPresIso_gpu.o CompPresAd_gpu.o 		\
	   compute_slopes_gpu.o compute_star_gpu.o		\
	   compute_emf_gpu.o update_magnetic_gpu.o		\
	   substep1_x_gpu.o substep1_y_gpu.o substep1_z_gpu.o 	\
           substep2_a_gpu.o substep2_b_gpu.o substep3_gpu.o 	\
	   divrho_gpu.o stockholm_gpu.o 			\
	   vanleer_x_a_gpu.o vanleer_x_b_gpu.o 			\
	   vanleer_y_a_gpu.o vanleer_y_b_gpu.o 			\
	   vanleer_z_a_gpu.o vanleer_z_b_gpu.o 			\
	   momenta_x_gpu.o momenta_y_gpu.o momenta_z_gpu.o	\
	   update_a_x_gpu.o update_a_y_gpu.o update_a_z_gpu.o   \
	   update_dens_a_x_gpu.o update_dens_a_y_gpu.o		\
	   update_dens_a_z_gpu.o multrho_gpu.o				\
	   update_b_x_gpu.o update_b_y_gpu.o update_b_z_gpu.o   \
	   newvel_x_gpu.o newvel_y_gpu.o newvel_z_gpu.o 	\
	   advect_shift_gpu.o compute_residual_gpu.o 		\
	   change_frame_gpu.o potential_gpu.o 			\
	   correct_vtheta_gpu.o cfl_gpu.o compute_force_gpu.o	\
	   lorentz_force_gpu.o CompPresPoly_gpu.o		\
	   fargo_ppa_a_gpu.o fargo_ppa_b_gpu.o			\
	   fargo_ppa_c_gpu.o fargo_ppa_d_gpu.o			\
	   fargo_ppa_d_2d_gpu.o resist_gpu.o 			\
	   integrate_emf_gpu.o visctensor_cart_gpu.o 		\
	   addviscosity_cart_gpu.o addviscosity_cyl_gpu.o	\
	   visctensor_cyl_gpu.o addviscosity_sph_gpu.o		\
	   mon_dens_gpu.o mon_momx_gpu.o mon_momy_gpu.o         \
	   mon_momz_gpu.o mon_torque_gpu.o checkmute_y_gpu.o	\
	   checkmute_z_gpu.o mon_reynolds_gpu.o                 \
	   mon_maxwell_gpu.o visctensor_sph_gpu.o               \
	   mon_bxflux_gpu.o mon_torque_excl_gpu.o                                    \
	   fargo_ppa_steep_gpu.o setuphook1_gpu.o		\
	   RT_TemperatureToEnergy_gpu.o RT_CompTemp_gpu.o       \
	   RT_MatrixIterate_gpu.o RT_NormRHS_gpu.o              \
	   RT_MatrixElements_gpu.o RT_ComputeEtas_gpu.o         \
	   RT_ComputeDiffCoef_gpu.o RT_CompOpacity_gpu.o        \
	   RT_boundaries_gpu.o RT_RadEnergyToTemperature_gpu.o  \
	   RT_heating_gpu.o RT_GuessErad_gpu.o			\
	   RT_DiffErad_gpu.o PlanetHeating_gpu.o                \
           fill_ghostsx_gpu.o fill_ghostsx_erad_gpu.o           \
	   remap_gpu.o remap_copy_gpu.o J_comm_kernels.o	\
	   J_comm_exec_kernels.o J_correct_kernels.o 		\
	   J_faceflux_kernels.o	substep4_a_gpu.o substep4_b_gpu.o\
	   LumExponentialSrc_gpu.o LumPoissonSrc_gpu.o		 \
	   mon_force_gpu.o LumPointSrc_gpu.o LumSphericalSrc_gpu.o 

ifeq (-DIRRADIATION,${findstring -DIRRADIATION,${FARGO_OPT}})
GPU_OBJBLOCKS += RT_IRR_OffsetPatches_gpu.o		        \
	   RT_IRR_TauStellar_gpu.o RT_IRR_getbgslice_gpu.o	\
	   RT_IRR_setbgslice_gpu.o RT_IRR_StellarHeating_gpu.o
endif
#</CUDA OBJECTS>
#---------------------------------------------------------------------

MAINOBJ += ${CPU_OBJECTS}
GPUOBJ += ${GPU_OBJ} ${GPU_OBJBLOCKS}

-include ${SETUPSDIR}/${SETUP}/${SETUP}.objects

ifeq (-DVISCOUSDRIFT,${findstring -DVISCOUSDRIFT,${FARGO_OPT}})
	LIBS += -lgsl -lgslcblas
endif

ifeq (-DPARALLEL,${findstring -DPARALLEL,${FARGO_OPT}})
	OPTIONS += ${PARAOPT}
	COMPILER = ${PARACC}
	LIBS += ${PARALIB}
else
	OPTIONS  += ${SEQOPT}
	MAINOBJ += ${MPIDUMMY}
	COMPILER = ${CC}
endif
	LINKER = ${COMPILER}

#----------------------------------------------------------------------
# If Stellar irradiation is included, the library cudpp is required to
# perform scans on the GPU. In this case, one needs the C++ compiler
# to link the object files rather than the C compiler, in order to
# properly include the template library 'thrust', needed by cudpp.
#----------------------------------------------------------------------

ifeq (${FARGO_GPU}, 1)
ifeq (-DIRRADIATION,${findstring -DIRRADIATION,${FARGO_OPT}})
	LINKER = ${CXX}
	LIBS += -lcudpp
ifeq (-DPARALLEL,${findstring -DPARALLEL,${FARGO_OPT}})
	LINKER = ${CXXPARA} 
endif
endif
endif

OBJECTS = ${MAINOBJ}

ifeq (${FARGO_GPU}, 1)
	OBJECTS += ${GPUOBJ}
	LIBS += ${CUDALIB}
endif

#---------------------------------------------------------------------
#Adding includes and objects if MATPLOTLIB is activated.
#These lines prevent compatibility problems 
#if matplotlib is not installed
ifeq (MATPLOTLIB,${findstring MATPLOTLIB,${FARGO_DISPLAY}})
	FARGO_OPT += -DMATPLOTLIB
	INCLUDE   += -I${PYTHON_INCLUDE}
	OBJECTS   +=  matplotlib.o
	LIBS      += -lpython2.7
endif
# Debug mode if FARGO_DEBUG = 1, Note that this mode is only 
# meant to execute gdb (or lldb, etc.) or valgrind.
ifeq (${FARGO_DEBUG}, 1)
OPTIONS = -g -fcommon
endif

ifeq (${MPICUDA}, 1)
	FARGO_OPT += -DENVRANK='${ENVRANK}'
endif

all:
	@echo "You must be in the upper level directory to invoke 'make'"
	@echo "Do 'cd ..' then 'make' to proceed."

allp: ${EXENAME}
	@echo ""
ifdef FARGO_COLOR
	@echo "${BLUE}"
endif
	@echo "          FARGO3D SUMMARY:           "
	@echo "          ===============            "
	@echo ""
ifeq (-DPARALLEL,${findstring -DPARALLEL,${FARGO_OPT}})
	@echo 'This built is PARALLEL (MPI). Use "make seq" to change that'
ifeq (${FARGO_GPU}, 1)
	@echo ""
	@echo "This built can be launched on"
	@echo "a cluster of GPUs."
endif
else
	@echo 'This built is SEQUENTIAL. Use "make para" to change that'
ifeq (${FARGO_GPU}, 1)
	@echo ""
	@echo "This built can be launched on"
	@echo "a CPU with a GPU card (1 GPU only)."
endif
endif
ifeq (MATPLOTLIB,${findstring MATPLOTLIB,${FARGO_DISPLAY}})	
	@echo ""
	@echo "This built has a graphical output,"
	@echo "which uses Python's matplotlib library."
	@echo 'Use "make noview" to change that.'
endif
	@echo ""
	@echo ""
	@echo "SETUP:      '"${SETUP}"'         "
	@echo '(Use "make SETUP=[valid_setup_string]" to change set up)'
	@echo '(Use "make list" to see the list of setups implemented)'
	@echo '(Use "make info" to see the current sticky build options)'
	@echo ""
ifeq (${PROFILING}, 1)
	@echo ""
	@echo "==================================="
	@echo "|   PROFILING MODE ACTIVATED!     |"
	@echo "==================================="
	@echo '(Use "make noprof" to change that)'
	@echo ""
endif
ifeq (${MPICUDA}, 1)
	@echo ""
	@echo "==================================="
	@echo "|   CUDA AWARE MPI!               |"
	@echo "==================================="
	@echo "Direct GPU-GPU MPI communications"
	@echo '(Use "make nompicuda" to change that)'
	@echo ""
endif
ifeq (${DEBUG}, 1)
	@echo ""
	@echo "==================================="
	@echo "|       DEBUG MODE ACTIVATED!     |"
	@echo "==================================="
	@echo '(Use "make nodebug" to change that)'
	@echo ""
endif
ifeq (${FULLDEBUG}, 1)
	@echo ""
	@echo "==================================="
	@echo "|  FULL DEBUG MODE ACTIVATED!     |"
	@echo "==================================="
	@echo "In FULL DEBUG mode buffer zones (ghosts) are included in output files"
	@echo '(Use "make nofulldebug" to change that)'
	@echo ""
endif
ifdef FARGO_COLOR
	@echo "${END_COLOR}"
endif

#---------------------------------------------------------------------
#------------------------MAIN RULE------------------------------------

${EXENAME} : ${OBJECTS}
	@echo ""
	@echo "All objects are OK. Linking stage"
# You may remove the '@' sign in the instructions below if you need a verbose output
	@${LINKER} ${OBJECTS} -o ${EXENAME} ${LIBS}
	@/bin/cp -f ${SRCDIR}/*.h $(ARCHDIR)
	@/bin/cp -f ${BINDIR}/*.h $(ARCHDIR)
	@cd $(ARCHDIR)/..;/bin/rm -f temparch; tar jcf temparch $(SHORTARCH) 2> /dev/null
	@mv $(ARCHDIR)/../temparch $(ARCHDIR)/$(ARCHNAME)
	@cat ${SETUPSDIR}/${SETUP}/${SETUP}.bound >> $(EXENAME)
	@wc -c ${SETUPSDIR}/${SETUP}/${SETUP}.bound | sed -e "s|${SETUPSDIR}/${SETUP}/${SETUP}.bound.*$$||" | xargs printf "%012d" >> ${EXENAME}
	@cat $(ARCHDIR)/$(ARCHNAME) >> $(EXENAME)
	@wc -c $(ARCHDIR)/$(ARCHNAME) | sed -e "s|$(ARCHDIR).*$$||" | xargs printf "%012d" >> ${EXENAME}
	@cat $(STDDIR)/.lastflags >> ${EXENAME}
	@wc -c $(STDDIR)/.lastflags | sed -e "s|$(STDDIR).*$$||" | xargs printf "%012d" >> ${EXENAME}
ifdef FARGO_COLOR
	@echo "${GREEN}"
	@echo "The executable file has been created."
	@echo "${END_COLOR}"
endif
#---------------------------------------------------------------------

var.c: ${PARFILE} stdpar.par #$(VARFILE) default_variables.dat
	@echo ""
	@echo "Generating var.c file"
	@python ${SCRIPTSDIR}/par.py $(PARFILE) ${STDDIR}/stdpar.par

${SETUP}.units:
	@touch ${SETUPSDIR}/${SETUP}/${SETUP}.units

rescale.c: ${PARFILE} ${SCRIPTSDIR}/unitparser.py ${STDDIR}/standard.units ${SETUP}.units
	@echo ""
	@echo "Generating rescale.c file"
	@python ${SCRIPTSDIR}/unitparser.py ${SETUP}

.PHONY: clean mrproper force
gitversion: force
	@echo '$(GIT_VERSION)' | cmp -s - $@ || echo '$(GIT_VERSION)' > $@
clean:
	@-rm -f *.o *~ *.s *.il global_ex.h param.h var.c xmin_bound.c xmax_bound.c ymin_bound.c ymax_bound.c zmin_bound.c zmax_bound.c rescale.c
	@echo
	@echo "Everything has been cleaned"
	@echo

main.o: gitversion

%_gpu.o : %.cu ${GLOBAL} #Rule: test_gpu.o is created from test.cu
	@echo NVCC $< "==>" $*_gpu.o
	@/bin/cp -f $< $(ARCHDIR)
	@${NVCC} $< -c ${CUDAOPT} ${CUDAINC} ${INCLUDE} -o $*_gpu.o ${FARGO_OPT}

%.cu: %.c
	@echo C2CUDA $< "==>" $*.cu
	@python ${SCRIPTSDIR}/c2cuda.py -i $< -o $*.cu -s ${SETUP}	

reduction_sum_device.o: reduction_sum_device.cu
	@echo NVCC $< "==>" $*.o
	@/bin/cp -f $< $(ARCHDIR)
	@${NVCC} $^ -c ${CUDAOPT} ${CUDAINC} ${INCLUDE}	${FARGO_OPT}

reduction_min_device.o: reduction_min_device.cu
	@echo NVCC $< "==>" $*.o
	@/bin/cp -f $< $(ARCHDIR)
	@${NVCC} $^ -c ${CUDAOPT} ${CUDAINC} ${INCLUDE}	${FARGO_OPT}

J_comm_kernels.o: J_comm_kernels.cu
	@echo NVCC $< "==>" $*.o
	@/bin/cp -f $< $(ARCHDIR)
	@${NVCC} $^ -c ${CUDAOPT} ${CUDAINC} ${INCLUDE} ${FARGO_OPT} 

J_comm_exec_kernels.o: J_comm_exec_kernels.cu
	@echo NVCC $< "==>" $*.o
	@/bin/cp -f $< $(ARCHDIR)
	@${NVCC} $^ -c ${CUDAOPT} ${CUDAINC} ${INCLUDE} ${FARGO_OPT} 

J_correct_kernels.o: J_correct_kernels.cu
	@echo NVCC $< "==>" $*.o
	@/bin/cp -f $< $(ARCHDIR)
	@${NVCC} $^ -c ${CUDAOPT} ${CUDAINC} ${INCLUDE}	${FARGO_OPT}

J_faceflux_kernels.o: J_faceflux_kernels.cu
	@echo NVCC $< "==>" $*.o
	@/bin/cp -f $< $(ARCHDIR)
	@${NVCC} $^ -c ${CUDAOPT} ${CUDAINC} ${INCLUDE}	${FARGO_OPT}

%_bound.o: xmin_bound.c ${GLOBAL}
	@/bin/cp -f $*_bound.c $(ARCHDIR)
	@/bin/cp -f $< $(ARCHDIR)
	@${COMPILER} $*_bound.c -c ${INCLUDE} ${OPTIONS} ${FARGO_OPT}

xmin_bound.c: ${BOUNDARIES} ${GLOBAL}
	@echo PARSING BOUNDARIES...
	@-python ${SCRIPTSDIR}/boundparser.py $^

xmax_bound.c ymin_bound.c ymax_bound.c zmin_bound.c zmax_bound.c: xmin_bound.c

# You may remove the '@' sign in the instructions below if you need a verbose output
%.o : %.c ${GLOBAL}
	@echo CC $< "==>" $*.o
	@/bin/cp -f $< $(ARCHDIR)
	@${COMPILER} $< -c ${INCLUDE} ${OPTIONS} ${CUDAINC} ${FARGO_OPT} -DMONITOR_SCALAR="${MONITOR_SCALAR}" -DMONITOR_2D="${MONITOR_2D}" -DMONITOR_Y="${MONITOR_Y}" -DMONITOR_Y_RAW="${MONITOR_Y_RAW}" -DMONITOR_Z="${MONITOR_Z}" -DMONITOR_Z_RAW="${MONITOR_Z_RAW}" -DSETUPNAME="${SETUP}" -DVERSION=$(GIT_VERSION) -DPATHARCH=$(CURDIR)/$(ARCHDIR)/$(ARCHNAME) -DOPTIONS="${FARGO_OPT}"
