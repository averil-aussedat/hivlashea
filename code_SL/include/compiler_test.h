#ifndef SELA_COMPILER_TEST
#define SELA_COMPILER_TEST

// TODO: Some ICC versions overlap in time, so a simple >= is not enough.
// https://software.intel.com/en-us/articles/intel-compiler-and-composer-update-version-numbers-to-compiler-version-number-mapping

// Check for compiler version, to be able to use OpenMP simd pragmas.
// http://www.openmp.org/resources/openmp-compilers/
// OpenMP 4.0 C/C++/Fortran supported in version 15.0 and 16.0 compilers
// From GCC 4.9.1, OpenMP 4.0 is fully supported
#if defined(__ICC) || defined(__INTEL_COMPILER) // Intel ICC/ICPC
#    if (__INTEL_COMPILER_BUILD_DATE >= 20140726)
#        define SELA_OPENMP_4_0
#    endif
#elif defined(__GNUC__)                         // GCC
#    if ((__GNUC__ > 4) || (__GNUC__ == 4 && (__GNUC_MINOR__ >= 9)))
#        define SELA_OPENMP_4_0
#    endif
#else
#    define SELA_OPENMP_4_0
#endif

// Check for compiler version, to be able to use features from C11.
// Pic-Vert currently uses _Generic and _Alignas.
// https://software.intel.com/en-us/articles/c11-support-in-intel-c-compiler
// https://gcc.gnu.org/wiki/C11Status
#if defined(__ICC) || defined(__INTEL_COMPILER) // Intel ICC/ICPC
#    if (__INTEL_COMPILER_BUILD_DATE >= 20150805)
#        define SELA_C11
#    endif
#elif defined(__GNUC__)                         // GCC
#    if ((__GNUC__ > 4) || (__GNUC__ == 4 && (__GNUC_MINOR__ >= 9)))
#        define SELA_C11
#    endif
#else
#    define SELA_C11
#endif

// Check for MPI version, to be able to use MPI_Iallreduce from MPI-3.
// MPICH v3.0 implements the MPI-3.0 standard.
#if defined(OPEN_MPI)                  // OpenMPI
#    if (MPI_VERSION >= 3)
#        define SELA_MPI_3_0
#    endif
#elif defined(MPICH)                   // MPICH
#    if (MPICH_NUMVERSION >= 30000000)
#        define SELA_MPI_3_0
#    endif
#else
#    define SELA_MPI_3_0
#endif

#endif // ifndef SELA_COMPILER_TEST
