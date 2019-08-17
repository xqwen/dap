## require at least autoconf 2.63
AC_PREREQ([2.63])

## Process this file with autoconf to produce a configure script.
AC_INIT(dap, m4_esyscmd_s([awk '/^Version:/ {print $2}' DESCRIPTION]))

## Set R_HOME, respecting an environment variable if one is set 
: ${R_HOME=$(R RHOME)}
if test -z "${R_HOME}"; then
    AC_MSG_ERROR([Could not determine R_HOME.])   
fi

## Use R to set CXX and CXXFLAGS
CXX=$(${R_HOME}/bin/R CMD config CXX)
CXXFLAGS=$("${R_HOME}/bin/R" CMD config CXXFLAGS)

## We are using C++
AC_LANG(C++)
AC_REQUIRE_CPP

## Check if GSL is available (required)
AC_CHECK_LIB([m], [cos], [], [echo "m library not found. Please install m library before proceeding."; exit -1])
AC_CHECK_LIB([gslcblas], [cblas_dgemm], [], [echo "gsl library not found. Please install gsl library before proceeding."; exit -1])
AC_CHECK_LIB([gsl], [gsl_blas_dgemm], [], [echo "gsl library not found. Please install gsl library before proceeding."; exit -1])

## Check if -fopenmp is available
allldflags=$(${R_HOME}/bin/R CMD config --ldflags)
hasOpenMP=$(echo ${allldflags} | grep -- -fopenmp)

### Check if is Apple Clang++
apple_clang="no"
RSysinfoName=$("${R_HOME}/bin/Rscript" --vanilla -e 'cat(Sys.info()[["sysname"]])')
if test x"${RSysinfoName}" == x"Darwin"; then
    if test x"${CXX}" == x"clang++"; then
        apple_clang="yes"
    fi
fi
            

## OpenMP
OPENMP_FLAG=""
OPENMP_LIB=""
if test x"${hasOpenMP}" != x""; then
    if test x"${apple_clang}" == x"yes"; then
        OPENMP_FLAG="-Xpreprocessor -fopenmp"
        OPENMP_LIB="-Xpreprocessor -fopenmp -lomp"
    else
        OPENMP_FLAG="-fopenmp"
        OPENMP_LIB="-fopenmp"
    fi
fi


AC_SUBST(LIBS)
AC_SUBST(OPENMP_FLAG)
AC_SUBST(OPENMP_LIB)
AC_CONFIG_FILES([src/Makevars])
AC_OUTPUT

