/* This is the general preprocess include file. It should be included in every fortran code at
   the place where the USE statements are located. This ensures the correct working of all the 
   macros and defines located here and supplied by the make files... 
*/

/* formatted string with date and time */
#define __TIME_VARIABLE__ date
#define __TIME_STRING__ character(len=30) :: __TIME_VARIABLE__
#define __GET_TIME__ CALL fdate(__TIME_VARIABLE__)

/*==============================================================================*
 * Define alias for log file management.                                        *
 *==============================================================================*/

#if defined(LOG_FILE) 
#define __LOG_UNIT              19
#define __OPEN_LOG_FILE         OPEN( UNIT=__LOG_UNIT, FILE=LOG_FILE, \
                                      ACTION="write", POSITION="append" )
#define __CLOSE_LOG_FILE        CLOSE( UNIT=__LOG_UNIT )
#define __WRITE_TO_LOG          WRITE(__LOG_UNIT,*) __FILE__," | ",

#define __INIT_LOG_FILE \
OPEN(UNIT=__LOG_UNIT, FILE=LOG_FILE, ACTION="write"); \
WRITE(__LOG_UNIT,"(/,A)") " Log file for LEPORELLO run"; \
WRITE(__LOG_UNIT,*) "Version of the code: ", VERSIONTAG; \
__GET_TIME__; \
WRITE(__LOG_UNIT,*) "Program started on ",__TIME_VARIABLE__; \
CLOSE(UNIT=__LOG_UNIT)

#define __END_LOG_FILE \
OPEN(UNIT=__LOG_UNIT, FILE=LOG_FILE, ACTION="write", POSITION="append"); \
__GET_TIME__; \
WRITE(__LOG_UNIT,"(A,/)") " Program ended on "//__TIME_VARIABLE__; \
CLOSE(UNIT=__LOG_UNIT)
#endif

/*==============================================================================*
 * Define/undefine printing of spectral density.                                *
 *==============================================================================*/

#undef __PRINT_SPECTRAL_DENSITY

/*==============================================================================*
 * Define alias for OMP parallelization.                                        *
 *==============================================================================*/

#if defined(WITH_OPENMP)
#define __OMP_TotalNrOfThreads      OMP_GET_NUM_THREADS()
#define __OMP_CurrentThreadNum      OMP_GET_THREAD_NUM() + 1
#define __OMP_OnlyMasterBEGIN       IF ( OMP_GET_THREAD_NUM() == 0 ) THEN;
#define __OMP_OnlyMasterEND         ; END IF
#else
#define __OMP_TotalNrOfThreads      1
#define __OMP_CurrentThreadNum      1
#define __OMP_OnlyMasterBEGIN   
#define __OMP_OnlyMasterEND     
#endif

/*==============================================================================*
 * Load common modules.                                                         *
 *==============================================================================*/

   /* These are modules that are used in most parts of the code */
   USE ErrorTrap
   USE MyConsts
   USE MyLinearAlgebra
#if defined(WITH_OPENMP)
   USE omp_lib
#endif
   
