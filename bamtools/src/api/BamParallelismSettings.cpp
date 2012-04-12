#include "BamParallelismSettings.h"

#include <unistd.h>

int BamParallelismSettings::m_configured_threads = 0;
bool BamParallelismSettings::m_multithreading_enabled = false;

int BamParallelismSettings::availableCores()
{

#ifdef OPENMP
  // the openmp function is included here as it is known to be a more reliable way to
  // detect the number of available cores. However, it is only available when compiling
  // with OMP turned on.
	return omp_get_num_procs();
#endif
	
	return sysconf(_SC_NPROCESSORS_ONLN);
}

void BamParallelismSettings::setNumberThreads(int threads)
{
  m_configured_threads = threads;
}

int BamParallelismSettings::getNumberThreads()
{
  if(m_configured_threads == 0)
    return availableCores();
  else
    return m_configured_threads;
}

void BamParallelismSettings::disableMultithreading()
{
  m_multithreading_enabled = false;
}

void BamParallelismSettings::enableMultithreading()
{
  m_multithreading_enabled = true;
}
