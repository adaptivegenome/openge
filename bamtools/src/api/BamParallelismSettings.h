#ifndef BAMPARALLELISMSETTINGS_H
#define BAMPARALLELISMSETTINGS_H

class BamParallelismSettings
{
public:
  // return the number of cores availble in the current system
  static int availableCores();

  // Configure the number of threads to use for each thread pool. Set
  // 0 to use all threads on the system.
  static void setNumberThreads(int threads);
  static int getNumberThreads();

  // Enable or disable the use of multithreading
  static void disableMultithreading();
  static void enableMultithreading();
  static bool isMultithreadingEnabled() { return m_multithreading_enabled; }
  
protected:
  static int m_configured_threads;
  static bool m_multithreading_enabled;
};

#endif
