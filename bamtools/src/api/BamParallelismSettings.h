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

  
class BamSpinlock
{
public:
    BamSpinlock()
    : lock_holder(0) {}
    // for an explanation of the locks used here, see
    // http://stackoverflow.com/questions/1383363/is-my-spin-lock-implementation-correct-and-optimal
    void lock() {
        while (__sync_lock_test_and_set(&lock_holder, 1)) while(lock_holder);
        
    }
    void unlock() {
        __sync_synchronize();
        lock_holder = 0;
    }
protected:
    volatile char lock_holder;
};

#endif
