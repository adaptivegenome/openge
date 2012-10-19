#ifndef BAMPARALLELISMSETTINGS_H
#define BAMPARALLELISMSETTINGS_H
  
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
