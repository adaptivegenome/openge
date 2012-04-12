#ifndef COMMANDS_H
#define COMMANDS_H

#include <vector>
#include <string>

#include <boost/program_options.hpp>
using namespace boost;

class HelpCommand;
class OpenGECommand
{
    friend class HelpCommand;
public:
    OpenGECommand();
    ~OpenGECommand();
    int runWithParameters(int argc, const char ** argv);
protected:
    // Inheriting classes can override this function to add command line options, using
    // the 'options' member variable.
    virtual void getOptions() {};
    
    // Inheriting classes should implement this option if they wish to be notified when options have been processed.
    virtual void processOptions() {};
    
    // This function performs the work in each command, and must be implemented. The return value
    // is returned as the application's return value.
    virtual int runCommand() = 0;
    
    program_options::positional_options_description options_positional;
    program_options::options_description options;
    program_options::variables_map vm;
    
    // Automatically set- when true, include extra description of what is happening. Eg- progress indicators,
    // extra detail.
    bool verbose;
    
    // Automatically filled with filenames to be used as input. Can also include stdin.
    std::vector<std::string> input_filenames;
    
    // Multithreading is disabled
    bool nothreads;
    
    // Number of threads per threadpool.
    unsigned int num_threads;
};
class CommandMarshall
{
public:
    static OpenGECommand * commandWithName(const std::string name);
protected:
    friend class HelpCommand;
};

/////////////////////////////////
// Actual commands

class ConvertCommand : public OpenGECommand
{
protected:
    void getOptions();
    virtual int runCommand();
};

class CountCommand : public OpenGECommand
{
protected:
    void getOptions();
    virtual int runCommand();
};

class HeadCommand : public OpenGECommand
{
protected:
    void getOptions();
    virtual int runCommand();
};

class HelpCommand: public OpenGECommand
{
protected:
    virtual int runCommand();
    void getOptions();
};

class MergeSortCommand : public OpenGECommand
{
protected:
    void getOptions();
    virtual int runCommand();
    
    class MergeSortCommandImplementation;
    MergeSortCommandImplementation * impl;
};

#endif