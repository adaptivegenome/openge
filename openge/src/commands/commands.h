#ifndef COMMANDS_H
#define COMMANDS_H
/*********************************************************************
 *
 * commands.h: Base class for an OpenGE command.
 * Open Genomics Engine
 *
 * Author: Lee C. Baker, VBI
 * Last modified: 21 May 2012
 *
 *********************************************************************
 *
 * This file is released under the Virginia Tech Non-Commercial 
 * Purpose License. A copy of this license has been provided in 
 * the openge/ directory.
 *
 *********************************************************************/

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
    std::string output_filename;
    
    // Directory in which to put temporary files. Recommended to be fast storage.
    std::string tmpdir;
    
    // Multithreading is disabled
    bool nothreads;
    
    // Number of threads per threadpool.
    unsigned int num_threads;
    
public:
    bool isVerbose() const { return verbose; }
    bool isNothreads() const { return nothreads; }
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

class DedupCommand: public OpenGECommand
{
protected:
    virtual int runCommand();
    void getOptions();
};

class HelpCommand: public OpenGECommand
{
protected:
    virtual int runCommand();
    void getOptions();
};

class LocalRealignCommand: public OpenGECommand
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

class StatsCommand: public OpenGECommand
{
protected:
    virtual int runCommand();
    void getOptions();
};

class VersionCommand: public OpenGECommand
{
protected:
    virtual int runCommand();
};

class ViewCommand: public OpenGECommand
{
protected:
    virtual int runCommand();
    void getOptions();
};

#endif