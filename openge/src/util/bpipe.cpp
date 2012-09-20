/*********************************************************************
 *
 * bpipe.cpp: A bpipe script object
 * Open Genomics Engine
 *
 * Author: Lee C. Baker, VBI
 * Last modified: 12 September 2012
 *
 *********************************************************************
 *
 * This file is released under the Virginia Tech Non-Commercial
 * Purpose License. A copy of this license has been provided in
 * the openge/ directory.
 *
 *********************************************************************/

#include "bpipe.h"

#include <iostream>
#include <fstream>
#include <vector>
#include <algorithm>
#include <iterator>
#include <map>
#include <ctime>

#include <boost/spirit/include/qi.hpp>

#include <boost/fusion/include/adapt_struct.hpp>
#include <boost/spirit/include/phoenix_core.hpp>
#include <boost/spirit/include/phoenix_operator.hpp>
#include <boost/spirit/include/phoenix_stl.hpp>
#include <boost/spirit/include/phoenix_bind.hpp>
#include <boost/spirit/include/phoenix.hpp>
#include <boost/spirit/home/phoenix/object/new.hpp>
using namespace std;

namespace qi = boost::spirit::qi;
namespace ascii = boost::spirit::ascii;
namespace phoenix = boost::phoenix;

struct stage {
    string name;
    string exec;
    
    void setName(const string & name) { this->name = name; }
    
    stage(const stage & s)
    : name(s.name)
    , exec(s.exec)
    {}
    stage(const string & exec) : exec(exec) {}
    stage() {}
};

class StageQueue {
public:
    virtual bool check() { return true; }
    virtual bool execute(string input, string & output) = 0;
    virtual void print() { cerr << "(no description)"; }
};

class ParallelStageQueue : public StageQueue {
    StageQueue * q1, * q2;
public:
    ParallelStageQueue(StageQueue * q1, StageQueue * q2)
    : q1(q1)
    , q2(q2)
    {}
    
    virtual bool check() { return q1->check() && q2->check(); }
    virtual bool execute(string input, string & output) { return q1->execute(input, output) && q2->execute(input, output); }
    virtual void print() { cerr << "Parallel(";q1->print();cerr << ","; q2->print(); cerr << ")"; }
};

class SerialStageQueue : public StageQueue {
    StageQueue * q1, * q2;
public:
    SerialStageQueue(StageQueue * q1, StageQueue * q2)
    : q1(q1)
    , q2(q2)
    {}
    virtual bool check() { return q1->check() && q2->check(); }
    virtual bool execute(string input, string & output) { string intermediate_filename; return q1->execute(input, intermediate_filename) && q2->execute(intermediate_filename, output); }
    virtual void print() { cerr << "Serial(";q1->print();cerr << ","; q2->print(); cerr << ")"; }
};

class StageReference : public StageQueue {
    string name, command, filter, input, output;
    vector<stage> & stages;
    stage * s;
public:
    StageReference(string name, vector<stage> & stages)
    : name(name)
    , stages(stages)
    , s(NULL)
    {}
    
    virtual bool check();
    virtual bool execute(string input, string & output);
    virtual void print() { cerr << name; }
};

bool StageReference::check() {
    for(vector<stage>::iterator si = stages.begin(); si != stages.end(); si++)
        if(si->name == name) {
            s = &*si;
            break;
        }
    
    if(!s) {
        cerr << "BPipe file error: stage name '" << name << "' didn't match any known stages." << endl;
        return false;
    }
    
    command = s->exec;
    
    output = input + "." + name;

    if(!filter.empty()) {
        string::iterator dot = find(input.begin(), input.end(), '.');
        if(dot != input.end()) {
            output = string(input.begin(), dot) + "." + filter + string(dot, input.end());
        }
    }
    
    //TODO replace $input and $output

    return true;
}

bool StageReference::execute(string input, string & output) {
    
    time_t now = time(NULL);
    cerr << "=== Stage " << name << " " << ctime(&now) << " ===" << endl;
    int ret = system(command.c_str());
    
    if(0 != ret)
        cerr << "Execution of stage failed." << endl;
    
    return 0 == ret;
}

template <typename Iterator>
struct BpipeParser : qi::grammar<Iterator, vector<stage>(), ascii::space_type>
{
    vector<stage> stages;
    map<string, stage> stage_names;
    StageQueue * run_task;
    
    static stage & setStageName(string name, stage & s) {
        s.name = name;
        return s;
    }
    
    BpipeParser() : BpipeParser::base_type(start)
    {
        using qi::lit;
        using qi::lexeme;
        using ascii::char_;
        using qi::_val;
        using qi::_1;
        using qi::_2;
        using qi::space;
        using qi::alnum;
        using phoenix::ref;
        using phoenix::bind;
        using phoenix::construct;
        using phoenix::val;
        using phoenix::push_back;
        using phoenix::new_;
        
        run_task = NULL;

        quoted_string.name("quoted_string");
        unquoted_string.name("unquoted_string");
        stage_block.name("stage_block");
        run_block.name("run_block");
        bpipe_file.name("bpipe_file");

        quoted_string = lexeme['"' >> +(char_ - '"') >> '"'];
        unquoted_string %= +(lit("\\\"")[_val='"'] | lit("\\\\")[_val='\\'] | alnum);
        exec_statement %= lit("exec") >> quoted_string ;
        stage_block %= '{' >> exec_statement[_val = construct<stage>(_1)] >>  '}' ;
        stage_assignment = (unquoted_string >> lit("=") >> stage_generator)[_val = bind(&setStageName, _1, _2)];
        stage_generator %= stage_block | stage_assignment;
        stage_definition = stage_generator[push_back(ref(stages), _1)];
        
        stage_reference.name("StageReference");
        stage_serial_queue.name("StageSerialQueue");
        stage_parallel_queue.name("StageParallelQueue");
        stage_queue.name("StageQueue");
        run_block.name("RunBlock");

        stage_reference = unquoted_string [_val = new_<StageReference>(_1, ref(stages))];
        stage_serial_queue = stage_reference[_val = _1] >> -(('+' >> stage_queue)[_val = new_<SerialStageQueue>(_val, _1)]);
        stage_parallel_queue = '[' >> stage_queue[_val = _1] >> *((',' >> stage_queue)[_val = new_<ParallelStageQueue>(_val, _1)]) >> ']';
        stage_queue %= (stage_serial_queue | stage_parallel_queue | stage_reference)[_val = _1];
        run_block = (lit("Bpipe.run") >> '{' >> stage_queue >> '}')[ref(run_task) = _1] ;
        bpipe_file %= *(stage_definition) >> run_block;
        start %= bpipe_file;
    }

    //stage definitions
    qi::rule<Iterator, string(), ascii::space_type> quoted_string;
    qi::rule<Iterator, string(), ascii::space_type> unquoted_string;
    qi::rule<Iterator, string(), ascii::space_type> exec_statement;
    qi::rule<Iterator, stage(), ascii::space_type> stage_generator;
    qi::rule<Iterator, ascii::space_type> stage_definition;
    qi::rule<Iterator, stage(), ascii::space_type> stage_block;
    qi::rule<Iterator, stage(), ascii::space_type> stage_assignment;
    
    //run
    qi::rule<Iterator, ascii::space_type> run_block;
    qi::rule<Iterator, StageQueue*(), ascii::space_type> stage_reference;
    qi::rule<Iterator, StageQueue*(), ascii::space_type> stage_serial_queue;
    qi::rule<Iterator, StageQueue*(), ascii::space_type> stage_parallel_queue;
    qi::rule<Iterator, StageQueue*(), ascii::space_type> stage_queue;
    
    //file
    qi::rule<Iterator, vector<stage>(), ascii::space_type> bpipe_file;
    qi::rule<Iterator, vector<stage>(), ascii::space_type> start;
    
};

BPipe::BPipe() {
    
}

bool BPipe::load(const string & filename) {
    this->filename = filename;
    
    ifstream file(filename.c_str());
    
    if(file.fail()) {
        cerr << "Error opening file " << filename << ". Aborting." << endl;
        exit(-1);
    }
    
    //read in script
    string line;
    while(getline(file, line))
        script_text += line + "\n";

    file.close();
    
    
    //remove  /**/ comments
    while(true) {
        size_t cstart = script_text.find("/*");
        size_t cend = script_text.find("*/");
        
        if(cstart != string::npos && cend != string::npos)
            script_text.erase(cstart, cend - cstart +2);
        else
            break;
    }
    
    // remove // comments
    while(true) {
        size_t cstart = script_text.find("//");
        size_t cend = script_text.find("\n", cstart);
        
        if(cstart != string::npos && cend != string::npos)
            script_text.erase(cstart, cend - cstart);
        else
            break;
    }
    
    return true;
}

bool BPipe::check() {
    using qi::double_;
    using qi::phrase_parse;
    using ascii::space;
    
    string::const_iterator first = script_text.begin();
    string::const_iterator last = script_text.end();
    parser = new BpipeParser<string::const_iterator>();
    vector<stage> result;
    bool r = phrase_parse( first, last, *parser, space, result);
    
    string::const_iterator p = first;
    string::const_iterator q = script_text.begin();
    int d = distance(q, p);
    if (first != last) { // fail if we did not get a full match
        cerr << "Parsed up to " << string(script_text, d) << endl;
        return false;
    }
    else {
        return parser->run_task->check();
    }
    return false;
}

void BPipe::print() {
    
    parser->run_task->print();
}

bool BPipe::execute() {
    time_t start_time = time(NULL);
    cerr << "=== Starting pipeline at " << ctime(&start_time) << " ===" << endl;
    string output;
    bool ret = parser->run_task->execute("debug_asdf.txt", output);
    
    time_t stop_time = time(NULL);
    
    if(!ret)
        cerr << "=== Pipeline FAILED at " << ctime(&stop_time) << " ===" << endl;
    else {
        cerr << "=== Finished successfully at " << ctime(&stop_time) << " ===" << endl;
    }
    return ret;
}