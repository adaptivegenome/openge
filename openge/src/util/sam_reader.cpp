//
//  SamReader.cpp
//  BamTools
//
//  Created by Lee Baker on 3/16/12.
//  Copyright (c) 2012 LCB. All rights reserved.
//

#include <iostream>
#include <sstream>
#include <cstdio>
#include "sam_reader.h"

#include <api/BamParallelismSettings.h>

#ifdef __linux__
#include <sys/prctl.h>
#endif

#include <cstring>
using namespace std;
using namespace BamTools;

const int MAX_LINE_QUEUE_SIZE = 6000;

/*
class SamParseJob : public ThreadJob
{
protected:
    SamLine_t & data;
public:
    SamParseJob(SamLine_t & data) : data(data) { }
    virtual void runJob();
};

void SamParseJob::runJob()
{
#ifdef __linux__
    prctl(PR_SET_NAME,"SAM_line_worker",0,0,0);
#endif
    data.al = data.reader->ParseAlignment(data.line);
    data.parsed = true;
}
 */

SamReader::SamReader()
: file(NULL)
, loaded(false)
, finished(false)
{ }


void * SamReader::LineWorkerThread(void * reader_p)
{
#ifdef __linux__
    prctl(PR_SET_NAME,"SAM_line_worker",0,0,0);
#endif
    SamReader * reader = (SamReader *) reader_p;
    while(1) {
        SamLine * data = NULL;

        //if there is a job, take it. Otherwise sleep for a bit.
        while(!reader->workers_finished) {
            reader->worker_jobs_lock.lock();
            if(reader->jobs_for_workers.size() > 0) {
                data = reader->jobs_for_workers.pop();
                reader->worker_jobs_lock.unlock();
                break;
            } else {
                reader->worker_jobs_lock.unlock();
                usleep(2000);
            }
            
        }
        
        if(reader->workers_finished)
            break;
        
        data->al = reader->ParseAlignment(data->line);
        data->parsed = true;
    }

    return 0;
}

void * SamReader::LineGenerationThread(void * data)
{
#ifdef __linux__
    prctl(PR_SET_NAME,"SAM_line_master",0,0,0);
#endif
    SamReader * reader = (SamReader *) data;
    
    FILE * fp = fopen(reader->filename.c_str(), "r");
    if(!fp) {
        cerr << "Error opening file for line parser" << endl;
        abort();
    }
    
    //open our C file handle to the same position that the ifstream was at.
    size_t position = reader->file.tellg();
    fseek(fp, position, SEEK_SET);
    
    char line_s[10240];
    
    while(!reader->finished) {
        while( reader->jobs.size() > MAX_LINE_QUEUE_SIZE)
            usleep(20000);
        
        char * read = fgets(line_s, 1024, fp);
        
        if(!read || strlen(read) < 10)  // if line is invalid
            break;
        
        SamLine_t * lt = new SamLine_t;
        
        //create a buffer for the worker thread to use.
        lt->line = (char *) malloc(strlen( line_s) + 1);
        strcpy(lt->line, line_s);

        lt->reader = reader;
        reader->jobs.push(lt);
        //SamParseJob * parse_job = new SamParseJob(*lt);
        //reader->pool.addJob(parse_job);
        reader->jobs_for_workers.push(lt);
    }
    
    fclose(fp);

    reader->finished = true;

    return NULL;
}

bool SamReader::Open(const string & filename)
{
    this->filename = filename;
    file.open(filename.c_str(), ios::in);
    
    LoadHeaderData();
    
    loaded = true;
    workers_finished = false;
    
    for(int i = 0; i < BamParallelismSettings::getNumberThreads(); i++)
    {
        pthread_t t;
        pthread_create(&t, NULL, LineWorkerThread, this);
        worker_threads.push_back(t);
    }
    pthread_create(&line_generation_thread, 0, LineGenerationThread, this);
    
    return true;
}

bool SamReader::Close()
{
    workers_finished = true;
    for(int i = 0; i < BamParallelismSettings::getNumberThreads(); i++)
        pthread_join(worker_threads[i], NULL);

    worker_threads.clear();
    
    finished = true;
    pthread_join(line_generation_thread, NULL);
    file.close();
    return true;
}

SamHeader SamReader::GetHeader() const
{
    return header;
}

bool SamReader::IsLoaded()
{
    return loaded;
}

// retrieves header text from BAM file
void SamReader::LoadHeaderData(void)
{
    stringstream header_txt("");
    
    while(file.peek() == '@')
    {
        string line;
        getline(file, line);
        header_txt << line << endl;
    }
    
    header.SetHeaderText(header_txt.str());
    
    m_refData.reserve(header.Sequences.Size());
    
    for(SamSequenceConstIterator it = header.Sequences.Begin(); it != header.Sequences.End(); it++)
    {
        //convert length to an int from a string.
        int length;
        stringstream length_ss(it->Length);
        length_ss >> length;
        
        m_refData.push_back(RefData(it->Name, length));
    }
}

bool AddHeaderAttributeArray(BamAlignment & alignment, const string & tag,const string & value);

// retrieves BAM alignment under file pointer
// (does no overlap checking or character data parsing)
BamAlignment * SamReader::LoadNextAlignment()
{
#if 1 
    
    //wait for there to be something at the back of the queue that is done.
    while(true) {
        if(jobs.size() && jobs.front() && jobs.front()->parsed)
            break;

        //if there is nothing left in the queue, we are done
        if(finished && 0 == jobs.size())    
            return NULL;
        usleep(20);
    }
    //cerr << "pop" << endl;
    SamLine_t * s = jobs.pop();
    BamAlignment * ret = s->al;
    free(s->line);
    delete s;
    return ret;
#else
    
    if(!loaded)
        return NULL;
    
    if(file.eof())
        return NULL;
    
    string line_s;
    getline(file, line_s);
    return ParseAlignment(line_s);
#endif
}

BamAlignment * SamReader::ParseAlignment(const string & line_s)
{
    BamAlignment * al = new BamAlignment();
    BamAlignment & alignment = *al;
    
    string & qname = alignment.Name;
    uint32_t & flag = alignment.AlignmentFlag;
    //string rname;
    int & pos = alignment.Position;
    uint16_t & mapq = alignment.MapQuality;
    //string cigar;
    //string rnext;
    int & pnext = alignment.MatePosition;
    int & tlen = alignment.InsertSize;
    string & seq = alignment.QueryBases;
    string & qual = alignment.Qualities;
    
    char * line = (char *) malloc(sizeof(char) * (line_s.length() + 1));
    memcpy(line, line_s.c_str(), line_s.length());
    line[line_s.length()] = 0;  //null terminate
    char * field_starts[12] = {0};
    
    for(int ct = 0; ct < 11; ct++)  //for 11 columns, null terminate each field and put in array of field
    {
        field_starts[ct] = line;
        line = strchr(line, '\t');
        if(ct < 10) {
            if(!line || line >= line + line_s.length())
                return NULL;    
            if(*line) {
                *line = 0;
                line++;
            }
        }
    }
    
    qname.assign(field_starts[0]);
    flag = atoi(field_starts[1]);
    const char * rname = field_starts[2];
    pos = atoi(field_starts[3]);
    mapq = atoi(field_starts[4]);
    const char * cigar = field_starts[5];
    const char * rnext = field_starts[6];
    pnext = atoi(field_starts[7]);
    tlen = atoi(field_starts[8]);
    seq.assign(field_starts[9]);
    qual.assign(field_starts[10]);

    // zero based indexes:
    alignment.Position--;
    alignment.MatePosition--;

    // rname
    if(rname[0] == '*')
        alignment.RefID = -1;
    else if(header.Sequences.Contains(rname))
        alignment.RefID = header.Sequences.IndexOfString(rname);
    else {
        cerr << "Rname " << rname << " missing from sequence dictionary" << endl;
        alignment.RefID = -1;
    }

    // rnext
    if(rnext[0] == '=')
        alignment.MateRefID = alignment.RefID;
    else if(rnext[0] == '*')
        alignment.MateRefID = -1;
    else if(header.Sequences.Contains(rnext))
        alignment.MateRefID = header.Sequences.IndexOfString(rnext);
    else {
        cerr << "RNext " << rnext << " missing from sequence dictionary" << endl;
        alignment.MateRefID = -1;
    }

    //CIGAR ops
    if(cigar[0] != '*')
    {
        int cigar_ct = 0;
        char * cigar_p = (char *) cigar;    //cast away const-ness, because we will be changing the ptr (not the data)
        while(isdigit(*cigar_p))
        {
            int num;
            char c;
            num = strtol(cigar_p, &cigar_p, 10);
            c = *cigar_p;
            cigar_p++;
            alignment.CigarData.push_back(CigarOp(c,num));
            cigar_ct ++;
            assert(cigar_ct < 100);
        }
    }
    
    if(!line)
        return NULL;
    
    vector<char *> additional_parameters;
    
    while(line && *line)
    {
        additional_parameters.push_back(line);
        line = strchr(line, '\t');
        if(!line)
            break;
        *line = 0;
        line++;
    }
    
    // cache string objects for some types to ensure they don't have to be constructed
    // every time they are encountered.
    static const string typeZ("Z");
    static const string typei("i");
    //optional attributes
    for (int i = 1; i < additional_parameters.size(); i++) {
        char * segment = additional_parameters[i];
        
        //null terminate field separators (:) in the segment so we can use to construct strings
        segment[2] = 0;
        segment[4] = 0;
        
        const char * tag(segment);
        const char * type(&segment[3]);
        const char * value(&segment[5]);
        
        bool retval = false;
        switch(*type)
        {
            case 'i': //int
            {
                int i = atoi(value);
                retval = alignment.AddTag(tag, typei, i);
                break;
            }
            case 'f': //floating point
            {
                float f = atof(value);
                retval = alignment.AddTag(tag, type, f);
                break;
            }
            case 'A': //single char
            {
                retval = alignment.AddTag(tag, type, *value);
                break;
            }
            case 'Z': //string, with spaces
                retval = alignment.AddTag(tag, typeZ, string(value));
                break;
            case 'B': //byte array
                retval = alignment.AddTag(tag, type, value);
                break;
            case 'H': //array type
                retval = AddHeaderAttributeArray(alignment, tag, value);
                break;
            default:
                cerr << "Invalid attribute type " << type << endl;
                break;
        }
        
        if(true != retval)
            cerr << "Error parsing optional attribute " << tag << " in line " << segment << endl;
    }
    
    free( field_starts[0] );
    
    return al;
}

template <class T>
vector<T> ParseAttributeArray(const string & value)
{
    vector<T> array;
    T var;
    
    //cut off the first letter, it indicates the type and is handled below
    stringstream ss(value);
    
    string element;
    getline(ss, element, ',');  //read in initial letter and comma
    
    while(!ss.eof()) //now add each element to the array
    {
        getline(ss, element, ',');
        stringstream sse(element);
        sse >> var;
        array.push_back(var);
    }
    return array;
}

bool AddHeaderAttributeArray(BamAlignment & alignment, const string & tag,const string & value)
{
    switch(value[0])
    {
        case 'i':
            return alignment.AddTag(tag, ParseAttributeArray<int>(value));
        case 'I':
            return alignment.AddTag(tag, ParseAttributeArray<unsigned int>(value));
        case 's':
            return alignment.AddTag(tag, ParseAttributeArray<short>(value));
        case 'S':
            return alignment.AddTag(tag, ParseAttributeArray<unsigned short>(value));
        case 'c':
            return alignment.AddTag(tag, ParseAttributeArray<char>(value));
        case 'C':
            return alignment.AddTag(tag, ParseAttributeArray<unsigned char>(value));
        case 'f':
            return alignment.AddTag(tag, ParseAttributeArray<float>(value));
        default:
            cerr << "Invalid array type for optional parameter " << tag << endl;
            return false;
    }
}

const RefVector & SamReader::GetReferenceData(void)
{
    return m_refData;
}

// retrieves next available alignment
bool SamReader::GetNextAlignment(BamAlignment& alignment)
{
    BamAlignment * ret = LoadNextAlignment();
    
    if(!ret)
        return false;
    alignment = *ret;
    return true;
}
BamAlignment * SamReader::GetNextAlignment()
{
    return LoadNextAlignment();
}