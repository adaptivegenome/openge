/*********************************************************************
 *
 * sam_reader.cpp: A multithreaded sequential SAM file reader.
 * Open Genomics Engine
 *
 * Author: Lee C. Baker, VBI
 * Last modified: 16 March 2012
 *
 *********************************************************************
 *
 * This file is released under the Virginia Tech Non-Commercial 
 * Purpose License. A copy of this license has been provided in 
 * the openge/ directory.
 *
 *********************************************************************/


#include <iostream>
#include <sstream>
#include <cstdio>
#include <cassert>
#include "sam_reader.h"
#include <errno.h>
#include <climits>

#include <fcntl.h>           /* For O_* constants */
#include <sys/stat.h>        /* For mode constants */
#include <semaphore.h>

#include <cstring>
using namespace std;

#define SAM_READER_MT

const int MAX_LINE_QUEUE_SIZE = 6000;

SamReader::SamReader()
: loaded(false)
, finished(false)
, lines_since_last_sem_unlock(0)
{ }

#ifdef SAM_READER_MT

// The worker thread parses line by line. There is a queue
// distributing jobs to the worker threads, and this function also
// throttles the worker pool using semaphores to block the threads that are
// simply polling an empty queue. This ensures that most of the time, 
// there is work for the workers in the queue to grab without needing to block
// the thread (and wasting time). One worker (the first) never acquires the semaphore-
// this way we ensure that there is always at least one worker present.
void * SamReader::LineWorkerThread(void * reader_p)
{
    ogeNameThread("SAM_line_worker");

    SamReader * reader = (SamReader *) reader_p;
    while(1) {
        SamLine * data = NULL;
        bool post_a_semaphore = false;

        //if there is a job, take it. Otherwise sleep for a bit.
        while(!reader->workers_finished) {
            bool wait_a_semaphore = false;
            {
                reader->worker_jobs_lock.lock();
                
                size_t size = reader->jobs_for_workers.size();
                reader->lines_since_last_sem_unlock++;
                if(0 == reader->lines_since_last_sem_unlock % 10000 && size > 2000 && reader->active_workers < reader->num_worker_threads) {
                    post_a_semaphore = true;
                    reader->active_workers++;
                    reader->lines_since_last_sem_unlock = 0;
                }

                //reverify size just in case the last element was popped between the size() and lock() above
                if(0 != size)
                    data = reader->jobs_for_workers.pop();
                wait_a_semaphore = (data == NULL && reader->active_workers > 1);

                if(wait_a_semaphore)
                    reader->active_workers--;

                reader->worker_jobs_lock.unlock();
            }
            
            if(data)    //success, got a job.
                break;

            // ensure the first worker thread cannot block- that way we ensure that not all
            // threads block if the queue is empty for a while.
            if(wait_a_semaphore) {
                if(0 != sem_wait(reader->sam_worker_sem))
                    perror("Error waiting for sam_worker semaphore");
            }
        }
        
        // if the queue gets too big, allow more workers to work.
        // since a partially populated queue is much better than an empty queue
        // we need to make sure this can't happen often. This means making it harder to
        // unblock a worker than to block one.
        //
        // The best way without locking that I could think of is to only unlock if the
        // write counter is divisible by some value- so a 1 in N probability that another 
        // thread will be unblocked.
        if(post_a_semaphore) {
            if(0 != sem_post(reader->sam_worker_sem))
                perror("Error posting sam_worker semaphore");
        }
        
        if(reader->workers_finished)
            break;

        assert(data);
        assert(data->line);

        data->al = reader->ParseAlignment(data->line);
        data->parsed = true;
    }
    
    // must post a semaphore in case all remaining workers are waiting.
    if(0 != sem_post(reader->sam_worker_sem))
        perror("Error posting sam_worker semaphore while quitting.");

    return 0;
}

void * SamReader::LineGenerationThread(void * data)
{
    ogeNameThread("SAM_line_master");

    SamReader * reader = (SamReader *) data;
    
    FILE * fp = NULL;
    
    if(reader->filename == "stdin")
        fp = stdin;
    else
        fp = fopen(reader->filename.c_str(), "r");
    if(!fp) {
        cerr << "Error opening file for line parser (" << reader->filename << ")." << endl;
        exit(-1);
    }

    stringstream header_txt("");
    bool parsing_header = true;

    reader->finished = false;

    char line_s[10240];
    while(!reader->finished) {
        while( reader->jobs.size() > MAX_LINE_QUEUE_SIZE)
            usleep(20000);
        
        char * read = fgets(line_s, sizeof(line_s)-1, fp);
        
        if(!read)
            break;

        if(ferror(fp)) {
            cerr << "Error reading lines from SAM file." << endl;
            if(errno)
                perror("Error reading line: ");
            break;
        }

        if(line_s[0] == '@') {
            header_txt << string(line_s);
        } else {
            if(parsing_header) {
                reader->header = BamHeader(header_txt.str());

                parsing_header = false;
                reader->loaded = true;
            }

            size_t read_len = strlen(line_s);
            line_s[read_len-1] = 0;   //remove newline from string

            if(read_len < 10) { // if line is invalid
                cerr << "read len end" << endl;
                break;
            }
            
            SamLine * lt = new SamLine;
            
            //create a buffer for the worker thread to use.
            lt->line = lt->line_static;
            
            if(read_len >= sizeof(lt->line_static))
                lt->line = (char *) malloc(read_len + 1);

            assert(lt->line);
            strcpy(lt->line, line_s);

            reader->jobs.push(lt);
            reader->jobs_for_workers.push(lt);
        }
        
        line_s[0] = 0;
        if(feof(fp))
            reader->finished = true;
    }
    
    if(fp != stdin)
        fclose(fp);

    reader->finished = true;

    return NULL;
}
#endif

bool SamReader::open(const string & filename)
{
    this->filename = filename;
    
#ifdef SAM_READER_MT
    workers_finished = false;
    num_worker_threads = OGEParallelismSettings::getNumberThreads();
    active_workers = num_worker_threads;
    
    int32_t sem_id = 0xffffffff & (int64_t) this;
    sprintf(sam_worker_sem_name, "sam_wkr_%x",sem_id);
    if(-1 == sem_unlink(sam_worker_sem_name) && errno != ENOENT);   //success of this function doesn't matter; there's nothing we can do about it.

    sam_worker_sem = sem_open(sam_worker_sem_name, O_CREAT | O_EXCL,0700,0);
    
	if(sam_worker_sem == SEM_FAILED &&  0 != errno) {
		perror("Error opening SAM worker semaphore");
		assert(0);
	}
    
    for(int i = 0; i < OGEParallelismSettings::getNumberThreads(); i++)
    {
        pthread_t t;
        pthread_create(&t, NULL, LineWorkerThread, this);
        worker_threads.push_back(t);
    }
    pthread_create(&line_generation_thread, 0, LineGenerationThread, this);

#endif
    return true;
}

void SamReader::close()
{
#ifdef SAM_READER_MT
    workers_finished = true;
    for(int i = 0; i < OGEParallelismSettings::getNumberThreads(); i++)
        pthread_join(worker_threads[i], NULL);

    worker_threads.clear();
    
    finished = true;
    pthread_join(line_generation_thread, NULL);
    
    if(0 != sem_close(sam_worker_sem))
        perror("Error closing sam_worker semaphore");
	if(0 != sem_unlink(sam_worker_sem_name))
        perror("Error unlinking sam_worker semaphore");
#endif
}

const BamHeader & SamReader::getHeader() const
{
    while (!loaded)
        usleep(20000);

    return header;
}

bool AddHeaderAttributeArray(OGERead & alignment, const string & tag,const string & value);

// retrieves BAM alignment under file pointer
// (does no overlap checking or character data parsing)
OGERead * SamReader::LoadNextAlignment()
{
#ifdef SAM_READER_MT
    
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
    SamLine * s = jobs.pop();
    OGERead * ret = s->al;
    if(s->line != s->line_static)
        free(s->line);

    s->line = NULL;
    delete s;
    return ret;
#else
    //LCB - this single threaded code should be kept for backwards compat and --nothreads, but needs to be verified for correctness.
    if(!loaded)
        return NULL;
    
    if(file.eof())
        return NULL;
    
    string line_s;
    getline(file, line_s);
    return ParseAlignment(line_s.c_str());
#endif
}

OGERead * SamReader::ParseAlignment(const char * line_s) const
{
    OGERead * al = OGERead::allocate();
    OGERead & alignment = *al;

    size_t line_length = strlen(line_s);
    if(line_length < 10)    //if the line is shorter than 10 chars, it is definitely not a full SAM line.
        return NULL;
    char * line = (char *) malloc(sizeof(char) * (line_length + 1));
    memcpy(line, line_s, line_length);
    line[line_length] = 0;  //null terminate
    char * field_starts[12] = {0};
    
    for(int ct = 0; ct < 11; ct++)  //for 11 columns, null terminate each field and put in array of field
    {
        field_starts[ct] = line;
        line = strchr(line, '\t');
        if(ct < 10 && (!line || line >= line + line_length))
            return NULL;    
        if(line && *line) {
            *line = 0;
            line++;
        }
    }
    
    alignment.setName(field_starts[0]);
    alignment.setAlignmentFlag(atoi(field_starts[1]));
    const char * rname = field_starts[2];
    alignment.setPosition(atoi(field_starts[3]));
    alignment.setMapQuality( atoi(field_starts[4]));
    const char * cigar = field_starts[5];
    const char * rnext = field_starts[6];
    alignment.setMatePosition( atoi(field_starts[7]));
    alignment.setInsertSize( atoi(field_starts[8]));
    alignment.setQueryBases( field_starts[9]);
    alignment.setQualities(field_starts[10]);

    // zero based indexes:
    alignment.setPosition(alignment.getPosition() - 1);
    alignment.setMatePosition(alignment.getMatePosition() - 1);

    // rname
    if(rname[0] == '*')
        alignment.setRefID(-1);
    else if(header.getSequences().contains(rname))
        alignment.setRefID(header.getSequences().indexOfString(rname));
    else {
        cerr << "Rname " << rname << " missing from sequence dictionary" << endl;
        alignment.setRefID(-1);
        exit(-1);
    }

    // rnext
    if(rnext[0] == '=')
        alignment.setMateRefID(alignment.getRefID());
    else if(rnext[0] == '*')
        alignment.setMateRefID(-1);
    else if(header.getSequences().contains(rnext))
        alignment.setMateRefID(header.getSequences().indexOfString(rnext));
    else {
        cerr << "RNext " << rnext << " missing from sequence dictionary" << endl;
        alignment.setMateRefID(-1);
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
            vector<CigarOp> cigar = alignment.getCigarData();
            cigar.push_back(CigarOp(c,num));
            alignment.setCigarData(cigar);
            cigar_ct ++;
            assert(cigar_ct < 100);
        }
    }
    
    if(!line)
        return al;
    
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
    static const string types("s");
    static const string typec("c");
    static const string typeI("I");
    static const string typeS("S");
    static const string typeC("C");
    //optional attributes
    for (int i = 0; i < additional_parameters.size(); i++) {
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
                if(i <= UCHAR_MAX && i >= 0)
                    retval = alignment.AddTag(tag, typeC, (uint8_t)i);
                else if(i <= SCHAR_MAX && i >= SCHAR_MIN)
                    retval = alignment.AddTag(tag, typec, (int8_t)i);
                else if(i <= USHRT_MAX && i >= 0)
                    retval = alignment.AddTag(tag, typeS, (uint16_t)i);
                else if(i <= SHRT_MAX && i >= SHRT_MIN)
                    retval = alignment.AddTag(tag, types, (int16_t)i);
                else if(i <= UINT_MAX && i >= 0)
                    retval = alignment.AddTag(tag, typeI, (uint32_t)i);
                else if(i <= INT_MAX && i >= INT_MIN)
                    retval = alignment.AddTag(tag, typei, (int32_t)i);
                else
                    cerr << "Error: SAM file contains unparseable int parameter " << tag << ". Skipping. " << endl;
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
                retval = alignment.AddTag<int8_t>(tag, type, *value);
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
    field_starts[0] = NULL;
    
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

bool AddHeaderAttributeArray(OGERead & alignment, const string & tag,const string & value)
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

OGERead * SamReader::read()
{
    return LoadNextAlignment();
}
