//
//  SamReader.cpp
//  BamTools
//
//  Created by Lee Baker on 3/16/12.
//  Copyright (c) 2012 LCB. All rights reserved.
//

#include <iostream>
#include <sstream>
#include "SamReader.h"
#include <cstring>
using namespace std;

SamReader::SamReader()
: file(NULL)
, loaded(false)
{ }

bool SamReader::Open(const string & filename)
{
  file.open(filename.c_str(), ios::in);

  LoadHeaderData();

  loaded = true;

  return true;
}

bool SamReader::Close()
{
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
bool SamReader::LoadNextAlignment(BamAlignment& alignment)
{
  if(!loaded)
    return false;

  if(file.eof())
    return false;

  string line_s;
  getline(file, line_s);
    


  string & qname = alignment.Name;
  uint32_t & flag = alignment.AlignmentFlag;
  string rname;
  int & pos = alignment.Position;
  uint16_t & mapq = alignment.MapQuality;
  string cigar;
  string rnext;
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
    rname.assign(field_starts[2]);
    pos = atoi(field_starts[3]);
    mapq = atoi(field_starts[4]);
    cigar.assign(field_starts[5]);
    rnext.assign(field_starts[6]);
    pnext = atoi(field_starts[7]);
    tlen = atoi(field_starts[8]);
    seq.assign(field_starts[9]);
    qual.assign(field_starts[10]);
    
    // zero based indexes:
    alignment.Position--;
    alignment.MatePosition--;
    
    // rname
    if(rname == "*")
        alignment.RefID = -1;
    else if(header.Sequences.Contains(rname))
        alignment.RefID = header.Sequences.IndexOfString(rname);
    else {
        cerr << "Rname " << rname << " missing from sequence dictionary" << endl;
        alignment.RefID = -1;
    }
    
    // rnext
    if(rnext == "=")
        alignment.MateRefID = alignment.RefID;
    else if(rnext == "*")
        alignment.MateRefID = -1;
    else if(header.Sequences.Contains(rnext))
        alignment.MateRefID = header.Sequences.IndexOfString(rnext);
    else {
        cerr << "RNext " << rnext << " missing from sequence dictionary" << endl;
        alignment.MateRefID = -1;
    }
    
    //CIGAR ops
    if(cigar.c_str()[0] != '*')
    {
        int cigar_ct = 0;
        char * cigar_p = (char *) cigar.c_str();    //cast away const-ness, because we will be changing the ptr (not the data)
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
        return true;
    
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

    //optional attributes
    for (int i = 1; i < additional_parameters.size(); i++) {
        char * segment = additional_parameters[i];

        //null terminate field separators (:) in the segment so we can use to construct strings
        segment[2] = 0;
        segment[4] = 0;

        string tag(segment);
        string type(&segment[3]);
        string value(&segment[5]);
        
        bool retval = false;
        switch(type[0])
        {
            case 'i': //int
            {
                int i = atoi(value.c_str());
                retval = alignment.AddTag(tag, type, i);
                break;
            }
            case 'f': //floating point
            {
                float f = atof(value.c_str());
                retval = alignment.AddTag(tag, type, f);
                break;
            }
            case 'A': //single char
            {
                char a = value.c_str()[0];
                retval = alignment.AddTag(tag, type, a);
                break;
            }
            case 'Z': //string, with spaces
                retval = alignment.AddTag(tag, type, value);
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

  return true;
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
  return LoadNextAlignment(alignment);
}
