
#ifndef OpenGE_AlignmentBlock_h
#define OpenGE_AlignmentBlock_h
/*********************************************************************
 *
 * AlignmentBlock.cpp: Port of GATK's AlignmentBlock.
 * Open Genomics Engine
 *
 * Author: Lee C. Baker, VBI
 * Last modified: 9 May 2012
 *
 *********************************************************************
 *
 * This file has been ported from GATK's implementation in Java, and
 * is released under the Virginia Tech Non-Commercial Purpose License.
 * A copy of this license has been provided in  the openge/ directory.
 * 
 * The original file, AlignmentBlock.java, was released 
 * under the following license:
 *
 * Copyright (c) 2009 The Broad Institute. 
 * Ported to C++ by Lee C. Baker, Virginia Bioinformatics Institute
 *
 * Permission is hereby granted, free of charge, to any person
 * obtaining a copy of this software and associated documentation
 * files (the "Software"), to deal in the Software without
 * restriction, including without limitation the rights to use,
 * copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the
 * Software is furnished to do so, subject to the following
 * conditions:
 *
 * The above copyright notice and this permission notice shall be
 * included in all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
 * EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
 * OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
 * NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
 * HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
 * WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR
 * THE USE OR OTHER DEALINGS IN THE SOFTWARE.
 */
//package net.sf.samtools;

/**
 * Represents the contiguous alignment of a subset of read bases to a reference
 * sequence. Simply put an alignment block tells you that read bases from
 * readStart are aligned to the reference (matching or mismatching) from
 * referenceStart for length bases.
 *
 * @author Tim Fennell
 */
class AlignmentBlock {
private:
    int readStart;
    int referenceStart;
    int length;
    
public:
    /** Constructs a new alignment block with the supplied read and ref starts and length. */
    AlignmentBlock(int readStart, int referenceStart, int length)
    : readStart(readStart)
    , referenceStart(referenceStart)
    , length(length)
    { }
    
    /** The first, 1-based, base in the read that is aligned to the reference reference. */
    int getReadStart() { return readStart; }
    
    /** The first, 1-based, position in the reference to which the read is aligned. */
    int getReferenceStart() { return referenceStart; }
    
    /** The number of contiguous bases aligned to the reference. */
    int getLength() { return length; }
};

#endif
