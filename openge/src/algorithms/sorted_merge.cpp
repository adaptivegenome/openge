/*********************************************************************
 *
 * sorted_merge.cpp:  Merge multiple read streams, keeping them 
 *                    sorted.
 * Open Genomics Engine
 *
 * Author: Lee C. Baker, VBI
 * Last modified: 7 June 2012
 *
 *********************************************************************
 *
 * This file is released under the Virginia Tech Non-Commercial 
 * Purpose License. A copy of this license has been provided in 
 * the openge/ directory.
 *
 *********************************************************************/

#include "sorted_merge.h"

#include <algorithm>
#include <numeric>
#include <pthread.h>

#include "../util/bamtools/Sort.h"

using namespace std;

SortedMerge::~SortedMerge()
{
    for( int i = 0; i < input_proxies.size(); i++)
        delete input_proxies[i];
}

SortedMerge::SortedMergeInputProxy::SortedMergeInputProxy(SortedMerge * parent) 
: merge_class(parent)
{ }

int SortedMerge::SortedMergeInputProxy::runInternal()
{
    merge_class->done_signal_mutex.lock();
    while(!merge_class->done.isSet())
        merge_class->done_signal.wait(merge_class->done_signal_mutex);
    merge_class->done_signal_mutex.unlock();

    return 0;
}

void SortedMerge::addSource(AlgorithmModule * source)
{
    SortedMergeInputProxy * proxy = new SortedMergeInputProxy(this);
    input_proxies.push_back(proxy);
    
    //set parent to be the first source, so parsing around the tree keeps working.
    if(input_proxies.size() == 1)
        proxy->addSink(this);

    source->addSink(proxy);
}

bool SortedMerge::SortedMergeElement::operator<(const SortedMergeElement & t) const
{
    BamTools::Algorithms::Sort::ByPosition cmp = BamTools::Algorithms::Sort::ByPosition();
    return cmp(this->read, t.read);
}

int SortedMerge::runInternal()
{
    ogeNameThread("am_merge_sorted");

    multiset<SortedMergeElement> reads;
    
    // first, get one read from each queue
    // make sure and deal with the case where one chain will never have any reads. TODO LCB

    for(int ctr = 0; ctr < input_proxies.size(); ctr++)
    {
        OGERead * read = input_proxies[ctr]->getInputAlignment();

        if(!read) {
            done_signal_mutex.lock();
            done_signal.notify_all();
            done_signal_mutex.unlock();
            continue;
        }

        reads.insert(SortedMergeElement(read, input_proxies[ctr]));
    }

    //now handle the steady state situation. When sources are done, We
    // won't have a read any more in the reads pqueue.
    while(!reads.empty()) {
        SortedMergeElement el = *reads.begin();
        reads.erase(reads.begin());

        putOutputAlignment(el.read);

        el.read = el.source->getInputAlignment();
        if(!el.read) {
            done_signal_mutex.lock();
            done_signal.notify_all();
            done_signal_mutex.unlock();
            continue;
        }

        reads.insert(el);
    }

    return 0;
}