Open Genomics Engine

README
Bamtools combined merge/sort command, enhancement 1

Author: Lee Baker, lee@leecbaker.com for David Mittelman
22 March 2012

1. Description
This patch builds on the previously sent 'merge sort' patch, with three significant changes:
* The -forceCompression and -forceNoCompression flags for the mergesort command have been replaced with a single, optional -compressionLevel <#> flag, allowing selection of the level of compression of the saved BAM file.
* A new parameter, -threads <#>, sets the number of threads in each thread pool. This parameter is optional; when not present, the number of threads is set equal to the number of cores.
* A new command has been added to convert sam files to bam files called "compress".

2. Usage
For the merge sort command, options are as described in the previous README, with exceptions outlined below:

-threads <num_threads>: Specifies the number of threads to use in each thread pool. When omitted, the number is automatically set equal to the number of processing cores, as reported by the OS. Note that setting this to 1 does not disable multithreading, but only creates one worker thread per pool. Use -noThreads to disable thread pools completely.

-compressionLevel <level> Specifies the level of compression to apply the the BAM file being saved. 0 corresponds to no compression, and 9 to most compression. The default level  is 6. See the zlib documentation for more details on the compression level.

Run "openge help mergesort"  for a brief description of each command.

The 'compress' command takes in a single SAM file, and produces a single BAM file containing the same information. The following parameters may be specified:

-in <filename>: The sam file to read from. When not specified, stdin is used.
-out <filename>: The location to write the bam file. When omitted, stdout is used.

3. How to use patch

The following commands will download a fresh copy of bamtools, apply the necessary patch, and compile bamtools. More details can be found in the bamtools documentation at:
https://github.com/pezmaster31/bamtools

As I update the patch, the filename is changed to reflect a version number (eg. mergesort_009.diff). The patch can be found in the same zip file as this README.

mkdir openge
git clone git://github.com/pezmaster31/bamtools.git
mkdir openge/build
cd openge
git apply patch_025.diff
cd build
cmake ..
make

4. Example commands:

To convert a SAM file to a BAM file:
openge compress -in a.sam -out a.bam

To convert a BAM file to a SAM file (already existed in code):
openge convert -format sam -in a.bam -out a.sam

To merge and sort bam files, saving the result at a lower level of compression than the default:
openge mergesort -v -compressionLevel 3  -in a.bam -in b.bam -in c.bam -out d.bam

To merge and sort bam files, using 12 cores:
openge mergesort -v -threads 12 -in a.bam -in b.bam -in c.bam -out d.bam

Additional information about commands and what they do can be found in the more complete PDF documentation included with this README.
