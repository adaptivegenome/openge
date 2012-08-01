#!/bin/bash
# This file takes two parameters- a filename, and a command to run. Since CTest can't use pipes, we need this to test piping abilities of OpenGE.

eval cat $1 | $2
retval=$?

exit $retval
