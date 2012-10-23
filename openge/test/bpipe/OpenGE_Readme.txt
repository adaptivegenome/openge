Bpipe unit tests - OpenGE

OpenGE uses Bpipe's unit tests in order to (1) save work and to (2) ensure compatibility with Bpipe when possible. Some changes have been made to tests in order to reflect features that we don't plan on supporting, or to improve tests. These changes are outlined here, so in the future when importing changes to the Bpipe tests, we can ensure that we don't overwrite these changes.

* slow test: decreased wait time from 60 seconds to 6 seconds (so that test execution doesn't take as long)
* filter test: changed Java annotation style to function style (annotation style isn't supported in OpenGE) 
* produce test: changed Java annotation style to function style (annotation style isn't supported in OpenGE) 
* transform test: changed Java annotation style to function style (annotation style isn't supported in OpenGE)

Other tests test features that we will have no plans of supporting- for example, inline Java/Groovy code:

* nested_parallel: inline Java code
* produce_glob: inline Java code
* rscript: inline R code

In many cases, these tests should be rewritten or replaced.
