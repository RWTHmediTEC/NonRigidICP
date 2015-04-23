According to the profiler, 'pointTriangleDistance' is the computationally most
expensive function. You can speed it up by a factor of 10 by transforming it
into a mex file. This can be done automatically using 'coder'. The project file
has already been created. Thus it suffice to just run

>> coder -build pointToTriangleDistance.prj