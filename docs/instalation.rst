

Installation/Compiling
----------------------

The suggested way for using CFDWind is by using Docker https://www.docker.com/
Docker works under most of the operating systems and allows deplying CFDWind from a single packadge.


Quick start using Docker Hub
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

A pre-build image of CFDWind is available from the Docker Hub. Once you have the Docker installed ( https://store.docker.com/search?type=edition&offering=community ) you can run the CFDWind environment with a single command::

    $ docker run -it windbench/cfdwind /bin/bash

The first time the command might take some time downloading the image.

The above command will open a terminal that is ready to use.
You can test it (takes 1 minute) by running::

    $ ./test/pre-run-quickTest.sh


.. todo:: To share input/output with the OS, and avoid storing everything in the virtual drive, we should mount workdir to the host. 


Development workflow
^^^^^^^^^^^^^^^^^^^^

1. You will need:
	- Docker - https://store.docker.com/search?type=edition&offering=community
	- and git - https://gist.github.com/derhuerst/1b15ff4652a867391f03

2. Clone the CFDWind repository::

    $ git clone https://github.com/iat-cener/CFDWind.git cfdwind

3. Build your local CFDWind Docker image (for the first time it will take around 2 hours)::

    $ cd cfdwind
    $ docker build -t cfdwind ./

At this point all is ready for development

4. Make some changes in the code
5. Rebuild the image (it should take only few minutes)::

    $ docker build -t cfdwind ./

.. todo:: The build step rebuilds the whole model. Better if recompiled inside the image.

6. Run a quick model test (takes around 1 minute)::

    $ docker run -it cfdwind /bin/bash -c "source ~/OpenFOAM/OpenFOAM-2.4.0/etc/bashrc && ./test/pre-run-quickTest.sh"

7. Get involved, improove the code and repeat steps 4-6.


Without Docker
^^^^^^^^^^^^^^

For a custom instalation, it is easiest to follow the build script for the Docker image https://github.com/iat-cener/CFDWind/blob/master/Dockerfile which lists all the dependencies.

