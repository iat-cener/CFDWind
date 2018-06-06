FROM unidata/python


USER root

RUN apt-get update


RUN apt-get  install -y --force-yes apt-transport-https
RUN apt-get install -y --force-yes make libpng-dev ssh gnuplot
RUN apt-get update && apt-get install -y --no-install-recommends apt-utils
RUN apt-get install -y --force-yes build-essential binutils-dev cmake flex bison zlib1g-dev qt4-dev-tools libqt4-dev libqtwebkit-dev gnuplot
RUN apt-get install -y --force-yes libreadline-dev libncurses-dev libxt-dev libopenmpi-dev openmpi-bin libboost-system-dev libboost-thread-dev libgmp-dev
RUN apt-get install -y --force-yes libmpfr-dev python python-dev libcgal-dev

USER $NB_UID

RUN mkdir OpenFOAM
WORKDIR OpenFOAM
RUN wget -q "http://downloads.sourceforge.net/foam/OpenFOAM-2.4.0.tgz?use_mirror=mesh" -O OpenFOAM-2.4.0.tgz
RUN wget -q "http://downloads.sourceforge.net/foam/ThirdParty-2.4.0.tgz?use_mirror=mesh" -O ThirdParty-2.4.0.tgz
 
RUN tar -xzf OpenFOAM-2.4.0.tgz 
RUN tar -xzf ThirdParty-2.4.0.tgz

RUN sed -i -e 's/^\(cgal_version=\).*/\1cgal-system/' OpenFOAM-2.4.0/etc/config/CGAL.sh


RUN echo "source \$HOME/OpenFOAM/OpenFOAM-2.4.0/etc/bashrc $FOAM_SETTINGS" >> $HOME/.bashrc     ###This one does nothing usefull in the Dockerfile commands

WORKDIR OpenFOAM-2.4.0
RUN find src applications -name "*.L" -type f | xargs sed -i -e 's=\(YY\_FLEX\_SUBMINOR\_VERSION\)=YY_FLEX_MINOR_VERSION < 6 \&\& \1='

RUN mkdir -p /home/jovyan/OpenFOAM/ThirdParty-2.4.0/platforms/linux64Gcc/cgal-system #mkdir -p $CGAL_ARCH_PATH
RUN mkdir -p /home/jovyan/OpenFOAM/ThirdParty-2.4.0/platforms/linux64Gcc/boost-system #mkdir -p $BOOST_ARCH_PATH

# This next command will take a while... somewhere between 30 minutes to 3 hours.
RUN /bin/bash -c "source $HOME/OpenFOAM/OpenFOAM-2.4.0/etc/bashrc WM_NCOMPPROCS=4 &&\ 
    export QT_SELECT=qt4 &&\
    ./Allwmake > log.make 2>&1 &&\
	./Allwmake" #Run it a second time for getting a summary of the installation


WORKDIR /home/jovyan

USER root
COPY ./applications ./CFDWind/applications
COPY ./exampleCases ./CFDWind/exampleCases
COPY ./src ./CFDWind/src
COPY ./Allwclean ./CFDWind/Allwclean
COPY ./Allwmake ./CFDWind/Allwmake

# Only used for GABLS3 example case
ADD http://hdl.handle.net/11304/59d8d6d7-8300-4aff-aeb4-9812b6153ee9 ./CFDWind/exampleCases/GABLS3/inputData/GABLS3_tendencies_d02_YSU_w60_L9000.nc

RUN chown -R jovyan ./CFDWind 


USER $NB_UID
RUN  mkdir -p ./CFDWind/platforms/linux64GccDPOpt/bin
RUN  mkdir -p ./CFDWind/platforms/linux64GccDPOpt/lib
RUN  mkdir -p ./CFDWind/run
#RUN  mkdir -p ./CFDWind/windmesh/libs
#RUN  mkdir -p ./CFDWind/windmesh/applications

WORKDIR CFDWind

RUN /bin/bash -c "source $HOME/OpenFOAM/OpenFOAM-2.4.0/etc/bashrc WM_NCOMPPROCS=4 &&\ 
    export QT_SELECT=qt4 &&\
	./Allwclean &&\
    ./Allwmake > log.make 2>&1 &&\
	./Allwmake" #Run it a second time for getting a summary of the installation




