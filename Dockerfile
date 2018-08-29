FROM unidata/python


USER root


RUN apt-get update && \
    apt-get install -y --no-install-recommends apt-transport-https make libpng-dev ssh gnuplot apt-utils && \
	apt-get install -y --no-install-recommends build-essential binutils-dev cmake flex bison zlib1g-dev qt4-dev-tools libqt4-dev libqtwebkit-dev gnuplot && \
	apt-get install -y --no-install-recommends libreadline-dev libncurses-dev libxt-dev libopenmpi-dev openmpi-bin libboost-system-dev libboost-thread-dev libgmp-dev && \
	apt-get install -y --no-install-recommends libmpfr-dev python python-dev libcgal-dev && \
    apt-get clean && \
	rm -rf /var/lib/apt/lists/*

USER $NB_UID

RUN mkdir OpenFOAM
WORKDIR OpenFOAM
RUN wget -q "http://downloads.sourceforge.net/foam/OpenFOAM-2.4.0.tgz?use_mirror=mesh" -O OpenFOAM-2.4.0.tgz  && \
	wget -q "http://downloads.sourceforge.net/foam/ThirdParty-2.4.0.tgz?use_mirror=mesh" -O ThirdParty-2.4.0.tgz  && \
	tar -xzf OpenFOAM-2.4.0.tgz && \
	tar -xzf ThirdParty-2.4.0.tgz  && \
	sed -i -e 's/^\(cgal_version=\).*/\1cgal-system/' OpenFOAM-2.4.0/etc/config/CGAL.sh && \
	echo "source \$HOME/OpenFOAM/OpenFOAM-2.4.0/etc/bashrc $FOAM_SETTINGS" >> $HOME/.bashrc && \
	find ./OpenFOAM-2.4.0/src ./OpenFOAM-2.4.0/applications -name "*.L" -type f | xargs sed -i -e 's=\(YY\_FLEX\_SUBMINOR\_VERSION\)=YY_FLEX_MINOR_VERSION < 6 \&\& \1=' && \
	mkdir -p /home/jovyan/OpenFOAM/ThirdParty-2.4.0/platforms/linux64Gcc/cgal-system  && \
	mkdir -p /home/jovyan/OpenFOAM/ThirdParty-2.4.0/platforms/linux64Gcc/boost-system && \
	/bin/bash -c "cd OpenFOAM-2.4.0 &&\
		source $HOME/OpenFOAM/OpenFOAM-2.4.0/etc/bashrc WM_NCOMPPROCS=4 &&\ 
	    export QT_SELECT=qt4 &&\
	    ./Allwmake > log.make 2>&1 &&\
		./Allwmake" &&\
	rm OpenFOAM-2.4.0.tgz && \
	rm ThirdParty-2.4.0.tgz &&\
	rm -R ./OpenFOAM-2.4.0/doc &&\
	rm -R ./OpenFOAM-2.4.0/tutorials &&\
	rm -R ./ThirdParty-2.4.0/cmake-2.8.12.1 &&\
	rm -R ./ThirdParty-2.4.0/ParaView-4.1.0 &&\
	rm -R ./ThirdParty-2.4.0/scotch_6.0.3 &&\
	rm -R ./ThirdParty-2.4.0/CGAL-4.6 &&\
	rm -R ./ThirdParty-2.4.0/etc &&\
	rm -R ./ThirdParty-2.4.0/openmpi-1.8.5 &&\
	rm ./ThirdParty-2.4.0/* || true

WORKDIR /home/jovyan

USER root
COPY ./applications ./CFDWind/applications
COPY ./exampleCases ./CFDWind/exampleCases
COPY ./src ./CFDWind/src
COPY ./test ./CFDWind/test
COPY Allwclean Allwmake CFDWind/


# Only used for GABLS3 example case
RUN mkdir -p ./CFDWind/exampleCases/GABLS3/inputData  &&\
	wget -q "http://hdl.handle.net/11304/59d8d6d7-8300-4aff-aeb4-9812b6153ee9" -O ./CFDWind/exampleCases/GABLS3/inputData/GABLS3_tendencies_d02_YSU_w60_L9000.nc

RUN chown -R jovyan ./CFDWind 


USER $NB_UID
RUN mkdir -p ./CFDWind/platforms/linux64GccDPOpt/bin  &&\
	mkdir -p ./CFDWind/platforms/linux64GccDPOpt/lib  &&\
	mkdir -p ./CFDWind/run


WORKDIR CFDWind

RUN /bin/bash -c "source $HOME/OpenFOAM/OpenFOAM-2.4.0/etc/bashrc WM_NCOMPPROCS=4 &&\ 
    	export QT_SELECT=qt4 &&\
		./Allwclean &&\
	    ./Allwmake > log.make 2>&1 &&\
		./Allwmake" 
# &&\
		#rm -R $HOME/OpenFOAM/OpenFOAM-2.4.0/src


