FROM centos:7
RUN rpm --import https://repo.anaconda.com/pkgs/misc/gpgkeys/anaconda.asc
RUN mkdir /etc/yum/repos.d
RUN echo -e "[conda]\nname=Conda\nbaseurl=https://repo.anaconda.com/pkgs/misc/rpmrepo/conda\nenabled=1\ngpgcheck=1\ngpgkey=https://repo.anaconda.com/pkgs/misc/gpgkeys/anaconda.asc" > /etc/yum/repos.d/conda.repo
RUN yum -y update
RUN yum -y install conda
RUN yum -y install less
RUN source /opt/conda/etc/profile.d/conda.sh; conda create -y --name AmberTools20; conda activate AmberTools20; conda install -y -c conda-forge ambertools=20
RUN mkdir /deploy
RUN mkdir /deploy/utils /deploy/artifacts /deploy/src /deploy/lib /deploy/config
RUN mkdir /deploy/artifacts/base /deploy/artifacts/usr
RUN mkdir /deploy/src/pyAOMlite /deploy/src/STO-proj-AOM-overlap /deploy/src/fobsh-cp2k
RUN yum -y install nano
RUN yum group install -y "Development Tools"
RUN yum -y install lapack-devel
RUN yum -y install fftw-devel
ARG CP2K_VERSION
ENV CP2K_VERSION=${CP2K_VERSION}
COPY src_cp2k/${CP2K_VERSION}.zip /deploy/artifacts/base
COPY config/cp2k.local.sopt_centos_libs /deploy/config/cp2k.local.sopt
RUN unzip /deploy/artifacts/base/${CP2K_VERSION}.zip -d /deploy/src/fobsh-cp2k/
RUN cp /deploy/config/cp2k.local.sopt /deploy/src/fobsh-cp2k/${CP2K_VERSION}/cp2k/arch/local.sopt
RUN sed -i -e 's/python/python2/g' /deploy/src/fobsh-cp2k/${CP2K_VERSION}/cp2k/tools/build_utils/*.py
RUN cd /deploy/src/fobsh-cp2k/${CP2K_VERSION}/cp2k/makefiles/ ; make ARCH=local VERSION=sopt
ENV CP2K_EXEC=/deploy/src/fobsh-cp2k/${CP2K_VERSION}/cp2k/exe/local/cp2k.sopt
COPY utils/*.sh /deploy/utils/
COPY utils/*.c /deploy/utils/
COPY utils/*.py /deploy/utils/
RUN chmod +x /deploy/utils/*.sh ; chmod +x /deploy/utils/*.py
RUN /deploy/utils/build_C.sh
RUN yum install -y wget
RUN wget https://github.com/oziogos/pyAOMlite/archive/refs/heads/main.zip ; mv main.zip /deploy/artifacts/base/pyAOMlite.zip ; 
RUN unzip /deploy/artifacts/base/pyAOMlite.zip -d /deploy/src/pyAOMlite/
# HAB79 test
RUN mkdir /work/ ; mkdir /work/HAB79 ; mkdir /work/HAB79/single_molecules
COPY test/HAB79/single_molecules/*.mol2 /work/HAB79/single_molecules/
COPY test/HAB79/README /work/HAB79/
COPY test/HAB79/test.sh /work/HAB79/
RUN chmod +x /work/HAB79/test.sh
COPY test/fobsh_input/*template* /deploy/config/
# start!
CMD source /opt/conda/etc/profile.d/conda.sh; conda activate AmberTools20; bash
