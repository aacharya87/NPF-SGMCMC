FROM continuumio/miniconda

RUN apt-get update
RUN apt-get -y install \
    build-essential \
    libgsl0ldbl \
    gsl-bin \
    git \
    libgsl0-dev \
    libblas-dev \
    liblapack-dev \
    libatlas-base-dev \
    gfortran \
    cmake \
    libarmadillo-dev
RUN conda update conda
RUN conda install jupyter -y --quiet \
    && mkdir /opt/notebooks
RUN conda install pip pandas scipy matplotlib seaborn pymongo ipython scikit-learn nltk gensim
RUN pip install lda pyldavis
RUN apt-get -y install libboost-all-dev nano vim
RUN apt -y install libgl1-mesa-glx
RUN ["apt-get", "update"]
RUN ["apt-get", "install", "-y", "zsh"]
RUN wget https://github.com/robbyrussell/oh-my-zsh/raw/master/tools/install.sh -O - | zsh || true
