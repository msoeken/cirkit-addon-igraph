# CirKit Addon: igraph

## Install and build

Inside the CirKit main directory:

    cd addons
    git clone https://github.com/msoeken/cirkit-addon-igraph.git
    cd ..
    cd build
    cmake -Denable_cirkit-addon-igraph=ON ..
    make cirkit
