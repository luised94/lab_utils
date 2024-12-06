# Update and install necessary dependencies
sudo apt-get update
sudo apt-get install -y build-essential libcurl4-gnutls-dev libxml2-dev libssl-dev gfortran
sudo apt-get install -y libx11-dev libxt-dev libpng-dev libjpeg-dev libcairo2-dev libxext-dev libxrender-dev libxmu-dev libxmuu-dev x11-apps xauth libreadline-dev libbz2-dev liblzma-dev
sudo apt-get install -y default-jdk
# Install related libraries for R packages.
sudo apt-get install libharfbuzz-dev libfribidi-dev libfreetype6-dev libpng-dev libtiff5-dev libfontconfig1-dev
sudo apt-get install libmagick++-dev

# Install LaTeX and related packages
sudo apt-get install -y texlive texlive-fonts-extra texlive-latex-extra texinfo
sudo apt-get install -y texlive-science texlive-extra-utils texlive-bibtex-extra
sudo apt-get install -y texlive-fonts-recommended texlive-plain-generic
sudo apt-get install -y texlive-fonts-extra

# Set up environment variables
export JAVA_HOME=/usr/lib/jvm/java-11-openjdk-amd64
export PATH=$JAVA_HOME/bin:$PATH

# Download and extract R source
R_VERSION="R-4.2.0"
DIR_TO_INSTALL=$HOME
curl --output "$HOME/${R_VERSION}.tar.gz" "https://cran.r-project.org/src/base/$(echo ${R_VERSION} | cut  -d. -f1)/${R_VERSION}.tar.gz"
tar -xzvf ${R_VERSION}.tar.gz
cd ${R_VERSION}

# Configure and compile R with X11 and Cairo support
./configure --enable-R-shlib --with-blas --with-lapack --with-cairo --with-x
make
sudo make install

# Set up environment
export PATH=/usr/local/bin:$PATH
mkdir -p $HOME/R/library
export R_LIBS_USER=$HOME/R/library
alias R='R --no-save'

# Start R and install packages
R --vanilla << EOF
dir.create(Sys.getenv("R_LIBS_USER"), recursive = TRUE)
.libPaths(Sys.getenv("R_LIBS_USER"))
options(repos = c(CRAN = "https://cloud.r-project.org"))
install.packages(c("renv", "xml2", "lintr", "roxygen2", "languageserver"), dependencies = TRUE, INSTALL_opts = '--no-lock')
q()
EOF
