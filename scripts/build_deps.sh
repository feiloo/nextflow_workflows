set -euo pipefail

# install libbz2-devel gsl-devel nasm yasm

mkdir -p deps
export GRADLE_USER_HOME=$(pwd)/deps/gradle
cd deps

###buildstep
if [ ! -d gatk ] ; then
git clone --recursive --depth 1 https://github.com/broadinstitute/gatk.git
fi
cd gatk
#git pull
#git checkout 4.5.0.0

./gradlew localJar
#cp -r libs ..
#cp gatk ..
#./gradlew clean
#rm -rf src
#rm -rf ~/.gradle

cd ..

###buildstep
if [ ! -d htslib ] ; then
git clone --recursive https://github.com/samtools/htslib.git
fi

cd htslib
autoreconf -i
./configure
echo $(pwd)
make -j

cd ..

###buildstep
if [ ! -d samtools ] ; then
git clone --recursive https://github.com/samtools/samtools.git
fi

cd samtools
autoheader
autoconf -Wno-syntax
./configure --without-curses
make -j


cd ..

###buildstep
git clone --recursive https://github.com/bwa-mem2/bwa-mem2
cd bwa-mem2

make -j
cd ..

###buildstep
# assumes we already have htslib
git clone --recursive https://github.com/samtools/bcftools.git
cd bcftools
autoheader && autoconf && ./configure --enable-libgsl --enable-perl-filters
make -j

cd ..

###buildstep
# fastp

## libisal
git clone --recursive https://github.com/intel/isa-l.git
cd isa-l
./autogen.sh
./configure --prefix=/usr --libdir=/usr/lib64
make
cd ..

###buildstep
# libdeflate
git clone https://github.com/ebiggers/libdeflate.git
cd libdeflate
cmake -B build
cmake --build build
cmake --install build
cd ..

###buildstep
# fastp
git clone https://github.com/OpenGene/fastp.git
cd fastp
