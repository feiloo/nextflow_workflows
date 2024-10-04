set -euo pipefail

cd
cd gatk
cp $(pwd)/build/libs/* /usr/local/lib
export GATK_LOCAL_JAR=/usr/local/lib/gatk.jar
echo "export GATK_LOCAL_JAR=/usr/local/lib/gatk.jar" >> /etc/bash.bashrc.local
cp gatk /usr/local/bin
cd ..

cd htslib
make install
cd ..

cd samtools
make install
cd ..

cd bwa-mem2
cp bwa-mem2 bwa-mem2.avx bwa-mem2.avx2 bwa-mem2.avx512bw bwa-mem2.sse41 bwa-mem2.sse42 /usr/local/bin
cd ..

cd bcftools
make install
cd ..

cd isa-l
make install
cd ..

cd libdeflate
cmake --install build
cd ..

cd fastp
#export LD_LIBRARY_PATH=$LD_LIBRARY_PATH
#echo "export LD_LIBRARY_PATH=$LD_LIBRARY_PATH" >> /etc/bash.bashrc.local
make install
cd ..



