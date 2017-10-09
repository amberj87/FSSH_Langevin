SCR_DIR=/scratch/$PBS_JOBID
mkdir $SCR_DIR
cd $PBS_O_WORKDIR

cp aout $SCR_DIR
cp *.inp $SCR_DIR

cd $SCR_DIR
./aout

cp output* $PBS_O_WORKDIR
cp average_result* $PBS_O_WORKDIR
cp fort.* $PBS_O_WORKDIR
cp *.out $PBS_O_WORKDIR

cd ..
rm -r $SCR_DIR

