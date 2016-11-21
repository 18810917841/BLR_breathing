########################################################
#PBS -N run3.sh
#PBS -o /ihepbatch/mbhd01/user/liyanrong/Yuyang/RAD_V1/mpites.out
#PBS -e /ihepbatch/mbhd01/user/liyanrong/Yuyang/RAD_V1/mpites.err
#PBS -l nodes=1:ppn=16


DIR_HERE=/ihepbatch/mbhd01/user/liyanrong/Yuyang/RAD_V1

cat $PBS_NODEFILE > $DIR_HERE/hostfile
NCPU=`wc -l < $PBS_NODEFILE`
NNODES=`uniq $PBS_NODEFILE | wc -l`
HOST=$(uniq $PBS_NODEFILE)

echo "$HOST" >> $DIR_HERE/hostfile

cd $DIR_HERE

mpiexec -f ${PBS_NODEFILE} -np $NCPU $DIR_HERE/Radiation
#########################################################
