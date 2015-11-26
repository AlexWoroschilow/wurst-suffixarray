#$ -clear
#$ -S /bin/sh
#$ -w w
#$ -cwd
#$ -q 4c.q
#$ -j y
# get rid of this #$ -v PATH=$path

./wurstimbiss.x -a alphabet_t500 -s all/pdb_seq.list -b $1 -l 9 -l 15 -n 50 >/scratch/margraf/$1.out

