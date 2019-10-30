#/bin/sh

# SGE: the job name
#$ -N  plot_drifter
#
#$ -l hostname=science-bs35 #make sure it is in the right node
#$ -V
#
# SGE: your Email here, for job notification
#$ -M a.l.vries@students.uu.nl
#
# SGE: when do you want to be notified (b : begin, e : end, s : error)?
#$ -m e 
#$ -m s
#
# SGE: ouput in the current working dir
#$ -wd /home/students/6252699/thesis/working_directory   
#
cd /home/students/6252699/thesis/parcels2
python plot_drifter.py
