#/bin/sh

# SGE: the job name
#$ -N peeken_1part
#
#$ -V
#
#$ -q all.q@science-bs35      # dont make it longer than 24 hours
#
# The requested run-time, expressed as (xxxx sec or hh:mm:ss)
#$ -l s_rt=160:02:00
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
python advect_ice_peeken.py
