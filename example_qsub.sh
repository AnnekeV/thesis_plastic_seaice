#/bin/sh

# SGE: the job name
#$ -N particles_everywhere_trial_less_particles
#
# The requested run-time, expressed as (xxxx sec or hh:mm:ss)
#$ -l s_rt=10:02:00
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
python advect_glorys_data-all.py
