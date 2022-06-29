source activate /home/osmalama/softwares/conda_envs/MOODS

export PATH=/home/osmalama/softwares/MOODS-python-1.9.4.1/scripts:$PATH

M are PWM matrices

S is the sequence

n=300000 the number of best hits

moods-dna.py -S M1 M2 -s S -B $n -o outfile
