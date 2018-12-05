# .bashrc

# Source global definitions
if [ -f /etc/bashrc ]; then
	. /etc/bashrc
fi

# User specific aliases and functions
alias tmp='cd /tmpdir/blaised'
alias t1='cd /tmpdir/blaised/to_transfer'
alias t2='cd /tmpdir/blaised/to_transfer2'
alias t3='cd /tmpdir/blaised/to_transfer3'
alias t4='cd /tmpdir/blaised/to_transfer4'
alias t5='cd /tmpdir/blaised/to_transfer5'
alias t6='cd /tmpdir/blaised/to_transfer6'
alias t7='cd /tmpdir/blaised/to_transfer7'
alias t8='cd /tmpdir/blaised/to_transfer8'
alias t9='cd /tmpdir/blaised/to_transfer9'
alias t10='cd /tmpdir/blaised/to_transfer10'

alias solclean='rm *.bin *.end'
alias fluclean='rm uf* vf* wf*'
alias alloc='salloc -N 1 -n 20'
alias qu='squeue -u $USER'
alias cnode='ssh $SLURM_NODELIST'
alias sub='sbatch script.slurm'
