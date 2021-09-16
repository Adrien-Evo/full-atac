#rm -r out 
#rm snakejob*
snakemake --snakefile snakefile.smk --use-conda --cluster "sleep 1s && qsub -q max-7d.q -e ./logs/ -o ./logs/" --jobs 50 --jobscript job.sh -pn
