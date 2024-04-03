export CPUS_PER_TASK=16
export N_PARAMS=7
echo "Running on ${N_PARAMS} nodes with ${CPUS_PER_TASK} cores per node."
echo ""

mkdir logs

echo "Submit results job."
jid1=$(sbatch  --array=1-$N_PARAMS -c $CPUS_PER_TASK --export=NUM_THREADS=$CPUS_PER_TASK --parsable ./get_results.sh)
echo "DONE."
echo ""

echo "Submit output job"
jid2=$(sbatch --time=1:30:00 --depend=afterany:$jid1 ./make_output.sh)
#jid2=$(sbatch --time=1:30:00 ./make_output.sh)
echo "DONE."
echo ""
