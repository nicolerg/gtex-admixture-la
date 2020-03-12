for chr in {1..22}
do
	taskset -c 36-46 Rscript /users/nicolerg/gtex-admix/scripts/f2-analysis/gtex-f2.R ${chr} &

	running=`ps -ef | grep "nicolerg" | grep "gtex-f2.R" | wc -l`
	while [ $running -gt 5 ] # run max 5 jobs at a time
	do
		sleep 180
		running=`ps -ef | grep "nicolerg" | grep "gtex-f2.R" | wc -l`
	done
done
