result_execution="result-orkut_pagerank.txt"

for s in `seq 0 2`
do
	for k in `seq 15 15`
	do
		let c=k+1
		for rep in `seq 0 2`
		do
			taskset -c 0-$k ./pagerank_pthreads data/com-orkut_prepared_cn_naive.txt 3072441 0.00000000001 0.85 $c 1 1 0 $s >> $result_execution
		done
	done
done

for s in `seq 0 2`
do
	for k in `seq 15 15`
	do
		let c=k+1
		for rep in `seq 0 2`
		do
			taskset -c 0-$k ./pagerank_pthreads data/com-orkut_prepared.txt 3072441 0.00000000001 0.85 $c 1 1 0 $s >> $result_execution
		done
	done
done
