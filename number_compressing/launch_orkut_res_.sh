result_execution="result-lj_pagerank.txt"

	for k in `seq 0 2`
	do
		for n in `seq 0 2`
		do
			taskset -c 0-31 ./pagerank_pthreads data/com-orkut_prepared_gorder.txt 3072441 0.00000000001 0.85 32 1 1 0 $k >> $result_execution
		done
	done
	for k in `seq 0 2`
	do
		for n in `seq 0 2`
		do
			taskset -c 0-31 ./pagerank_pthreads data/com-orkut_prepared_rabbit.txt 3072441 0.00000000001 0.85 32 1 1 0 $k >> $result_execution
		done
	done
	for k in `seq 0 2`
	do
		for n in `seq 0 2`
		do
			taskset -c 0-31 ./pagerank_pthreads data/com-orkut_prepared_cn_rabbit.txt 3072441 0.00000000001 0.85 32 1 1 0 $k >> $result_execution
		done
	done
	for k in `seq 0 2`
	do
		for n in `seq 0 2`
		do
			taskset -c 0-31 ./pagerank_pthreads data/com-orkut_prepared_cn_louvain.txt 3072441 0.00000000001 0.85 32 1 1 0 $k >> $result_execution
		done
	done
	for k in `seq 0 2`
	do
		for n in `seq 0 2`
		do
			taskset -c 0-31 ./pagerank_pthreads data/com-orkut_prepared_numbaco.txt 3072441 0.00000000001 0.85 32 1 1 0 $k >> $result_execution
		done
	done
