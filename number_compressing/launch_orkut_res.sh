echo "Taille du graphe: "
read n

result_perf="com-orkut_pagerank_perf.txt"
result_execution="result-orkut_pagerank.txt"
sudo-g5k

for graphTxt in `ls data_com-orkut/*.txt`
do
for s in `seq 2 2`
do
	for k in `seq 20 31`
	do
		let c=k+1
		for rep in `seq 0 2`
		do
			sudo	perf stat -B -e cache-references,cache-misses,L1-dcache-load-misses,L1-dcache-store-misses,L1-dcache-prefetch-misses,LLC-loads,LLC-stores,dTLB-load-misses,dTLB-store-misses,cycles -o $result_perf --append taskset -c 0-$k ./pagerank_pthreads $graphTxt $n 0.00000000001 0.85 $c 7 1 0 $s >>$result_execution
		done
	done
done
done
