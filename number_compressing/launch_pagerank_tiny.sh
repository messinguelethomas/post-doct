echo "Taille du graphe: "
read n

result_perf="com-lj_pagerank_perf.txt"


for graphTxt in `ls data_com-lj/*.txt`
do
	for s in `seq 0 2`
	do
			for i in `seq 1 2`
			do
		sudo	perf stat -B -e cache-references,cache-misses,L1-dcache-load-misses,L1-dcache-store-misses,L1-dcache-prefetch-misses,LLC-loads,LLC-stores,dTLB-load-misses,dTLB-store-misses,cycles -o $result_perf --append taskset -c 0-0 ./pagerank_pthreads $graphTxt $n 0.00000000001 0.85 1 1 1 0 $s
			done
	done
done
