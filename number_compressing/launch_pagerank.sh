echo "Nombre de coeur : "
read nbcoeur
echo "Taille du graphe: "
read n

result_perf="result_perf.txt"
filename="graph_file_orkut.txt"
graph_file=`cat $filename`
#echo Start
for graphTxt in `cat $filename`
do
	for scd in `seq 1 12`
	do
		#if [ $scd -eq 1 ] || [ $scd -eq 6 ] || [ $scd -eq 7 ] || [ $scd -eq 12 ]
		if [ $scd -eq 1 ] || [ $scd -eq 7 ]
		then
			for k in `seq 0 $nbcoeur`
			do
				if [ $k -eq 9 ] || [ $k -eq 19 ] || [ $k -eq 31 ]
				then
					let c=k+1
					for i in `seq 1 2`
					do
					perf stat -B -e cache-references,cache-misses,L1-dcache-load-misses,L1-dcache-store-misses,L1-dcache-prefetch-misses,LLC-loads,LLC-stores,dTLB-load-misses,dTLB-store-misses,cycles -o $result_perf --append taskset -c 0-$k ./pagerank_pthreads $graphTxt $n 0.00000000001 0.85 $c $scd 1 0
					done
				fi
			done
		fi
	done
done
