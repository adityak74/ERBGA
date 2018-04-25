#1/bin/bash
make #run make for compiling bfs
for i in {1..5}
do
	./erbga $1 $1.out | tee $1-console > /dev/null &
	sleep 1
done
