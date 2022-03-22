for thres in 2 4 6 8 10
do
	g++ -std=c++11 -O2 -DHOT_THRES=$thres main_SS.cpp -o SF -march=native
	
	if [ $? -ne 0 ]; then
	    echo "failed"
	else
	    echo "succeed"
	    ./SF
	fi
done
