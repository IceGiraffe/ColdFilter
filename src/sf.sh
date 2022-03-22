rm SF.out
g++ SF.cpp -o SF.out -march=native

if [ $? -ne 0 ]; then
    echo "failed"
else
    echo "succeed"
    ./SF.out
fi