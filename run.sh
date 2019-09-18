g++ *.cpp -o a.out

for i in $(seq 1 $1); do
    mkdir -p data$i
done

./a.out

