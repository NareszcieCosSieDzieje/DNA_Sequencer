echo "" > $1
for i in $(ls $2)
do
    echo $i
    echo $i >> $1
    python main.py $2/$i >> $1
done