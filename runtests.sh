echo "" > $1
for i in $(ls $2)
do
    echo $i
    echo $i >> $1
    python main.py $2/$i >> $1
done
# uncomment to print to file $3 in format:
# fileName sequenceLength \n
# cat $1 | awk '$1 ~ /^[0-9]+$/ {print nazwa_pliku, $3} $1 ~ /txt$/ {nazwa_pliku = $1}' > $3