RESULTFILE=""
NUM=0
SAVE=0
VERBOSE=0

for i in $@
do
    NUM=$(($NUM+1))
    if [ $i == '-f' ]
    then
        SAVE=$(($NUM+1))
    fi
    if [ $NUM == $SAVE ]
    then
        RESULTFILE=$i
    fi
    if [ $i == '-v' ]
    then
        VERBOSE=1
        SAVE=$(($NUM))
    fi
done

if (($VERBOSE == 1))
then
    echo "Running tests for $(($NUM - $SAVE)) files"
fi

NUM=0
echo "" > $RESULTFILE

for i in $@
do
    NUM=$(($NUM+1))
    if (($NUM > $SAVE))
    then
        if (($VERBOSE == 1))
        then
            echo $i
        fi
        python main.py $i >> $RESULTFILE
    fi
done

if (($VERBOSE == 1))
then
    echo "Results saved to $RESULTFILE"
fi
