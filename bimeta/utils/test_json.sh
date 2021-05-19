JSON_PATH=/home/dhuy237/thesis/code/bimetaReduce/bimeta/data/testMRS/21_05_19_09_42_33.json

value=($(jq -r '.steps.Step3' $JSON_PATH))

# a="${value}"

echo $value

STEP_1_3=($(jq -r '.steps.Step3' $JSON_PATH))

if [ "$STEP_1_3" = "true" ]
then
    echo "1"
else
    echo "2"
fi

NUM_OF_SPECIES=($(jq -r '.params.kNumber' $JSON_PATH))

if [ "$NUM_OF_SPECIES" = "false" ]
then
    echo "1"
else
    echo "2"
fi

DATA_PATH=bimeta/data/testMRS
OUTPUT_GRAPH=($(jq -r '.nodeGraph' $JSON_PATH))

echo $DATA_PATH/$OUTPUT_GRAPH