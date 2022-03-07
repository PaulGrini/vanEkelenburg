# Collect all counts.tsv that were renamed accordingly

for D in map_*/;
do
    cp ../$D/count_*.tsv .
done
echo DONE
