# Collect counts.tsv files in this folder

for D in map_*/;
do
    cp ../$D/*.tsv .
done
echo DONE
