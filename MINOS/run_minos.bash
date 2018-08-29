#!/bin/bash
#
# Run a command
#
file=$1

python write_card_file.py << EOF2
$1
EOF2

./minos_bran << EOF
$file.txt
out.txt
$file.tmp
1.d-8 100
1
1 30 1. 10. 0 30
EOF
cat out.txt | grep " s " | awk {'print $1, $2, $3, $5,  $6'} > $file.modes

#Clean up
#rm out.txt
rm $file.tmp

exit
