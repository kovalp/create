rm timecreate
echo ' Begin running create.com' >> timecreate
date >> timecreate
cat << EOF | ./create.x > output.log
pavel             
0
-99
EOF
echo ' Ending running create.com' >> timecreate
date >> timecreate
