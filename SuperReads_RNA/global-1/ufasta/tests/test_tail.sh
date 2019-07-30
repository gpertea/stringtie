
x=$(ufasta tail -n +1 $TEST1 | grep -c '^>')
EXPECT_EQ "50" $x "Output all 1"
x=$(ufasta tail -n 50 $TEST1 | grep -c '^>')
EXPECT_EQ "50" $x "Output all 2"

x=$(ufasta tail -n +2 $TEST1 | grep -c '^>')
EXPECT_EQ "49" $x "Output all but first 1"
x=$(ufasta tail -n 49 $TEST1 | grep -c '^>')
EXPECT_EQ "49" $x "Output all but first 2"

x=$(ufasta tail -n +49 $TEST1 | wc -l)
eh=$(tail -n $x $TEST1 | head -n 1)
EXPECT_EQ ">read48" "$eh" "Output last two headers"
x2=$(ufasta tail -n 2 $TEST1 | wc -l)
EXPECT_EQ $x $x2 "Complementary tail"

el=$(tail -n $x $TEST1)
l=$(ufasta tail -n +49 $TEST1)
EXPECT_EQ "$el" "$l" "Output last two entries 1"
l=$(ufasta tail -n 2 $TEST1)
EXPECT_EQ "$el" "$l" "Output last two entries 2"

l=$(ufasta tail -c +0 $TEST1)
EXPECT_EQ "$(cat $TEST1)" "$l" "-c +0"

test1_len=$(wc -c < $TEST1)
for i in $(seq 1 53 $(wc -c < $TEST1)); do
    l=$(ufasta tail -c +$i $TEST1 | wc -c)
    EXPECT_GE $((test1_len - i)) $l "-c +$i"

    l=$(ufasta head -c $i $TEST1; ufasta tail -c +$i $TEST1)
    EXPECT_EQ "$(cat $TEST1)" "$l" "tail + head -c +$i"
done
