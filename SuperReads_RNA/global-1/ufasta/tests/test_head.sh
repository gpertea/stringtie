
x=$(ufasta head $TEST1 | grep -c '^>')
EXPECT_EQ 10 $x "Output default number of entries"

x=$(ufasta head -n 2 $TEST1 | grep -c '^>')
EXPECT_EQ 2 $x "Output 2 entries"
x=$(ufasta head -n -2 $TEST1 | grep -c '^>')
EXPECT_EQ 48 $x "Output 48 entries"

x=$(ufasta head -n 3 $TEST1 | wc -l)
r=$(head -n $((x + 1)) $TEST1 | tail -n 1)
EXPECT_EQ ">read3" "$r" "Output correct number of lines 1"
x2=$(ufasta head -n -3 $TEST1 | wc -l)
r2=$(head -n $x2 $TEST1 | grep '^>' | tail -n 1)
EXPECT_EQ ">read46" "$r2" "Output correct number of lines 2"

l=$(ufasta head -n 3 $TEST1)
el=$(head -n $x $TEST1)
EXPECT_EQ "$el" "$l" "Output correct lines 1"
l=$(ufasta head -n -3 $TEST1)
el=$(head -n $x2 $TEST1)
EXPECT_EQ "$el" "$l" "Output correct lines 2"

l=$(ufasta head -v $TEST1 | head -n 1)
EXPECT_EQ "==> $TEST1 <==" "$l" "Output header verbose"
l=$(ufasta head $TEST1 $TEST1 | head -n 1)
EXPECT_EQ "==> $TEST1 <==" "$l" "Output header many"
l=$(ufasta head -q $TEST1 $TEST1 | head -n 1)
EXPECT_EQ ">read0" "$l" "Output header quiet"


# Test byte switch
l=$(ufasta head -c 0 $TEST1)
EXPECT_EQ "" "$l" "-c 0"
for i in $(seq 1 50 $(wc -c < $TEST1)); do
    l=$(ufasta head -c $i $TEST1 | grep -c '^>')
    EXPECT_LT 0 "$l" "-c $i at least one entry"
    l=$(ufasta head -c $i $TEST1 | wc -c)
    EXPECT_LE $i $l "-c $i at least $i bytes"
    
    l=$(ufasta head -c $i $TEST1; ufasta tail -c +$i $TEST1)
    EXPECT_EQ "$(cat $TEST1)" "$l" "tail + head -c +$i"
done
