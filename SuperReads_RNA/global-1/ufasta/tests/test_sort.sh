
SORT_SIZES=$(ufasta sort $TEST1 | ufasta sizes -H /dev/fd/0)
SIZES_SORT=$(ufasta sizes -H $TEST1 | sort -k1,1)
EXPECT_EQ "$SIZES_SORT" "$SORT_SIZES" "Sorted sizes"

SORT_SIZES=$(sed -e 's/>read\([0-9]\+\)/>blah read\1/' $TEST1 | ufasta sort -k 2 /dev/stdin | ufasta sizes -H /dev/fd/0)
SIZES_SORT=$(ufasta sizes -H $TEST1 | sort -k1,1 | sed -e 's/read[0-9]\+/blah/')
EXPECT_EQ "$SIZES_SORT" "$SORT_SIZES" "Sorted sizes, second column"

SORT_SIZES=$(sed -e 's/>read\([0-9]\+\)/>blah read\1/' $TEST1 | ufasta sort -C 5 -k 2 -n /dev/stdin | ufasta sizes -H /dev/fd/0)
SIZES_SORT=$(ufasta sizes -H $TEST1 | sort -n -k1.5,1.10 | sed -e 's/read[0-9]\+/blah/')
EXPECT_EQ "$SIZES_SORT" "$SORT_SIZES" "Sorted sizes, numeric, second column"


SORT_SORT=$(ufasta sort $TEST1 | sort)
SORT_ONCE=$(sort $TEST1)
EXPECT_EQ "$SORT_ONCE" "$SORT_SORT" "Sorted content"

SORT_SIZES=$(ufasta sort -C 5 -n $TEST1 | ufasta sizes -H /dev/fd/0)
SIZES_SORT=$(ufasta sizes -H $TEST1 | sort -n -k1.5,1.10)
EXPECT_EQ "$SIZES_SORT" "$SORT_SIZES" "Sorted numeric sizes"

SORT_SORT=$(cat $TEST1 | ufasta sort -C 5 -n /dev/fd/0 | sort)
SORT_ONCE=$(sort $TEST1)
EXPECT_EQ "$SORT_ONCE" "$SORT_SORT" "Sorted numeric content"

SORT_SORT=$(echo -e ">1 2\n>1 1" | ufasta sort /dev/fd/0)
EXPECT_EQ "$(echo -e ">1 2\n>1 1")" "$SORT_SORT" "Sort first word"
SORT_SORT=$(echo -e ">1 2\n>1 1" | ufasta sort -H /dev/fd/0)
EXPECT_EQ "$(echo -e ">1 1\n>1 2")" "$SORT_SORT" "Sort full header line"

