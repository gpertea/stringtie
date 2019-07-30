
trap "rm -f split1 split2" EXIT
ufasta split -i $TEST1 >(cat > split1) >(cat > split2)
diff <(cat split1 split2 | ufasta sort) <(ufasta sort $TEST1)
EXPECT_SUCCESS "Split join"
