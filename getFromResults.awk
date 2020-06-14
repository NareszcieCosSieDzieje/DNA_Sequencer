/^OptLen/ {total += $9; totalS += $9/$3; num += 1}
END {print total, totalS, num, totalS/num}