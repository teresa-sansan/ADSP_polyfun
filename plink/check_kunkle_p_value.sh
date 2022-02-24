cd ~/output/cT/kunkle/

cat kunkle.pvalue | grep 'A\|T\|C\|G'|wc     ##2287
sed '/A\|T\|C\|G/d' kunkle.pvalue >kunkle_check.pvalue
cat kunkle_check.pvalue | grep 'A\|T\|C\|G'|wc
