echo "EXAMPLE: RUN SEARCH"

echo -e "\n%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n"
echo "===COMPILING==="
echo -e "\n%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n"
make clean -C ..
make search -C ..

echo -e "\n%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n"
echo "===SEARCHING==="
echo -e "\n%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n"
../search \
  -q "G000307305.fq" \
  -i "G000307305_nbr_mapping" \
  -o "." \
  -c 1 \
  --thread-count 1 \
  --classify-reads --save-distances
