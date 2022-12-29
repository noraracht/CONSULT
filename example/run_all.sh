echo "EXAMPLE: RUN ALL"

echo -e "\n%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n"
echo "===COMPILING==="
echo -e "\n%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n"
make clean -C ..
make all -C ..

echo -e "\n%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n"
echo "===MINIMIZATION==="
echo -e "\n%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n"
../minimize \
  -i "k35C_bef_mininimization.fa" \
  -o "k32C_af_mininimization.fa"

echo -e "\n%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n"
echo "===MAPPING==="
echo -e "\n%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n"
[ ! -e G000307305_nbr_mapping ] || rm -Ir G000307305_nbr_mapping
../map \
  -i "k32C_af_mininimization.fa" \
  -o "G000307305_nbr_mapping" \
  -p "3"

echo -e "\n%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n"
echo "===SEARCHING==="
echo -e "\n%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n"
../search \
  -q "G000307305.fq" \
  -i "G000307305_nbr_mapping" \
  -o "." \
  -c 1 -t 1 \
  --classify-reads --save-distances
