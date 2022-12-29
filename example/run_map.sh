echo "EXAMPLE: RUN MAP"

echo -e "\n%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n"
echo "===COMPILING==="
echo -e "\n%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n"
make clean -C
make map -C

echo -e "\n%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n"
echo "===MAPPING==="
echo -e "\n%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n"
[ ! -e G000307305_nbr_mapping ] || rm -Ir G000307305_nbr_mapping
../map \
  -i "k32C_af_mininimization.fa" \
  -o "G000307305_nbr_mapping" \
  -p 3
