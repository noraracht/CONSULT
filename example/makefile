all: clean minimize map search

clean:
	@make clean -C ..
	[ ! -e k32C_af_mininimization.fa   ] || rm k32C_af_mininimization.fa
	[ ! -e G000307305_nbr_mapping      ] || rm -r G000307305_nbr_mapping
	[ ! -e classified-seq_G000307305   ] || rm classified-seq_G000307305
	[ ! -e unclassified-seq_G000307305 ] || rm unclassified-seq_G000307305
	[ ! -e kmer-distances_G000307305   ] || rm kmer-distances_G000307305

minimize:
	@echo -e "\n===COMPILING===\n"
	@make minimize -C ..
	@echo -e "\n%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n"
	@echo "===MINIMIZATION==="
	@echo -e "\n%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n"
	@../minimize \
		-i "k35C_bef_mininimization.fa" \
		-o "k32C_af_mininimization.fa"

search:
	@echo -e "\n===COMPILING===\n"
	@make search -C ..
	@echo -e "\n%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n"
	@echo "===SEARCHING==="
	@echo -e "\n%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n"
	@../consult_search \
		-q "G000307305.fq" \
		-i "G000307305_nbr_mapping" \
		-o "." \
		-c 1 \
		--thread-count 1 \
		--unclassified-out --classified-out \
		--save-distances --maximum-distance 5

map:
	@echo -e "\n===COMPILING===\n"
	@make map -C ..
	@echo -e "\n%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n"
	@echo "===MAPPING==="
	@echo -e "\n%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n"
	@[ ! -e G000307305_nbr_mapping ] || rm -Ir G000307305_nbr_mapping
	@../consult_map \
		-i "k32C_af_mininimization.fa" \
		-o "G000307305_nbr_mapping" \
		-p 3 \
		-t 2
