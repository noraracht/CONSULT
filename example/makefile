all: clean compile minimize map search init-ID update-ID save-matches classify

clean:
	@make clean -C ..
	[ ! -e G000307305-k32C_minimized.fa   ] || rm G000307305-k32C_minimized.fa
	[ ! -e G000307305_nbr_mapping         ] || rm -Ir G000307305_nbr_mapping
	[ ! -e classified-seq_query000        ] || rm classified-seq_query000
	[ ! -e unclassified-seq_query000      ] || rm unclassified-seq_query000
	[ ! -e kmer-distances_query000        ] || rm kmer-distances_query000
	[ ! -e match-info_query000            ] || rm match-info_query000
	[ ! -e predictions_query000           ] || rm predictions_query000

compile:
	@echo -e "\n%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n"
	@echo -e "\n===COMPILING===\n"
	@echo -e "\n%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n"
	@make minimize -C ..
	@make map -C ..
	@make search -C ..
	@make classify -C ..

minimize:
	@echo -e "\n%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n"
	@echo "===MINIMIZATION==="
	@echo -e "\n%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n"
	@../minimize \
		-i "G000307305-k35C.fa" \
		-o "G000307305-k32C_minimized.fa"

search:
	@echo -e "\n%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n"
	@echo "===SEARCHING==="
	@echo -e "\n%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n"
	@../consult_search \
		-q "query000.fq" \
		-i "G000307305_nbr_mapping" \
		-o "." \
		-c 1 \
		--thread-count 12 \
		--unclassified-out --classified-out --save-distances

init-ID:
	@echo -e "\n%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n"
	@echo "===INITIALIZING-IDs==="
	@echo -e "\n%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n"
	@../consult_search \
		-q "G000307305-k32C_minimized.fq" \
		-i "G000307305_nbr_mapping" \
		-o "." \
		--thread-count 12 \
		--taxonomy-lookup-path "taxonomy-lookup" \
		--filename-map-path "filename-map" \
		--init-ID

update-ID:
	@echo -e "\n%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n"
	@echo "===UPDATING-IDs==="
	@echo -e "\n%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n"
	@../consult_search \
		-q "G000307305-k32C_minimized.fq" \
		-i "G000307305_nbr_mapping" \
		-o "." \
		--thread-count 12 \
		--taxonomy-lookup-path "taxonomy-lookup" \
		--filename-map-path "filename-map" \
		--update-ID

save-matches:
	@echo -e "\n%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n"
	@echo "===SAVING-MATCHES==="
	@echo -e "\n%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n"
	@../consult_search \
		-q "query000.fq" \
		-i "G000307305_nbr_mapping" \
		-o "." \
		--thread-count 12 \
		--taxonomy-lookup-path "taxonomy-lookup" \
		--filename-map-path "filename-map" \
		--save-matches

classify:
	@echo -e "\n%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n"
	@echo "===CLASSIFICATION==="
	@echo -e "\n%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n"
	@../consult_classify \
		-i "match-info_query000" \
		-o "." \
		--taxonomy-lookup-path "taxonomy-lookup" \
		--thread-count 12

map:
	@echo -e "\n===COMPILING===\n"
	@echo -e "\n%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n"
	@echo "===MAPPING==="
	@echo -e "\n%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n"
	@[ ! -e G000307305_nbr_mapping ] || rm -Ir G000307305_nbr_mapping
	@../consult_map \
		-i "G000307305-k32C_minimized.fa" \
		-o "G000307305_nbr_mapping" \
		-p 3 -t 2
