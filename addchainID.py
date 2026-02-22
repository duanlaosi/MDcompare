def replace_chain_id_with_A(pdb_file, output_file):
    with open(pdb_file, 'r') as f:
        lines = f.readlines()

    with open(output_file, 'w') as f:
        for line in lines:
    
            if line.startswith(('ATOM', 'TER')):
                
                if line[21] == ' ':
                    
                    line = line[:21] + 'A' + line[22:]
            f.write(line)

#usage
replace_chain_id_with_A('/home/duanqi/dq/1BPI/1BPI_pdb_4_23/1BPI_pdb_str/0.15_80_10_pH6.5_1bpi_renumber_refine_noh.result.pdb', '1bpi_renumber_refine_noh_protonation_chainIDA.pdb')