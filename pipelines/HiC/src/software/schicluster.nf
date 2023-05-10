nextflow.enable.dsl=2

process dim_reduct_schicluster {
    container 'kaizhang/schicluster:1.3.2'
    tag "$name"
    cpus 8

    input:
      tuple val(name), val(_), path("data.txt"), path("config.JSON"), path("chrom.sizes"), path("label_info.pickle")
    output:
      tuple val(name), val('scHiCluster'), path("reduced_dim.tsv")

    """
    #!/usr/bin/env python3
    import os
    import subprocess
    from multiprocessing import Pool
    import tempfile
    import json
    import h5py
    import numpy as np

    with open('config.JSON', 'r') as file:
        config = json.load(file)
        resolution = str(config['resolution'])
        chromosomes = config['chrom_list']

    def split_cells(input, out_dir):
        cell_list = []
        os.makedirs(out_dir, exist_ok=True)
        current_cell = None
        buffer = []
        with open(input, 'r') as fl:
            header = fl.readline().strip().split()
            cell_idx = header.index('cell_id')
            chrom1_idx = header.index('chrom1')
            chrom2_idx = header.index('chrom2')
            pos1_idx = header.index('pos1')
            pos2_idx = header.index('pos2')
            count_idx = header.index('count')
            for line in fl:
                items = line.strip().split()
                chr1 = items[chrom1_idx]
                chr2 = items[chrom2_idx]
                pos1 = items[pos1_idx]
                pos2 = items[pos2_idx]
                count = int(items[count_idx])
                cell = items[cell_idx]
                record = f"{chr1}\\t{chr2}\\t{pos1}\\t{pos2}"
                if current_cell is not None and current_cell != cell:
                    with open(out_dir + current_cell, 'w') as out:
                        out.write('\\n'.join(buffer))
                    cell_list.append(current_cell)
                    buffer = []
                current_cell = cell
                for _ in range(count):
                    buffer.append(record)

        if len(buffer) > 0:
            with open(out_dir + current_cell, 'w') as out:
                out.write('\\n'.join(buffer))
            cell_list.append(current_cell)
        return cell_list

    def generate_matrix(dir, cell, out_dir):
        os.makedirs(out_dir, exist_ok=True)
        subprocess.run(["hicluster", "generatematrix-cell",
            "--infile", dir + cell,
            "--outdir", out_dir + '/',
            "--chrom_file", "chrom.sizes",
            "--res", resolution,
            "--cell", cell, "--chr1", "0", "--pos1", "2", "--chr2", "1", "--pos2", "3"
        ], check = True)

    def impute_matrix(dir, cells, chromosome, out_dir):
        input = dir + '/' + chromosome + '/'
        out = out_dir + '/' + chromosome + '/'
        if len(os.listdir(input)) > 30:
            os.makedirs(out, exist_ok=True)
            for cell in cells:
                file = input + cell + '_' + chromosome + '.txt'
                open(file, 'a').close()
                subprocess.run(["hicluster", "impute-cell",
                    "--indir", input,
                    "--outdir", out,
                    "--cell", cell,
                    "--chrom", chromosome,
                    "--res", resolution,
                    "--chrom_file", "chrom.sizes"
                ], check = True)
            with open(out + "cell_list.txt", 'w') as f:
                f.write('\\n'.join([f"{out}/{cell}_{chromosome}_pad1_std1_rp0.5_sqrtvc.hdf5" for cell in cells]))
            subprocess.run(["hicluster", "embed-concatcell-chr",
                "--cell_list", out + "cell_list.txt",
                "--outprefix", out_dir + "/impute_" + chromosome,
                "--res", resolution,
                "--dim", "30",
            ], check = True)

    with tempfile.TemporaryDirectory(dir='./') as temp_dir:
        contact_map_dir = temp_dir + "/contact_map/"
        cells = split_cells("data.txt", contact_map_dir)

        matrix_dir = temp_dir + "/matrix/"
        with Pool(processes=8) as pool:
            pool.starmap(generate_matrix, [(contact_map_dir, cell, matrix_dir) for cell in cells])

        imputed_matrix_dir = temp_dir + "/imputed_matrix/"
        with Pool(processes=8) as pool:
            pool.starmap(impute_matrix, [(matrix_dir, cells, chromosome, imputed_matrix_dir) for chromosome in chromosomes])

        subprocess.run(f"ls {imputed_matrix_dir + '*npy'} > file_list.txt", shell=True, check = True)
        subprocess.run(["hicluster", "embed-mergechr",
            "--embed_list", "file_list.txt",
            "--outprefix", temp_dir + "/reduced_dim",
            "--dim", "30",
        ], check = True)

        with h5py.File(temp_dir + '/reduced_dim.svd30.hdf5', 'r') as f:
            np.savetxt("reduced_dim.tsv", f['data'][()], delimiter="\t")
    """
}