import os
import csv
import sys

def main():
    if len(sys.argv) < 2:
        print("Usage: python3 update_mapping.py <outdir>")
        sys.exit(1)
        
    outdir = os.path.abspath(sys.argv[1])
    csv_path = os.path.join(outdir, "mapeo_muestras_parametros.csv")
    
    print(f"Actualizando mapeo en: {csv_path}")
    
    if not os.path.exists(outdir):
        print(f"Error: El directorio de salida {outdir} no existe.")
        sys.exit(1)

    csv_data = []

    # Escanear el directorio outdir
    for item in os.listdir(outdir):
        item_path = os.path.join(outdir, item)
        
        # Omitir carpetas grupales y archivos
        if item in ['summary', 'multiqc', 'allsorts', 'summary_low_mapping', 'multiqc_low_mapping']:
            continue
        if not os.path.isdir(item_path) and not os.path.islink(item_path):
            continue
            
        sample_id = item
        
        # 1. Buscar resultados estándar en outdir/<sample_id>/
        try:
            files_in_sample = os.listdir(item_path)
        except Exception as e:
            print(f"Advertencia: No se pudo leer {item_path}: {e}")
            continue
            
        if 'arriba_fusions.tsv' in files_in_sample:
            csv_data.append({
                'Sample_ID': sample_id,
                'Run_Type': 'standard',
                'Caller': 'arriba',
                'Filter_rRNA': 'false',
                'Parameters': '--method arriba --run_allsorts true',
                'Result_File': os.path.join(item_path, 'arriba_fusions.tsv')
            })
            
        if 'star-fusion.fusion_predictions.tsv' in files_in_sample:
            csv_data.append({
                'Sample_ID': sample_id,
                'Run_Type': 'standard',
                'Caller': 'starfusion',
                'Filter_rRNA': 'false',
                'Parameters': '--method starfusion --skip_fastqc true --run_allsorts false --annotate_common_fusions false',
                'Result_File': os.path.join(item_path, 'star-fusion.fusion_predictions.tsv')
            })
            
        if 'fusioncatcher_fusions.txt' in files_in_sample:
            csv_data.append({
                'Sample_ID': sample_id,
                'Run_Type': 'standard',
                'Caller': 'fusioncatcher',
                'Filter_rRNA': 'false',
                'Parameters': '--method fusioncatcher --fusioncatcher_dir /datos/migccl/fusioncatcher_data --skip_fastqc true --run_allsorts false --annotate_common_fusions false',
                'Result_File': os.path.join(item_path, 'fusioncatcher_fusions.txt')
            })
            
        # 2. Buscar resultados de bajo mapeo en outdir/<sample_id>/low_mapping_rrna_filter/
        low_mapping_dir = os.path.join(item_path, 'low_mapping_rrna_filter')
        if os.path.exists(low_mapping_dir) and (os.path.isdir(low_mapping_dir) or os.path.islink(low_mapping_dir)):
            try:
                files_in_low = os.listdir(low_mapping_dir)
            except Exception as e:
                print(f"Advertencia: No se pudo leer {low_mapping_dir}: {e}")
                files_in_low = []
                
            if 'arriba_fusions.tsv' in files_in_low:
                csv_data.append({
                    'Sample_ID': sample_id,
                    'Run_Type': 'low_mapping_rrna_filter',
                    'Caller': 'arriba',
                    'Filter_rRNA': 'true',
                    'Parameters': '--method arriba_starfusion --filter_rrna true --rrna_ref /datos/migccl/rRNA_refs/human_rRNA.fasta',
                    'Result_File': os.path.join(low_mapping_dir, 'arriba_fusions.tsv')
                })
                
            if 'star-fusion.fusion_predictions.tsv' in files_in_low:
                csv_data.append({
                    'Sample_ID': sample_id,
                    'Run_Type': 'low_mapping_rrna_filter',
                    'Caller': 'starfusion',
                    'Filter_rRNA': 'true',
                    'Parameters': '--method arriba_starfusion --filter_rrna true --rrna_ref /datos/migccl/rRNA_refs/human_rRNA.fasta',
                    'Result_File': os.path.join(low_mapping_dir, 'star-fusion.fusion_predictions.tsv')
                })

            if 'fusioncatcher_fusions.txt' in files_in_low:
                csv_data.append({
                    'Sample_ID': sample_id,
                    'Run_Type': 'low_mapping_rrna_filter',
                    'Caller': 'fusioncatcher',
                    'Filter_rRNA': 'true',
                    'Parameters': '--method fusioncatcher --fusioncatcher_dir /datos/migccl/fusioncatcher_data --filter_rrna true --skip_fastqc true',
                    'Result_File': os.path.join(low_mapping_dir, 'fusioncatcher_fusions.txt')
                })

    # Escribir CSV
    if csv_data:
        headers = ['Sample_ID', 'Run_Type', 'Caller', 'Filter_rRNA', 'Parameters', 'Result_File']
        with open(csv_path, 'w', newline='', encoding='utf-8') as f:
            writer = csv.DictWriter(f, fieldnames=headers)
            writer.writeheader()
            for row in csv_data:
                writer.writerow(row)
        print("Mapeo actualizado exitosamente.")
    else:
        print("No se encontraron resultados para mapear.")

if __name__ == "__main__":
    main()
