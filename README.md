Pipeline de Análise de Variantes Somáticas
=======================================================

1. Restaurar ambiente R (renv):
   Rscript -e 'renv::restore()'

2. Estrutura de pastas:
   - samples/ : arquivos de variantes em .txt e groups.txt
   - results/ : saídas (plots, rds, tables)
   - snvSignatures.R : script para análise de variantes de nucleotídeo único (SNVs)
   - indelSignatures.R : script para análise de variantes de inserção ou deleção (InDel)
   - config.yaml : configurações de caminhos

3. Configuração (config.yaml):
   vcf_dir: "samples"
   group_file: "groups.txt"

4. Resultados:
   - results/plots : gráficos gerados com o MutationalPatterns
   - results/rds   : objetos R intermediários
   - results/tables: resumos e conferências de contagem

